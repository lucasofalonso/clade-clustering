import re
from pathlib import Path

import pandas as pd


NUM_METADATA_COLUMNS = 9

ACCESSION_PATTERN = re.compile(r"GC[AF]_\d+\.\d+")

# Artifact file names inside a matrices directory
X_GENOTYPES_FILE = "X_genotypes.parquet"
SNP_INFO_FILE = "snp_info.parquet"
SAMPLE_METADATA_FILE = "sample_metadata.parquet"


def read_vcf(path: Path) -> pd.DataFrame:
    """
    Read a VCF file into a DataFrame with one row per variant.

    The "##" meta-information lines are skipped and the "#CHROM" header line
    provides the column names. The 9 fixed VCF columns are normalized to
    lowercase without the leading "#"; sample names are kept as they are. No
    variant is filtered and no genotype is parsed.

    Args:
        path: VCF file written by HarvestTools.

    Returns:
        DataFrame with the 9 fixed VCF columns (chrom, pos, id, ref, alt, qual,
        filter, info, format) followed by one column per sample, in file order.

    Raises:
        ValueError: If the file has no "#CHROM" header line.
    """
    header = None
    header_idx = None

    with open(path) as f:
        for i, line in enumerate(f):
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                header_idx = i
                break

    if header is None:
        raise ValueError(f"No #CHROM header line found in VCF: {path}")

    header[:NUM_METADATA_COLUMNS] = [col.lstrip("#").lower() for col in header[:NUM_METADATA_COLUMNS]]

    return pd.read_csv(
        path,
        sep="\t",
        skiprows=header_idx + 1,
        names=header,
        # Text fields stay text even if every value looks numeric (e.g. chrom "1")
        dtype={"chrom": str, "id": str, "ref": str, "alt": str, "filter": str},
        low_memory=False,
    )


def _has_single_base_alleles(alt: str) -> bool:
    """
    Check whether every allele in a VCF ALT field is a single base.

    Args:
        alt: ALT field of a VCF record, with one or more alleles separated by
            commas (e.g. "T" or "T,C").

    Returns:
        True if every allele has length 1 ("T", "T,C"); False if any allele
        is longer ("TA", "T,CG").
    """
    return all(len(allele) == 1 for allele in alt.split(","))


def filter_pass_snps(df_vcf: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only PASS sites whose REF and every ALT allele are single bases.

    Both biallelic ("A" > "G") and multiallelic ("A" > "T,C") sites are kept;
    any site with an indel allele is dropped. FILTER is compared
    case-insensitively.

    Args:
        df_vcf: Variants as returned by read_vcf.

    Returns:
        The selected rows, in file order, with a new 0..n-1 index.
    """
    is_pass = df_vcf["filter"].str.lower().eq("pass")

    is_snp = (
        df_vcf["ref"].str.len().eq(1)
        & df_vcf["alt"].apply(_has_single_base_alleles)
    )

    df_snps = df_vcf[is_pass & is_snp].reset_index(drop=True)

    return df_snps


def build_snp_info(df_snps: pd.DataFrame) -> pd.DataFrame:
    """
    Build the locus table that describes each column of the genotype matrix.

    Each locus gets a snp_id in the form "chrom:pos:ref>alt" (e.g.
    "CP043531.1:60467:A>G"). Alleles are kept exactly as written in the VCF,
    including the case of REF and the comma-separated ALT of multiallelic
    sites ("CP043531.1:60452:c>C,T").
    
    Args:
        df_snps: Selected sites as returned by filter_pass_snps.

    Returns:
        DataFrame with one row per locus, in the same order as df_snps, and
        the columns chrom, pos, id, ref, alt, qual, filter and snp_id.
    
    Raises:
        ValueError: If two sites produce the same snp_id.
    """
    df_snp_info = df_snps[["chrom", "pos", "id", "ref", "alt", "qual", "filter"]]

    df_snp_info["snp_id"] = (
        df_snp_info["chrom"]
        + ":"
        + df_snp_info["pos"].astype(str)
        + ":"
        + df_snp_info["ref"]
        + ">"
        + df_snp_info["alt"]
    )

    if df_snp_info["snp_id"].duplicated().any():
        raise ValueError("Duplicated snp_id found in the selected sites.")

    return df_snp_info


def build_genotype_matrix(df_snps: pd.DataFrame, df_snp_info: pd.DataFrame) -> pd.DataFrame:
    """
    Build the sample x locus genotype matrix from the selected sites.

    The VCF layout (one row per locus, one column per sample) is transposed so
    that each row is a sample and each column a locus named by its snp_id.
    Each value is the integer allele index called for that sample: 0 for REF,
    1 for the first ALT, 2 for the second, and so on. Every genotype must be
    called: missing calls (".") or any non-integer value are rejected rather
    than converted to NaN.

    Args:
        df_snps: Selected sites as returned by filter_pass_snps.
        df_snp_info: Locus table built from the same df_snps by build_snp_info;
            its rows are matched to df_snps by position.

    Returns:
        DataFrame of int8 indexed by "vcf_sample" (sample names as in the
        VCF header, including the reference sample) with one column per
        snp_id, in the same order as df_snp_info.

    Raises:
        ValueError: If df_snp_info and df_snps have a different number of
            rows, if FORMAT has any field other than GT, or if any sample
            column holds a non-integer genotype.
    """
    if len(df_snp_info) != len(df_snps):
        raise ValueError(
            f"Got {len(df_snp_info)} snp_ids for {len(df_snps)} sites; "
            "df_snp_info must be built from the same df_snps."
        )

    if not df_snps["format"].eq("GT").all():
        raise ValueError("Expected FORMAT = GT only.")

    df_samples_gt = df_snps.iloc[:, NUM_METADATA_COLUMNS:]

    # Any missing call (".") or extra field makes read_csv parse the column as text
    non_integer = [
        sample for sample, dtype in df_samples_gt.dtypes.items()
        if not pd.api.types.is_integer_dtype(dtype)
    ]
    if non_integer:
        raise ValueError(
            f"{len(non_integer)} sample column(s) have non-integer genotypes, "
            f"e.g. {non_integer[:3]}"
        )

    # Allele indexes are small (0..number of ALT alleles), so int8 holds them
    df_genotypes = df_samples_gt.astype("int8").T
    df_genotypes.columns = df_snp_info["snp_id"].to_numpy()
    df_genotypes.index.name = "vcf_sample"

    return df_genotypes


def _extract_accession(sample_name: str) -> str | None:
    """
    Extract the NCBI assembly accession from a VCF sample name.

    Args:
        sample_name: Sample name as in the VCF header
            (e.g. "GCA_044589395.1_ASM4458939v1_genomic.fna").

    Returns:
        The first GCA_/GCF_ accession found (e.g. "GCA_044589395.1"), or None
        if the name has none (e.g. "Ca9.fasta").
    """
    match = ACCESSION_PATTERN.search(sample_name)
    return match.group(0) if match else None


def map_samples_to_accessions(vcf_samples: pd.Index) -> pd.DataFrame:
    """
    Map each VCF sample name to the identifier used across the pipeline.

    NCBI samples are identified by the assembly accession in their name.
    Samples without one, such as the manually added genomes ("Ca9.fasta",
    "SRR9007776.pilon.fasta"), keep their VCF sample name as identifier. The
    reference sample (name ending in ".ref") is dropped.

    Args:
        vcf_samples: Sample names as in the VCF header, e.g. the index of the
            matrix returned by build_genotype_matrix.

    Returns:
        DataFrame with the columns vcf_sample and accession, one row per
        sample except the reference, in the same order as vcf_samples and with
        a new 0..n-1 index.

    Raises:
        ValueError: If there is not exactly one reference sample, or if two
            samples end up with the same accession.
    """
    sample_accessions = pd.DataFrame({"vcf_sample": vcf_samples})
    sample_accessions["accession"] = sample_accessions["vcf_sample"].apply(_extract_accession)

    is_reference = sample_accessions["vcf_sample"].str.endswith(".ref")
    n_reference = int(is_reference.sum())
    if n_reference != 1:
        raise ValueError(
            f"Expected exactly one reference sample ending in '.ref', found {n_reference}."
        )
    sample_accessions = sample_accessions.loc[~is_reference].reset_index(drop=True)

    sample_accessions["accession"] = sample_accessions["accession"].fillna(
        sample_accessions["vcf_sample"]
    )

    is_duplicated = sample_accessions["accession"].duplicated(keep=False)
    if is_duplicated.any():
        duplicated = sample_accessions.loc[is_duplicated, "vcf_sample"].tolist()
        raise ValueError(f"Samples with the same accession: {duplicated}")

    return sample_accessions


def merge_sample_with_metadata(df_sample_accessions: pd.DataFrame, df_meta_clean: pd.DataFrame) -> pd.DataFrame:
    """
    Attach the cleaned NCBI metadata to each VCF sample by accession.

    Every sample is kept (left join). Samples absent from df_meta_clean, such
    as the manually added genomes, get missing values in the metadata columns.

    Args:
        df_sample_accessions: Sample table as returned by
            map_samples_to_accessions.
        df_meta_clean: Cleaned sample metadata with an "accession" column
            (sample_metadata_clean.parquet).

    Returns:
        DataFrame with the columns vcf_sample and accession followed by the
        metadata columns, one row per sample, in the same order as
        df_sample_accessions.

    Raises:
        pandas.errors.MergeError: If an accession appears more than once in
            df_meta_clean, which would duplicate samples.
    """
    return df_sample_accessions.merge(
        df_meta_clean,
        on="accession",
        how="left",
        validate="one_to_one"
    )


def align_by_accession(X_genotypes: pd.DataFrame, df_sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Index the genotype matrix and the sample metadata by accession, row by row.

    Rows of X_genotypes are selected and ordered by label to match the
    vcf_sample column of df_sample_metadata, so samples without metadata rows
    (the reference) are dropped. The vcf_sample index of the matrix is then
    replaced by the accession; the original sample name stays available in the
    vcf_sample column of the metadata.

    Args:
        X_genotypes: Matrix indexed by vcf_sample, as returned by
            build_genotype_matrix.
        df_sample_metadata: Sample table as returned by
            merge_sample_with_metadata.

    Returns:
        A tuple (X_genotypes, df_sample_metadata), both indexed by "accession"
        with the same rows in the same order.

    Raises:
        KeyError: If a vcf_sample of df_sample_metadata is missing from the
            index of X_genotypes.
    """
    df_sample_metadata = df_sample_metadata.set_index("accession")
    # Order X as df_sample_metadata; the .ref is dropped because it has no metadata row
    X_genotypes = X_genotypes.loc[df_sample_metadata["vcf_sample"]]
    X_genotypes.index = df_sample_metadata.index

    return X_genotypes, df_sample_metadata


def build_matrices_from_vcf(
    vcf_path: Path,
    sample_metadata_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build the genotype matrix, the locus table and the sample metadata from a VCF.

    Runs, in order: read_vcf, filter_pass_snps, build_snp_info,
    build_genotype_matrix, map_samples_to_accessions,
    merge_sample_with_metadata and align_by_accession. Nothing is written to
    disk; use save_genotype_matrices for that.

    Args:
        vcf_path: VCF file written by HarvestTools (e.g. parsnp.vcf).
        sample_metadata_path: Cleaned sample metadata
            (sample_metadata_clean.parquet).

    Returns:
        A tuple (X_genotypes, df_snp_info, df_sample_metadata), ready for
        save_genotype_matrices. The matrix and the metadata are indexed by
        accession.
    """
    df_vcf = read_vcf(vcf_path)
    df_snps = filter_pass_snps(df_vcf)
    df_snp_info = build_snp_info(df_snps)
    X_genotypes = build_genotype_matrix(df_snps, df_snp_info)

    df_sample_accessions = map_samples_to_accessions(X_genotypes.index)
    df_meta_clean = pd.read_parquet(sample_metadata_path)
    df_sample_metadata = merge_sample_with_metadata(df_sample_accessions, df_meta_clean)

    X_genotypes, df_sample_metadata = align_by_accession(X_genotypes, df_sample_metadata)

    return X_genotypes, df_snp_info, df_sample_metadata


def save_genotype_matrices(
    X_genotypes: pd.DataFrame,
    df_snp_info: pd.DataFrame,
    df_sample_metadata: pd.DataFrame,
    out_dir: Path,
) -> None:
    """
    Save the genotype matrix, the locus table and the sample metadata as parquet.

    The directory is created if needed and existing files are overwritten.

    Args:
        X_genotypes: Matrix indexed by accession, as returned by
            align_by_accession.
        df_snp_info: Locus table as returned by build_snp_info.
        df_sample_metadata: Metadata indexed by accession, as returned by
            align_by_accession.
        out_dir: Directory that receives X_genotypes.parquet, snp_info.parquet
            and sample_metadata.parquet.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    X_genotypes.to_parquet(out_dir / X_GENOTYPES_FILE)
    # The 0..n-1 index carries no information
    df_snp_info.to_parquet(out_dir / SNP_INFO_FILE, index=False)
    df_sample_metadata.to_parquet(out_dir / SAMPLE_METADATA_FILE)


def load_genotype_matrices(in_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load the artifacts written by save_genotype_matrices.

    Only the three known files are read; any other file in the directory
    (e.g. artifacts from later stages) is ignored.

    Args:
        in_dir: Directory with X_genotypes.parquet, snp_info.parquet and
            sample_metadata.parquet.

    Returns:
        A tuple (X_genotypes, df_snp_info, df_sample_metadata), in the same
        order as the arguments of save_genotype_matrices.

    Raises:
        FileNotFoundError: If any of the three files is missing.
    """
    X_genotypes = pd.read_parquet(in_dir / X_GENOTYPES_FILE)
    df_snp_info = pd.read_parquet(in_dir / SNP_INFO_FILE)
    df_sample_metadata = pd.read_parquet(in_dir / SAMPLE_METADATA_FILE)

    return X_genotypes, df_snp_info, df_sample_metadata
