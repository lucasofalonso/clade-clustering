import pandas as pd
import pytest

from clade_clustering.genotypes import filter_pass_snps, read_vcf


VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_a\tsample_b"

VCF_RECORDS = [
    "chr1\t100\t.\tA\tG\t40\tPASS\tNA\tGT\t0\t1",
    "chr1\t200\t.\tC\tT,A\t40\tLCB\tNA\tGT\t1\t2",
]


def write_vcf(path, header, records):
    """Write a minimal VCF with two meta-information lines."""
    lines = [
        "##fileformat=VCFv4.1",
        '##FILTER=<ID=LCB,Description="LCB smaller than 200bp">',
        header,
        *records,
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


# C1: read_vcf

def test_read_vcf_columns_and_values(tmp_path):
    df = read_vcf(write_vcf(tmp_path / "test.vcf", VCF_HEADER, VCF_RECORDS))

    assert df.shape == (2, 11)
    # Fixed columns are lowercased; sample names are kept as they are
    assert df.columns.tolist() == [
        "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format",
        "sample_a", "sample_b",
    ]
    assert df["pos"].tolist() == [100, 200]
    assert df["alt"].tolist() == ["G", "T,A"]
    assert df["sample_b"].tolist() == [1, 2]
    # "NA" in INFO is parsed as missing, as in the Parsnp VCF
    assert df["info"].isna().all()


def test_read_vcf_without_header_raises(tmp_path):
    path = write_vcf(tmp_path / "no_header.vcf", "##no header here", VCF_RECORDS)

    with pytest.raises(ValueError, match="#CHROM"):
        read_vcf(path)


# C2: filter_pass_snps

def test_filter_pass_snps_keeps_single_base_sites():
    df_vcf = pd.DataFrame({
        "pos":    [1,      2,      3,     4,      5,      6,      7],
        "ref":    ["A",    "A",    "A",   "AT",   "A",    "A",    "A"],
        "alt":    ["G",    "G",    "G",   "A",    "AT",   "T,C",  "T,CG"],
        "filter": ["PASS", "pass", "LCB", "PASS", "PASS", "PASS", "PASS"],
    })
    # 1: kept    2: kept (case-insensitive)    3: not PASS
    # 4: deletion    5: insertion    6: kept (multiallelic)
    # 7: multiallelic with an indel allele

    df_snps = filter_pass_snps(df_vcf)

    assert df_snps["pos"].tolist() == [1, 2, 6]
    assert df_snps.index.tolist() == [0, 1, 2]
    assert df_snps.columns.equals(df_vcf.columns)
