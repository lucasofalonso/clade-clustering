# clade-clustering

Clade identification for *Candida auris* (*Candidozyma auris*) from core-genome SNPs, using unsupervised learning and per-clade allelic signatures.

## Introduction

*C. auris* is a multidrug-resistant fungal pathogen split into phylogeographic clades (I to VI) with distinct epidemiological and resistance profiles. Knowing an isolate's clade supports genomic surveillance and outbreak investigation.

This project develops, as a series of Jupyter notebooks, the building blocks of a reproducible pipeline that:

1. downloads every public *C. auris* assembly from NCBI and applies basic quality control;
2. aligns the approved genomes, plus 11 manually added genomes shipped with the repository, against the B11220 reference (clade II) with **Parsnp**, and extracts SNPs with **HarvestTools**;
3. converts the VCF into genotype matrices (sample × locus);
4. recovers the clade structure without supervision (PCA + KMeans and Ward hierarchical clustering) and names each cluster using genomes whose clade was assigned by experts (*anchoring*);
5. derives per-clade **allelic signatures** (one-vs-rest allele frequency differences), validates them by cross-validation and on an independent external panel, and annotates the markers against the reference GFF to find coding genes that recur across each clade's signature.

### Project status

The repository is at the prototyping stage. The notebooks are not the final product: each one in `notebooks/stage/` is a building block of a future end-to-end pipeline and genomic surveillance CLI. Each block has a single responsibility, reads and writes explicit files on disk and checks its inputs and outputs, so it can later become a script or a pipeline step. Moving the notebook logic into `src/clade_clustering/` has started (see [Python package](#python-package)): the Parsnp matrix construction of notebook 03 is already available as functions, but the notebooks themselves still run their original code. Orchestration of the whole pipeline and the command-line interface are still to be done; for now, the steps are run by hand, in order.

Notebook 05 already produces the project's main results (allelic signatures, cross-validation, hold-out and external classification, coding gene recurrence), but it has not been refined yet. It is still long and exploratory, with inspection cells kept alongside the analysis and some code copied from notebook 04 instead of shared through `src/`. Its results are the ones reported below; its structure will change.

### Key results (snapshot `v4`)

| Step | Result |
|---|---|
| NCBI genomes after QC (manifest) | 1,244 of 1,247 (including the reference) |
| Manually added genomes (hold-out) | 11 (Ca9–Ca14 and 5 SRR assemblies) |
| SNPs (`FILTER = PASS`, simple SNPs) | 97,034 loci → 105,317 allele indicators after encoding multiallelic loci |
| Genomes used to fit PCA/KMeans | 1,243 |
| KMeans (k = 6, 4 PCs) | silhouette 0.998; ARI = 1.0 against Ward |
| Cluster → clade anchoring | 140 expert-labelled genomes; each cluster maps to exactly one clade |
| Signature CV (clades I–IV, 3 folds × 5 repeats) | accuracy 1.0 in every fold; mean margin 0.94 |
| External panel (30 genomes, 5 per clade) | accuracy 0.967 (29/30) |

Read these numbers together with the [Limitations](#limitations): clades V and VI have only 1 and 3 genomes in the fitting set.

## Table of contents

- [Introduction](#introduction)
  - [Project status](#project-status)
  - [Key results](#key-results-snapshot-v4)
- [Repository layout](#repository-layout)
- [Pipeline](#pipeline)
- [Python package](#python-package)
- [Reproducibility](#reproducibility)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Alignment tools (notebook 02)](#alignment-tools-notebook-02)
  - [Input data](#input-data)
  - [Reference genome annotation](#reference-genome-annotation)
  - [Running the pipeline](#running-the-pipeline)
  - [Generated artifacts](#generated-artifacts)
- [Limitations](#limitations)

## Repository layout

```text
clade-clustering/
├── manifests/
│   └── v4/approved_genomes.csv        # NCBI accessions approved by QC (contract between nb01 and nb02)
├── data/
│   ├── manual/
│   │   └── incoming/                  # 11 manually added FASTAs (Ca9–Ca14, SRR*.pilon.fasta)
│   └── reference/
│       ├── reference_sample_meta.csv  # expert clade labels used for anchoring
│       └── annotations/               # reference GFF (not versioned; downloaded from NCBI)
├── notebooks/
│   └── stage/                         # official pipeline, run in order
│       ├── 01_acquire_and_qc_dataset.ipynb
│       ├── 02_generate_parsnp_vcf.ipynb
│       ├── 03_build_snp_matrices.ipynb
│       ├── 04_epidemiology_and_ml_modeling.ipynb
│       └── 05_validate_clade_signatures.ipynb
├── src/clade_clustering/
│   ├── paths.py                       # fixed project paths (data/, manifests/, ...)
│   ├── layout.py                      # versioned directories and artifacts shared between stages
│   ├── genotypes.py                   # VCF → genotype matrix, locus table and sample metadata (nb03)
│   ├── metadata.py                    # NCBI metadata column normalization
│   └── scripts/                       # plotting and metric helpers (PCA, clustering)
├── tests/                             # pytest unit tests for layout.py and genotypes.py
├── pyproject.toml
└── uv.lock
```

Apart from the 11 manually added FASTAs in `data/manual/incoming/`, genomes, VCFs, matrices and results are **not versioned**. They live under `data/` and are regenerated by the notebooks (see [Input data](#input-data)). The notebooks in `notebooks/stage/` are committed with the outputs of their last run, so they can be read without executing anything. The notebook text and comments are in Portuguese.

## Pipeline

```text
NCBI Datasets ──► 01 QC ──► manifests/v4/approved_genomes.csv
                                     │
                                     ▼
data/manual/incoming ─────────► 02 Parsnp + HarvestTools ──► parsnp.vcf
                                                                 │
data/external/mapping_panel ──┐ (optional)                       ▼
                              └────────────────────────► 03 genotype matrices
                                                                 │
data/reference/reference_sample_meta.csv ──────────────► 04 PCA, KMeans, Ward, anchoring
                                                                 │
Reference GFF (GCF_003013715.1) ───────────────────────► 05 signatures, CV, annotation, external validation
```

| Notebook | What it does | Main outputs |
|---|---|---|
| `01_acquire_and_qc_dataset` | Downloads the full *C. auris* snapshot from NCBI, organizes the files, normalizes the metadata and applies QC (assemblies under 11 Mb, partial or atypical assemblies). | `data/processed/v4/sample_metadata_clean.parquet`, `excluded_genomes_report.csv`, `manifests/v4/approved_genomes.csv` |
| `02_generate_parsnp_vcf` | Downloads only the manifest's genomes (NCBI Datasets *dehydrate/rehydrate* flow), checks that the 11 manual FASTAs are present and adds them, runs Parsnp against `GCA_003013715.2` (B11220) and produces the VCF with HarvestTools. | `data/results/v4/parsnp/parsnp.vcf` |
| `03_build_snp_matrices` | Keeps `PASS` variants and simple SNPs, builds `X_genotypes`, `snp_info` and `sample_metadata`; if the optional external panel is present, also builds its matrix and genotype-call QC. | `data/results/v4/matrices/parsnp/*.parquet`, `data/results/v4/matrices/mapping_panel/*` |
| `04_epidemiology_and_ml_modeling` | One-hot encodes multiallelic loci, sets the 11 manual genomes aside as hold-outs, fits PCA (4 PCs), KMeans and Ward, and maps clusters to clades using the expert labels. Hold-outs are projected onto the fitted PCA and assigned by the fitted KMeans. | `X_encoded*.parquet`, `sample_metadata_labeled.parquet`, `sample_metadata_holdout.parquet`, `clustering/` |
| `05_validate_clade_signatures` | Discovers per-clade allelic signatures (thresholds 0.8 / 0.2), cross-validates them without letting validation experts into discovery, classifies the hold-outs (and the external panel, if present), annotates markers against the GFF (coding / genic non-coding / intergenic), computes per-gene recurrence and relates PCA loadings to genes. | `data/results/v4/signature_cross_validation/*.csv` |

`DATA_VERSION = "v4"` and `SEED = 42` are set at the top of each notebook.

## Python package

The notebook logic is being moved into the `clade_clustering` package, one stage at a time, without changing the scientific results. Each extracted function was checked against the `v4` artifacts before being adopted.

**`layout.py`** resolves the versioned directories and the files that one stage hands to the next (e.g. `sample_metadata_clean.parquet` from nb01 to nb03, `parsnp.vcf` from nb02 to nb03). It only builds paths; it never reads, writes or creates anything.

```python
from clade_clustering.layout import Layout

layout = Layout("v4")
layout.parsnp_vcf             # data/results/v4/parsnp/parsnp.vcf
layout.parsnp_matrices_dir    # data/results/v4/matrices/parsnp
```

**`genotypes.py`** reproduces the Parsnp part of notebook 03 (section 3.5.1). Processing and persistence are separate: the functions take explicit paths and do not depend on `Layout`.

| Function | Step |
|---|---|
| `read_vcf` | Reads the VCF; the 9 fixed columns are lowercased |
| `filter_pass_snps` | Keeps `PASS` sites whose REF and every ALT allele are single bases |
| `build_snp_info` | Builds the locus table, with `snp_id = chrom:pos:ref>alt` |
| `build_genotype_matrix` | Builds the sample × locus matrix (`int8`); rejects missing or non-`GT` genotypes |
| `map_samples_to_accessions` | Maps VCF sample names to NCBI accessions and drops the reference |
| `merge_sample_with_metadata` | Attaches `sample_metadata_clean` to each sample (left join) |
| `align_by_accession` | Indexes the matrix and the metadata by accession, row by row |
| `build_matrices_from_vcf` | Runs all the steps above, in order |
| `save_genotype_matrices` / `load_genotype_matrices` | Writes / reads `X_genotypes`, `snp_info` and `sample_metadata` as parquet |

```python
from clade_clustering import genotypes

X, snp_info, metadata = genotypes.build_matrices_from_vcf(
    vcf_path=layout.parsnp_vcf,
    sample_metadata_path=layout.sample_metadata_clean,
)
genotypes.save_genotype_matrices(X, snp_info, metadata, out_dir)
```

The external mapping panel (section 3.5.2 of notebook 03) has not been extracted yet. Unit tests run with `uv run pytest`.

## Reproducibility

### Requirements

- Python 3.12 and [uv](https://docs.astral.sh/uv/)
- [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) (`datasets`) on `PATH`, used by notebooks 01 and 02 (tested with 18.37.0)
- [Parsnp](https://github.com/marbl/parsnp) and [HarvestTools](https://github.com/marbl/harvest-tools) on `PATH`, used by notebook 02 only (see [Alignment tools](#alignment-tools-notebook-02))
- Disk: the full NCBI snapshot is about 4.9 GB compressed, and the Parsnp output directory takes about 12 GB

Notebooks 03 to 05 use only Python and run on an ordinary machine, as long as their inputs exist.

### Installation

```bash
git clone https://github.com/lucasofalonso/clade-clustering.git
cd clade-clustering
uv sync                 # creates .venv with the dependencies pinned in uv.lock
uv run jupyter lab
```

`uv sync` installs the `clade_clustering` package in editable mode; the notebooks import project paths (`clade_clustering.paths`) and helpers from it.

### Alignment tools (notebook 02)

Notebook 02 is meant to run on a compute server or cluster. For snapshot `v4` it ran inside a Docker container on a server, with Parsnp using 64 threads (`PARSNP_THREADS`); aligning 1,255 genomes took about 1.4 hours. Nothing in the notebook depends on that container: it only needs `datasets`, `parsnp` and `harvesttools` on `PATH`, and `PARSNP_THREADS` can be lowered for a smaller machine.

| Tool | Role in the pipeline | Version used for `v4` |
|---|---|---|
| Parsnp | Core-genome alignment of all genomes against B11220; writes `parsnp.ggr` | bioconda `parsnp` 2.1.6 (the binary reports itself as "Parsnp 2.1.5" in its log) |
| HarvestTools | Converts the Parsnp alignment (`parsnp.ggr`) into `parsnp.vcf` (`harvesttools -i parsnp.ggr -V parsnp.vcf`) | bioconda `harvesttools` 1.3 |

Both versions come from the container's conda revision history: the two packages were installed together on 2026-10-20 and were not changed before the run.

Parsnp and HarvestTools are not Python packages and are not on PyPI, so `uv sync` does not install them. To reproduce the same versions, install them from bioconda in a separate environment and activate it before starting Jupyter, so the binaries are on `PATH`:

```bash
conda create -n cauris-align -c conda-forge -c bioconda parsnp=2.1.6 harvesttools=1.3
conda activate cauris-align
uv run jupyter lab
```

If you already have `parsnp.vcf`, notebook 02 can be skipped entirely: place the file in `data/results/v4/parsnp/` and continue from notebook 03.

### Input data

| Input | Location | Versioned? | How to obtain |
|---|---|---|---|
| Approved NCBI genome list | `manifests/v4/approved_genomes.csv` | Yes | Produced by notebook 01 |
| Expert clade labels | `data/reference/reference_sample_meta.csv` | Yes | Included in the repository |
| NCBI genomes | `data/raw/v4/` (nb01) and `data/work/v4/` (nb02) | No | Downloaded by notebooks 01 and 02 |
| Manually added genomes (`Ca9`–`Ca14`, `SRR*.pilon.fasta`) | `data/manual/incoming/` | Yes | Included in the repository (see below) |
| External mapping panel (`cauris_mapping.vcf`, `IDs Cauris.xlsx`), optional | `data/external/mapping_panel/` | No | Produced externally by mapping reads against B11220; not distributed (see below) |
| Reference GFF annotation | `data/reference/annotations/` | No | Downloaded from NCBI (see [Reference genome annotation](#reference-genome-annotation)) |

**Manually added genomes.** The repository ships 11 FASTA files (about 12 MB each) in `data/manual/incoming/`:

| Files | Origin |
|---|---|
| `Ca9.fasta`, `Ca10.fasta`, `Ca11.fasta`, `Ca12.fasta`, `Ca13.fasta`, `Ca14.fasta` | Genomes provided by the lab |
| `SRR9007776.pilon.fasta`, `SRR18325429.pilon.fasta`, `SRR18325430.pilon.fasta`, `SRR18325431.pilon.fasta`, `SRR18325432.pilon.fasta` | Pilon-polished assemblies of public SRA runs |

These genomes are part of the Parsnp alignment, so they shape the core genome and the SNP set. They are then held out: notebook 04 fits PCA, KMeans and Ward on the 1,243 NCBI genomes only and projects the 11 genomes afterwards, and notebook 05 keeps them out of signature discovery and classifies them with the final panel. Notebook 02 stops with an explicit message if any of the 11 files is missing or if an unexpected FASTA is present.

**External mapping panel (optional).** The 30-genome panel used for external validation (5 genomes per clade) is a separate dependency and is not distributed with the repository. If you have it, place both files in `data/external/mapping_panel/`:

```text
data/external/mapping_panel/
├── cauris_mapping.vcf
└── IDs Cauris.xlsx
```

The main pipeline runs without the panel. Section 3.5.2 of notebook 03 (external matrix) and section 5.8 of notebook 05 (external validation) are the last sections of their notebooks; when the files are missing, each one prints which files it expected and stops the notebook there, after everything else has run and been saved. With the files present, both sections run unchanged.

**NCBI snapshot.** NCBI is updated continuously, so running notebook 01 today downloads a different set of genomes than the one used here (snapshot `v4` was downloaded on 2026-10-07). To reproduce the published results, start from notebook 02, which downloads exactly the accessions in `manifests/v4/approved_genomes.csv`. Notebook 03 also reads `data/processed/v4/sample_metadata_clean.parquet`, which is produced by notebook 01, so that file has to be generated first.

### Reference genome annotation

Notebook 05 annotates the signature markers and the top PCA loadings using the reference genome annotation:

| | |
|---|---|
| Organism | *Candida auris*, strain B11220 |
| Assembly (RefSeq) | `GCF_003013715.1` (`ASM301371v2`) |
| Format | GFF3 |
| Expected path | `data/reference/annotations/GCF_003013715.1_ASM301371v2_genomic.gff` |

`GCF_003013715.1` is the RefSeq counterpart of the GenBank assembly `GCA_003013715.2` used as the Parsnp reference. The two use different chromosome identifiers (`CP0435xx` in the VCF, `NC_0728xx` in the GFF), and notebook 05 maps one to the other.

The file is not versioned. Download it from the repository root with:

```bash
datasets download genome accession GCF_003013715.1 \
    --include gff3 \
    --filename cauris_annotation.zip

unzip cauris_annotation.zip -d cauris_annotation

mkdir -p data/reference/annotations

cp cauris_annotation/ncbi_dataset/data/GCF_003013715.1/genomic.gff \
   data/reference/annotations/GCF_003013715.1_ASM301371v2_genomic.gff
```

Notebook 05 checks for the file in its configuration cell and stops with a message pointing here if it is missing.

### Running the pipeline

Run the notebooks in `notebooks/stage/` in order (01 → 05), each from top to bottom. Each notebook:

- checks up front that its inputs exist and that the external executables are on `PATH`;
- refuses to overwrite heavy outputs that already exist (NCBI snapshot, Parsnp output directory, `parsnp.vcf`). To redo a step, delete or rename the corresponding output first;
- validates the shapes and index alignment of the artifacts it saves.

Notebooks 04 and 05 expect exactly 6 `Ca` and 5 `SRR` hold-out genomes and 1,243 fitting genomes; with a different input set these checks stop execution on purpose.

### Generated artifacts

```text
data/
├── raw/v4/                          # full NCBI snapshot (nb01)
├── processed/v4/                    # clean metadata and exclusion report (nb01)
├── work/v4/                         # rehydrated genomes and Parsnp input list (nb02)
└── results/v4/
    ├── parsnp/                      # parsnp.ggr, parsnp.vcf (nb02)
    ├── matrices/
    │   ├── parsnp/                  # X_genotypes, X_encoded, snp_info, sample_metadata* (nb03, nb04)
    │   └── mapping_panel/           # X_mapping_genotypes, genotype_call_qc (nb03)
    ├── clustering/                  # pca_loadings, expert_clade_labels (nb04)
    └── signature_cross_validation/  # CV metrics, panels, predictions, per-gene recurrence (nb05)
```

## Limitations

- **Clades V and VI are underrepresented.** The fitting set has 1 clade V genome and 3 clade VI genomes. Cross-validation covers clades I–IV only, and the V and VI signatures reflect the divergence of a handful of genomes rather than a population estimate.
- **Perfect CV scores.** *C. auris* clades are separated by thousands of fixed SNPs, so an accuracy of 1.0 is expected and does not by itself show that the classifier generalizes to low-quality genomes.
- **REF imputation in the external panel.** Uncalled genotypes and loci missing from the external VCF are set to the reference allele. Because the reference is clade II, low-coverage genomes tend to look like clade II; the `called_fraction` column in `genotype_call_qc.csv` lets you check this. The only error on the external panel (a clade VI genome predicted as clade I) is a partially covered genome.
- **Supervised classifiers disabled.** Section 4.8 (logistic regression and random forest) is kept in *raw* cells: with a single clade V genome, a stratified split is not possible.
