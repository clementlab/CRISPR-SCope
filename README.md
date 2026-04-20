# CRISPRSCope

CRISPRSCope is a Python-based analysis pipeline for single-cell CRISPR DNA sequencing experiments. It takes paired-end FASTQ files, assigns reads to valid cell barcodes, maps reads to expected amplicons, runs CRISPResso2 on each target, summarizes editing outcomes across cells, generates QC and summary plots, and can export results to `.h5ad` for downstream analysis in Scanpy or related tools.

At a high level, the pipeline:

- reads paired-end FASTQ inputs
- validates and error-corrects cell barcodes
- assigns reads to amplicons using primer matching and genome alignment
- runs CRISPResso2 on per-amplicon read sets
- builds filtered per-cell editing summaries and QC reports
- optionally writes an `.h5ad` file for downstream single-cell analysis

## Installation

For collaborators, the recommended setup is the portable conda environment file in this repository:

```bash
conda env create -f environment.yml
conda activate crisprscope
```

If you want an editable local install of the package after activating the environment:

```bash
pip install -e .
```

A simple import sanity check is:

```bash
python -c "import CRISPRSCope; print(CRISPRSCope.__version__)"
```

## Running The Pipeline

CRISPRSCope expects a tab-delimited settings file as its main input:

```bash
CRISPRSCope path/to/run_settings.txt
```

The first positional argument must be the settings file. A log file is written next to that settings file as:

```text
path/to/run_settings.txt.log
```

## Required Inputs

Before running the pipeline, you should have:

- paired-end FASTQ files (`r1` and `r2`)
- a barcode whitelist file with one valid barcode per line
- an amplicon definition file
- a Bowtie2 genome index prefix
- a tab-delimited settings file that points to the above inputs

## Example Input Files

These are visual examples only. They are included here to show the expected structure and are not shipped as runnable project files.

### Example Settings File

The settings file must be tab-delimited, with one `key<TAB>value` entry per line.

```tsv
r1	data/sample_A_R1.fastq.gz,data/sample_B_R1.fastq.gz
r2	data/sample_A_R2.fastq.gz,data/sample_B_R2.fastq.gz
constant1	GTTTAAGAGCTATGCTGGAAACAG
constant2	GTTTTAGAGCTAGAAATAGCAAGT
barcodes	inputs/barcodes.txt
amplicons	inputs/amplicons.tsv
bowtie2_index	references/hg38/hg38
output_root	results/demo_run
processes	8
allowBarcodeMismatches	True
keep_intermediate_files	False
ignore_substitutions	False
assign_reads_to_all_possible_amplicons	False
suppress_sub_crispresso_plots	False
min_total_reads_per_barcode	10
min_reads_per_amplicon_per_cell	0
include_high_score_high_depth	True
include_high_score_low_depth	True
include_low_score_high_depth	False
include_low_score_low_depth	False
write_h5ad	True
h5ad_output	results/demo_run.h5ad
h5ad_wt_max_mod_pct	20
h5ad_het_max_mod_pct	80
h5ad_hom_min_mod_pct	80
h5ad_compound_het_min_allele2_pct	20
```

### Example Barcode File

The barcode file is one barcode per line.

```text
AACCGGTTAA
AACCGGTTAC
AACCGGTTAG
TTGCAACGTA
TTGCAACGTC
TTGCAACGTG
```

### Example Amplicon File

The amplicon file is tab-delimited. The first two columns are required:

1. amplicon name
2. amplicon sequence

Optional columns currently supported by the pipeline are:

3. guide sequence
4. reference allele count

```tsv
AMP_TARGET_1	ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG	GGACTGACTGACTGACTGA	2
AMP_TARGET_2	TGACCTGATCGATCGTAGCTAGCTAGCTAGCATCGATCGA	CTGATCGATCGTAGCTAGC	2
AMP_TARGET_3	GGCTAACCGGTTAACCGGTTAACCGGTTAACTGACTGACT	ACCGGTTAACCGGTTAACT	2
```

### Example Alternate Alleles File

This file is optional. The current implementation expects a header row and uses:

- column 1: amplicon name
- column 3: comma-separated alternate allele sequences

```tsv
amplicon_name	label	alternate_alleles
AMP_TARGET_1	edited	ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG,ACTGACTGACTGACTGACTG---TGACTGACTGACTG
AMP_TARGET_2	edited	TGACCTGATCGATCGTAGCTAGCTAGCTAGCATCGATCGA,TGACCTGATCGATCGTAGCTAG---GCTAGCATCGATCGA
```

## Settings Reference

### Required Settings

| Key | Description |
| --- | --- |
| `r1` | Comma-separated list of R1 FASTQ files. |
| `r2` | Comma-separated list of matching R2 FASTQ files. |
| `constant1` | First constant sequence used during barcode/read parsing. |
| `constant2` | Second constant sequence used during barcode/read parsing. |
| `barcodes` | Path to barcode whitelist file. |
| `amplicons` | Path to tab-delimited amplicon definition file. |
| `bowtie2_index` | Bowtie2 index prefix. |

You may also use `genome` instead of `bowtie2_index`; internally the pipeline resolves either key to the Bowtie2 index prefix.

### Optional Settings

| Key | Default | Description |
| --- | --- | --- |
| `output_root` | settings file path | Prefix used for generated outputs. |
| `processes` | all available CPUs | Number of processes to use. |
| `allowBarcodeMismatches` | off | Enables single-mismatch barcode rescue. |
| `keep_intermediate_files` | `False` | Keeps intermediate files instead of cleaning them up. |
| `ignore_substitutions` | `False` | Ignores substitutions during downstream editing summarization. |
| `assign_reads_to_all_possible_amplicons` | `False` | If `True`, assigns ambiguous reads to every plausible amplicon. |
| `suppress_sub_crispresso_plots` | `False` | Disables per-amplicon CRISPResso plot/report generation. |
| `alt_alleles_file` | not used | Optional alternate allele definition file. |
| `min_total_reads_per_barcode` | `10` | Minimum total reads required for a barcode to be considered downstream. |
| `min_reads_per_amplicon_per_cell` | `0` | Minimum reads per amplicon per cell for scoring/filtering. |
| `write_h5ad` | `True` | Enables `.h5ad` export after the main run. |
| `h5ad_output` | `<output_root>.h5ad` | Output path for the generated `.h5ad` file. |

### Cell-Quality Inclusion Flags

If none of these flags are provided, the pipeline defaults to including only `HQ_HI`.

| Key | Meaning |
| --- | --- |
| `include_high_score_high_depth` | Include high-score, high-depth cells (`HQ_HI`). |
| `include_high_score_low_depth` | Include high-score, low-depth cells (`HQ_LO`). |
| `include_low_score_high_depth` | Include low-score, high-depth cells (`LQ_HI`). |
| `include_low_score_low_depth` | Include low-score, low-depth cells (`LQ_LO`). |

### h5ad Zygosity Parameters

These parameters control how the `.h5ad` export encodes zygosity calls.

| Key | Default |
| --- | --- |
| `h5ad_wt_max_mod_pct` | `20.0` |
| `h5ad_het_max_mod_pct` | `80.0` |
| `h5ad_hom_min_mod_pct` | `80.0` |
| `h5ad_compound_het_min_allele2_pct` | `20.0` |

## Expected Outputs

Given `output_root = results/demo_run`, collaborators should expect outputs such as:

```text
results/demo_run.html
results/demo_run.log
results/demo_run.seq_by_amplicon/
results/demo_run.crispresso/
results/demo_run.crispresso.filtered/
results/demo_run.amplicon_score.txt
results/demo_run.filteredEditingSummary.txt
results/demo_run.filteredEditingSummaryPseudobulk.txt
results/demo_run.h5ad
```

The exact set of plot PDFs, PNGs, and intermediate files depends on settings and on whether intermediate files are retained.

## Notes For Collaborators

- The pipeline requires external command-line tools, especially `bowtie2` and `CRISPResso2`.
- The main workflow is designed for Linux-like environments such as Linux, WSL, or an HPC cluster.
- The settings file is strict about format: each non-comment line must contain exactly one key and one value separated by a tab.
- Multiple FASTQ pairs can be analyzed together by passing comma-separated file lists in `r1` and `r2`.
