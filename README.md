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

For installation, the recommended setup is the portable conda environment file in this repository:

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

## Docker

CRISPRSCope can also be run from a Docker image. The image contains the
software environment, but your sequencing data, settings file, barcode file,
amplicon file, Bowtie2 index, and output directory should stay outside the
image and be mounted at runtime.

Build the image from the repository root:

```bash
docker build -t crisprscope:local .
```

That local image is useful for testing on your own computer, but it only targets
your current Docker architecture.

Check that the command-line tools are available:

```bash
docker run --rm crisprscope:local python -c "import CRISPRSCope; print(CRISPRSCope.__version__)"
docker run --rm crisprscope:local bowtie2 --version
docker run --rm crisprscope:local samtools --version
docker run --rm crisprscope:local CRISPResso --version
```

Run an analysis by mounting the folder that contains your settings file and
input data. In this example, everything is under the current directory and is
available inside the container as `/data`:

```bash
docker run --rm -v "$PWD:/data" crisprscope:local CRISPRSCope /data/example/example_settings.txt
```

For Docker Desktop from PowerShell, use `${PWD}`:

```powershell
docker run --rm -v "${PWD}:/data" crisprscope:local CRISPRSCope /data/example/example_settings.txt
```

The paths in the settings file must point to files visible inside the
container. Relative paths are resolved relative to the settings file, so keeping
the settings file, inputs, references, and results under the same mounted
directory is the simplest approach.

To publish a Docker Hub image that supports both Intel/AMD and ARM machines,
use Docker Buildx:

```bash
docker login
docker buildx build --platform linux/amd64,linux/arm64 -t DOCKERHUB_USERNAME/crisprscope:latest --push .
```

To publish both a versioned tag and `latest`:

```bash
docker buildx build --platform linux/amd64,linux/arm64 \
  -t DOCKERHUB_USERNAME/crisprscope:0.1.4 \
  -t DOCKERHUB_USERNAME/crisprscope:latest \
  --push .
```

This repository also includes a GitHub Actions workflow at
`.github/workflows/dockerhub.yml` that publishes a multi-architecture image to
Docker Hub. Add repository secrets named `DOCKERHUB_USERNAME` and
`DOCKERHUB_TOKEN`, then run the workflow manually or push a version tag such as
`v0.1.4`. The workflow builds and smoke-tests the `linux/amd64` image on a
native Intel/AMD GitHub runner before publishing the multi-architecture image.

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
write_editing_rate_ci	True
editing_rate_ci_bootstrap_iterations	10000
editing_rate_ci_confidence_level	0.95
editing_rate_ci_seed	42
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
| `ignore_substitutions` | `False` | Ignores substitution annotations when parsing CRISPResso read-alignment output and summarizing downstream editing calls. |
| `assign_reads_to_all_possible_amplicons` | `False` | If `True`, assigns ambiguous reads to every plausible amplicon. |
| `suppress_sub_crispresso_plots` | `False` | Disables per-amplicon CRISPResso plot/report generation. |
| `alt_alleles_file` | not used | Optional alternate allele definition file. |
| `min_total_reads_per_barcode` | `10` | Minimum total reads required for a barcode to be considered downstream. |
| `min_reads_per_amplicon_per_cell` | `0` | Minimum reads per amplicon per cell for scoring/filtering. |
| `write_editing_rate_ci` | `False` | Enables pointwise bootstrap confidence intervals for first-pass cell/allele editing rates. |
| `editing_rate_ci_bootstrap_iterations` | `10000` | Number of bootstrap resamples per amplicon; must be at least 100. |
| `editing_rate_ci_confidence_level` | `0.95` | Pointwise confidence level; must be greater than 0 and less than 1. |
| `editing_rate_ci_seed` | `42` | Non-negative base seed used for reproducible per-amplicon resampling. |
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

### Editing-Rate Confidence Intervals

When `write_editing_rate_ci` is enabled, CRISPRSCope resamples cells with replacement and computes pointwise percentile-bootstrap intervals from the first-pass inferred allele percentages in `editingSummary.txt`. Calls with missing modification percentages or coverage below `min_reads_per_amplicon_per_cell` are excluded independently for each amplicon.

The output reports three estimates for each amplicon:

- all analyzable cells with an eligible first-pass call
- cells in the quality categories enabled by the `include_*` settings
- the paired difference between the high-quality and all-cell estimates

These intervals quantify cell-sampling uncertainty within the current run; they do not represent uncertainty across biological replicates.

### h5ad Zygosity Parameters

These parameters control how the `.h5ad` export encodes zygosity calls.

| Key | Default |
| --- | --- |
| `h5ad_wt_max_mod_pct` | `20.0` |
| `h5ad_het_max_mod_pct` | `80.0` |
| `h5ad_hom_min_mod_pct` | `80.0` |
| `h5ad_compound_het_min_allele2_pct` | `20.0` |

## Expected Outputs

Given `output_root = results/demo_run`, you should expect outputs such as:

```text
results/demo_run.html
results/demo_run.log
results/demo_run.seq_by_amplicon/
results/demo_run.crispresso/
results/demo_run.crispresso.filtered/
results/demo_run.amplicon_score.txt
results/demo_run.filteredEditingSummary.txt
results/demo_run.filteredEditingSummaryPseudobulk.txt
results/demo_run.editingRateConfidenceIntervals.txt
results/demo_run.10_EditingRateConfidenceIntervals.{png,pdf}
results/demo_run.11_EditingRateQualityDelta.{png,pdf}
results/demo_run.h5ad
```

The exact set of plot PDFs, PNGs, and intermediate files depends on settings and on whether intermediate files are retained.

## Additional Notes

- The pipeline requires external command-line tools, especially `bowtie2` and `CRISPResso2`.
- The main workflow is designed for Linux-like environments such as Linux, WSL, or an HPC cluster.
- The settings file is strict about format: each non-comment line must contain exactly one key and one value separated by a tab.
- Multiple FASTQ pairs can be analyzed together by passing comma-separated file lists in `r1` and `r2`.
