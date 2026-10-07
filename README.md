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
cache_mode	auto
allowBarcodeMismatches	True
keep_intermediate_files	False
ignore_substitutions	False
assign_reads_to_all_possible_amplicons	False
suppress_sub_crispresso_plots	False
min_total_reads_per_barcode	10
min_reads_per_amplicon_per_cell	0
min_reads_per_amplicon_for_genotype	8
min_allele_support	2
amplicon_score_min_reads_per_amplicon	5
amplicon_score_min_covered_fraction	0.6666666666666666
amplicon_score_max_barcode_rank	10000
include_high_score_high_depth	True
include_high_score_low_depth	True
include_low_score_high_depth	False
include_low_score_low_depth	False
write_editing_rate_ci	True
editing_rate_ci_bootstrap_iterations	10000
editing_rate_ci_permutation_iterations	10000
editing_rate_ci_confidence_level	0.95
editing_rate_ci_seed	42
editing_rate_ci_coverage_exact_max_reads	10
editing_rate_ci_coverage_bin_width_reads	5
write_h5ad	True
write_output_manifest	False
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
| `cache_mode` | `auto` | Processing-stage cache policy: `auto` validates and reuses records, `refresh` recomputes and replaces them, and `disabled` recomputes without reading or writing records. |
| `allowBarcodeMismatches` | off | Enables single-mismatch barcode rescue. |
| `keep_intermediate_files` | `False` | Keeps intermediate files instead of cleaning them up. |
| `ignore_substitutions` | `False` | Ignores substitution annotations when parsing CRISPResso read-alignment output and summarizing downstream editing calls. |
| `assign_reads_to_all_possible_amplicons` | `False` | If `True`, assigns ambiguous reads to every plausible amplicon. |
| `suppress_sub_crispresso_plots` | `False` | Disables per-amplicon CRISPResso plot/report generation. |
| `alt_alleles_file` | not used | Optional alternate allele definition file. |
| `min_total_reads_per_barcode` | `10` | Minimum total reads required for a barcode to be considered downstream. |
| `min_reads_per_amplicon_per_cell` | `0` | Optional stricter gate requiring this many reads at every usable amplicon before a barcode is scored; it also sets per-amplicon eligibility in editing-rate analyses. |
| `min_reads_per_amplicon_for_genotype` | `8` | Minimum accepted reads at one cell–amplicon for a genotype. Below this depth, genotype `modPct` is `NA` without removing the cell or its read evidence. |
| `min_allele_support` | `2` | A whole number sets minimum reads for a candidate allele (`2`); decimal syntax sets a fraction of accepted reads at that cell–amplicon (`0.1` means 10%). Only one mode applies. The unused `min_allele_count_cutoff` and `min_allele_pct_cutoff` names are retired. |
| `amplicon_score_min_reads_per_amplicon` | `5` | Reads required for an amplicon to count as supported in the breadth score; must be at least 1. |
| `amplicon_score_min_covered_fraction` | `0.6666666666666666` | Fraction of usable amplicons that must be supported for a high score; must be greater than 0 and no greater than 1. |
| `amplicon_score_max_barcode_rank` | `10000` | Largest total-read barcode rank classified as high depth; must be at least 1. |
| `write_editing_rate_ci` | `True` | Enables pointwise bootstrap confidence intervals for first-pass edited-cell rates; set to `False` to disable. |
| `editing_rate_ci_bootstrap_iterations` | `10000` | Number of bootstrap resamples per amplicon; must be at least 100. |
| `editing_rate_ci_permutation_iterations` | `10000` | Number of configured analysis-group label permutations used for each two-sided significance test; must be at least 100. |
| `editing_rate_ci_confidence_level` | `0.95` | Pointwise confidence level; must be greater than 0 and less than 1. |
| `editing_rate_ci_seed` | `42` | Non-negative base seed used for reproducible per-amplicon resampling. |
| `editing_rate_ci_coverage_exact_max_reads` | `10` | Highest per-amplicon read count kept as an exact coverage stratum for the coverage-controlled test; must be non-negative. |
| `editing_rate_ci_coverage_bin_width_reads` | `5` | Width of coverage strata above the exact-count ceiling; must be at least 1. With the defaults, the first binned strata are 11–15, 16–20, and 21–25 reads. |
| `write_editing_rate_depth_stability` | `False` | Enables the optional finite-cohort downsampling table and detailed plot 14 for AllCells and InGroup. |
| `editing_rate_depth_stability_iterations` | `1000` | Number of without-replacement subsamples at each retained-cell percentage; must be at least 100. |
| `editing_rate_depth_stability_percentages` | `10,25,50,75,90` | Strictly increasing, unique retained-cell percentages between 0 and 100; an exact 100% reference is added automatically. |
| `write_h5ad` | `True` | Enables `.h5ad` export after the main run. |
| `h5ad_output` | `<output_root>.h5ad` | Output path for the generated `.h5ad` file. |
| `write_output_manifest` | `False` | Writes `<output_root>.outputManifest.json`, an ordered diagnostic inventory of output status, paths, data links, and any pipeline failure. |

### Cell-Quality Inclusion Flags

If none of these flags are provided, the pipeline defaults to including only `HQ_HI`.

The amplicon score is the fraction of usable first-pass amplicons meeting
`amplicon_score_min_reads_per_amplicon`. A barcode is high score when it has
at least `ceil(amplicon_score_min_covered_fraction × usable amplicons)`
supported amplicons. Reads beyond the support threshold at one amplicon do not
increase its contribution, so isolated amplification jackpots cannot compensate
for missing coverage elsewhere in the panel. High versus low depth remains a
separate classification based on `amplicon_score_max_barcode_rank`.

The `.amplicon_score.txt` table reports `Amplicon Score` on a 0–1 scale together
with `Supported Amplicons`, `Usable Amplicons`, total read count, barcode rank,
and the four-way quality classification. It is regenerated on every run so
changed scoring settings cannot silently reuse stale classifications.

| Key | Meaning |
| --- | --- |
| `include_high_score_high_depth` | Include high-score, high-depth cells (`HQ_HI`). |
| `include_high_score_low_depth` | Include high-score, low-depth cells (`HQ_LO`). |
| `include_low_score_high_depth` | Include low-score, high-depth cells (`LQ_HI`). |
| `include_low_score_low_depth` | Include low-score, low-depth cells (`LQ_LO`). |

### Editing-Rate Confidence Intervals

By default, CRISPRSCope resamples cells with replacement and computes pointwise percentile-bootstrap intervals for the percentage of cells with at least one edited allele. In non-pseudobulk `editingSummary.txt`, `modPct` is a categorical genotype encoding: `0` (WT/WT) is unedited and `50` (WT/Mut) or `100` (Mut/Mut) is edited. Missing calls are excluded independently for each amplicon, as are calls below `min_reads_per_amplicon_per_cell`; any other finite `modPct` value is rejected.

Allele calls use a multinomial model with a fixed 1% combined noise category and the expected copy count for each reference. Genotype depth counts only reads with an unambiguous reference alignment, matching amplicon arms, and parsable edit annotation. Read-count ties at the selection boundary use the median of each supporting read's lowest Phred quality near its annotated quantification-window edit; unresolved ties prefer WT, then a stable allele ordering. The `alleleCallQC.txt` table records accepted depth, selected and competing allele support, tie resolution, and whether a genotype was withheld. Indel-flank Phred quality is a tie breaker, not an alignment-confidence estimate. Raw-read pseudobulk `modPct` remains available when genotype `modPct` is `NA`. The defaults of 8 reads and 2 allele-supporting reads are guard rails, not biologically validated thresholds; compare 6 versus 8 reads on a representative dataset before manuscript use.

The output reports estimates for each amplicon from:

- **AllCells:** all analyzable cells with an eligible first-pass call
- **InGroup:** eligible cells in the quality categories enabled by the `include_*` settings
- **OutGroup:** every other eligible cell
- the InGroup-minus-AllCells and InGroup-minus-OutGroup differences
- a coverage-controlled InGroup-minus-OutGroup difference

For each amplicon, CRISPRSCope first permutes InGroup labels among AllCells without replacement while preserving the observed InGroup size. This unconditional two-sided test asks whether the InGroup behaves differently from a random same-sized subset. The exact permutation draws used for the p-value are also written to `.editingRateUnconditionalPermutationSimulations.txt`, summarized in `.editingRateUnconditionalPermutation.txt`, and shown in `.12_EditingRateUnconditionalPermutation.{png,pdf}` for every estimable amplicon.

The same draws are also shown in `.13_EditingRateObservedCenteredPermutationSwarm.{png,pdf}`. Each point is a simulated InGroup-sized subset edited-cell rate shown relative to that amplicon's observed InGroup edited-cell rate, in percentage points; the dashed zero line marks the observed rate. This is a second view of the existing simulation output, not an additional resampling procedure.

The coverage-controlled follow-up runs for every testable amplicon. Read counts through `editing_rate_ci_coverage_exact_max_reads` define exact strata; higher counts are grouped into consecutive bins of `editing_rate_ci_coverage_bin_width_reads`. Labels are permuted only within strata containing both InGroup and OutGroup cells. The adjusted effect is an information-weighted average of the within-stratum InGroup-minus-OutGroup differences. Cells in single-cohort coverage strata remain in the unconditional estimates but cannot contribute to the controlled effect; common-support counts and retained percentages are reported explicitly. The table also reports common-depth standardized InGroup and OutGroup means using the same normalized overlap weights; their difference equals the coverage-adjusted effect. A separate within-stratum bootstrap supplies the adjusted effect's confidence interval.

The unconditional and controlled permutation p-values receive separate Benjamini-Hochberg adjustments across all testable amplicons. The confidence-interval figure and optional editing-stability figure are filtered using the coverage-controlled adjusted p-value at or below 0.05, while the raw-versus-adjusted comparison figure and full output tables retain all estimable amplicons. If no amplicon passes the controlled threshold, the significance-filtered figures are omitted but the comparison figure is retained when possible.

These intervals and significance tests quantify cell-sampling and cell-selection behavior within the current run. The controlled effect applies only to coverage ranges represented in both populations. It does not represent uncertainty across biological replicates or establish a biological or causal effect of cell-quality selection.

### Editing-Rate Cell-Depth Stability

Set `write_editing_rate_depth_stability` to `True` to request this optional detailed analysis; it is disabled by default. CRISPRSCope then repeatedly downsamples eligible cells without replacement, independently for AllCells and InGroup. Requested percentages are applied to each amplicon's eligible cohort, so the actual number of sampled cells is reported for every result.

Within each iteration, the percentage levels are nested: one random ordering of eligible cells supplies the first 10%, 25%, 50%, 75%, and 90%. The output reports the median edited-cell rate, a central interval controlled by `editing_rate_ci_confidence_level`, and absolute deviations from the full-cohort estimate. An exact 100% reference is appended automatically. Random sampling uses `editing_rate_ci_seed`, making identical inputs and settings reproducible across serial and parallel runs.

Plot 14 reports percentage-point deviations for coverage-controlled significant amplicons with a usable cohort and orders them by the full InGroup edited-cell rate, highest first. If no amplicon passes the coverage-controlled threshold, the table is retained and plot 14 is omitted.

The stability bands answer how much the inferred edited-cell rate changes as cells from this run are retained or removed. They are finite-cohort downsampling diagnostics, not confidence intervals across biological replicates.

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
results/demo_run.alleleCallQC.txt
results/demo_run.filteredEditingSummary.txt
results/demo_run.filteredEditingSummaryPseudobulk.txt
results/demo_run.editingRateConfidenceIntervals.txt
results/demo_run.editingRateUnconditionalPermutation.txt
results/demo_run.editingRateUnconditionalPermutationSimulations.txt
results/demo_run.10_EditingRateConfidenceIntervals.{png,pdf}
results/demo_run.11_EditingRateCoverageAdjustedEffects.{png,pdf}
results/demo_run.12_EditingRateUnconditionalPermutation.{png,pdf}
results/demo_run.13_EditingRateObservedCenteredPermutationSwarm.{png,pdf}
results/demo_run.h5ad
```

When `write_output_manifest=True`, CRISPRSCope also writes
`results/demo_run.outputManifest.json`. This compact diagnostic inventory
records each declared final artifact's path and lifecycle status, linked data
artifacts, intermediate-directory summaries, and the active stage/error if a
run fails. Manifest schema 2 also records the cache mode, decision summary,
and ordered per-stage cache events.

## Processing Cache and Resume Behavior

CRISPRSCope stores versioned cache records under `<output_root>.cache/` for
parse/alignment, read splitting, both CRISPResso passes, CRISPResso parsing,
and selected-cell FASTQ filtering. CRISPResso-related records are independent
per amplicon, so a missing or changed output reruns only that amplicon and its
dependents. Aggregate tables, resampling analyses, plots, reports, and `.h5ad`
exports are regenerated on every run.

`cache_mode=auto` is the default. Use `cache_mode=refresh` to force every
managed stage to recompute and replace its record. Use `cache_mode=disabled`
to recompute without consulting or changing cache records. Old `.finished`,
`.summ.finished`, split TSV, and barcode-digest files remain useful diagnostics
but cannot create a cache hit without a current JSON record. Runs created before
this cache format therefore incur one conservative recomputation.

Cache validation avoids payload-size work for sequencing artifacts. Raw input
FASTQs, BAMs, Bowtie2 indexes, and large CRISPResso outputs use resolved path,
byte size, and nanosecond modification time. Per-amplicon gzip FASTQs produced
by read splitting use resolved path, compressed size, and the gzip trailer's
CRC32 and uncompressed size. This constant-time signature allows an unchanged
amplicon to remain cached after a run-wide split rerun while reading only the
gzip header and trailer. BAMs additionally undergo `samtools quickcheck`, and
gzipped FASTQs receive a gzip magic-byte check.
CRISPResso output FASTQs are accepted in either gzip form or as plain FASTQ,
because supported CRISPResso releases may write plain text with a `.fastq.gz`
suffix; this check reads only the format signature.
Small control and summary files use SHA-256. Consequently, an external process
that changes a stat-fingerprinted large file while preserving both its size and
modification time can evade validation. The gzip trailer signature is also not
a cryptographic content digest. Run with `cache_mode=refresh` when either case
is possible.

Only one process may use a given output root at a time. CRISPRSCope acquires a
non-blocking `<output_root>.cache.lock`; contention fails immediately and does
not replace the active run's manifest. Cache decisions are also written to the
normal log as `CACHE HIT`, `CACHE MISS`, `CACHE INVALID`, `CACHE REFRESH`, or
`CACHE DISABLED` lines.

When `write_editing_rate_depth_stability=True`, the additional detailed outputs are:

```text
results/demo_run.editingRateDepthStability.txt
results/demo_run.14_EditingRateDepthStability.{png,pdf}
```

The exact set of plot PDFs, PNGs, and intermediate files depends on settings and on whether intermediate files are retained.

## Additional Notes

- The pipeline requires external command-line tools, especially `bowtie2` and `CRISPResso2`.
- The main workflow is designed for Linux-like environments such as Linux, WSL, or an HPC cluster.
- The settings file is strict about format: each non-comment line must contain exactly one key and one value separated by a tab.
- Multiple FASTQ pairs can be analyzed together by passing comma-separated file lists in `r1` and `r2`.
