# Minimal BaF3 integration fixture

This directory contains a deterministic, 30-amplicon fixture for the full
CRISPRSCope pipeline. It contains original read pairs selected from the BaF3
run, not synthetic reads. The fixture selects 16 cells from each `HQ_HI`,
`HQ_LO`, `LQ_HI`, and `LQ_LO` category using a fixture-only
`amplicon_score_max_barcode_rank` of 32. This deliberately falls below the
64-cell fixture size, so the run exercises both high- and low-depth categories.
The cell panel itself is balanced across the completed source run's categories
at its original rank cutoff (6,000), since the source and fixture rank scales
are different. Each selected cell contributes at
most 450 deterministic read pairs, preserving its original amplicon breadth.
The mini Bowtie2 index contains only the 30 amplicon reference sequences and
is valid only for this controlled fixture. CRISPResso sub-plots remain enabled
because CRISPResso 2.3.3 fails after analysis when plot generation is suppressed
alongside FASTQ output; these temporary plot files are not golden-checked.
The fixture pins the cell-amplicon genotype depth to 8 reads and minimum allele
support to 2 reads. Golden hashes include allele-call QC, genotype and
pseudobulk summaries, amplicon scores, and editing-rate results.

Regenerate from the workspace root with the project environment active:

```bash
python analysis/CRISPR-SCope/tests/data/baf3_minimal/generate_fixture.py --overwrite
```

After a clean fixture run, record intentional golden-baseline changes with:

```bash
python generate_fixture.py --record-golden run --overwrite
```

Ordinary test runs exclude this fixture. Run it explicitly with
`pytest -m integration`.

The fixture uses `cache_mode=auto`. The integration test performs a cold run
followed by an identical warm run and verifies that all managed processing
stages hit their records while the final deliverables are regenerated.
