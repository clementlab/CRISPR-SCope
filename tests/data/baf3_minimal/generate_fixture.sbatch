#!/bin/bash
#SBATCH --job-name=baf3-fixture
#SBATCH --partition=notchpeak-shared-short
#SBATCH --qos=notchpeak-shared-short
#SBATCH --account=notchpeak-shared-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/general/vast/%u/baf3-fixture.%j.out
#SBATCH --error=/scratch/general/vast/%u/baf3-fixture.%j.err

set -euo pipefail

PROJECT=/uufs/chpc.utah.edu/common/home/clementm-group1/projects/20230818_scCRISPR/analysis/CRISPR-SCope
WORKSPACE=/uufs/chpc.utah.edu/common/home/clementm-group1/projects/20230818_scCRISPR
SCRATCHDIR=/scratch/general/vast/${USER}/baf3-fixture-${SLURM_JOB_ID}
SOURCE_RUN=${WORKSPACE}/analysis/01_run_on_20200804_BaF3_revision
SOURCE_DATA=${WORKSPACE}/data/20200804_BaF3_revision/data
FIXTURE=${PROJECT}/tests/data/baf3_minimal

mkdir -p "${SCRATCHDIR}"
cp -r "${SOURCE_RUN}/settings.txt.seq_by_amplicon" "${SCRATCHDIR}/"
mkdir -p "${SCRATCHDIR}/data"
cp "${SOURCE_DATA}/BaF3-NSG_S1_L001_R1_001.fastq.gz" "${SCRATCHDIR}/data/"
cp "${SOURCE_DATA}/BaF3-NSG_S1_L001_R2_001.fastq.gz" "${SCRATCHDIR}/data/"
cp "${WORKSPACE}/data/20200804_BaF3_revision/amplicons.txt" "${SCRATCHDIR}/"
cp "${SOURCE_RUN}/settings.txt.amplicon_score.txt" "${SCRATCHDIR}/"

/uufs/chpc.utah.edu/common/home/clementm-group1/conda/mambaforge/bin/conda run --no-capture-output \
  -n sc_CRISPRSCope python "${FIXTURE}/generate_fixture.py" \
  --destination "${SCRATCHDIR}/fixture" \
  --source-run "${SCRATCHDIR}" \
  --source-data "${SCRATCHDIR}/data"

cp "${SCRATCHDIR}/fixture"/R1.fastq.gz "${FIXTURE}/"
cp "${SCRATCHDIR}/fixture"/R2.fastq.gz "${FIXTURE}/"
cp "${SCRATCHDIR}/fixture"/barcodes.txt "${FIXTURE}/"
cp "${SCRATCHDIR}/fixture"/amplicons.txt "${FIXTURE}/"
cp "${SCRATCHDIR}/fixture"/minigenome* "${FIXTURE}/"
cp "${SCRATCHDIR}/fixture"/settings.txt "${FIXTURE}/"
cp "${SCRATCHDIR}/fixture"/manifest.json "${FIXTURE}/"
