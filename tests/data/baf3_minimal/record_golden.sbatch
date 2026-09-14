#!/bin/bash
#SBATCH --job-name=baf3-golden
#SBATCH --partition=notchpeak-shared-short
#SBATCH --qos=notchpeak-shared-short
#SBATCH --account=notchpeak-shared-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/general/vast/%u/baf3-golden.%j.out
#SBATCH --error=/scratch/general/vast/%u/baf3-golden.%j.err

set -euo pipefail

PROJECT=/uufs/chpc.utah.edu/common/home/clementm-group1/projects/20230818_scCRISPR/analysis/CRISPR-SCope
FIXTURE=${PROJECT}/tests/data/baf3_minimal
SCRATCHDIR=/scratch/general/vast/${USER}/baf3-golden-${SLURM_JOB_ID}

mkdir -p "${SCRATCHDIR}"
cp -r "${FIXTURE}" "${SCRATCHDIR}/fixture"
cd "${SCRATCHDIR}/fixture"

/uufs/chpc.utah.edu/common/home/clementm-group1/conda/mambaforge/bin/conda run --no-capture-output \
  -n sc_CRISPRSCope python -m CRISPRSCope.cli settings.txt
/uufs/chpc.utah.edu/common/home/clementm-group1/conda/mambaforge/bin/conda run --no-capture-output \
  -n sc_CRISPRSCope python generate_fixture.py --destination . --record-golden run --overwrite

cp golden.json "${FIXTURE}/golden.json"
