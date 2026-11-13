#!/bin/bash
#SBATCH --job-name=pairformer512-prep
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --gres=lscratch:10
#SBATCH --time=01:00:00
#SBATCH --output=logs/prepare-%j.out
#SBATCH --error=logs/prepare-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
module load python/3.11
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1
unset PYTHONPATH
export TMPDIR="/lscratch/$SLURM_JOB_ID"
test -d "$TMPDIR"
exec 9>.prepare.lock
flock -n 9 || { echo "Preparation already running" >&2; exit 1; }
PILOT="${PAIRFORMER_PILOT:-$(dirname "$PWD")/msa_pairformer_pilot}"
"$PILOT/env/bin/python" scripts/prepare_full.py --root "$PWD" --source "${MSA_SOURCE_ROOT:-/data/AF_interpretation/msa_swaps}"
