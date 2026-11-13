#!/bin/bash
#SBATCH --job-name=pairformer512-off
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1,lscratch:20
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=00:45:00
#SBATCH --output=logs/inference-%j.out
#SBATCH --error=logs/inference-%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
module load python/3.11 hhsuite/3.3.0
export TMPDIR="/lscratch/$SLURM_JOB_ID"
test -d "$TMPDIR"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS" MKL_NUM_THREADS="$OMP_NUM_THREADS" NUMEXPR_NUM_THREADS="$OMP_NUM_THREADS"
export PYTHONNOUSERSITE=1 HF_HUB_OFFLINE=1
unset PYTHONPATH
export MPLCONFIGDIR="$TMPDIR/matplotlib" XDG_CACHE_HOME="$TMPDIR/cache" TORCH_EXTENSIONS_DIR="$TMPDIR/torch_extensions"
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
exec 9>.inference.lock
flock -n 9 || { echo "Inference already running" >&2; exit 1; }
PILOT="${PAIRFORMER_PILOT:-$(dirname "$PWD")/msa_pairformer_pilot}"
exec 8>"$PILOT/.environment.lock"
flock -n 8 || { echo "Pilot environment is being changed or used" >&2; exit 1; }
module list 2>&1
nvidia-smi
"$PILOT/env/bin/python" scripts/run_full.py --root "$PWD" --pilot "$PILOT"
