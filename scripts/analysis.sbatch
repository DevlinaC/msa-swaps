#!/bin/bash
#SBATCH --job-name=msa-swaps-analysis
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail
: "${PYTHON_BIN:?Set PYTHON_BIN to the absolute Python executable in a prepared /data environment}"
[[ "$PYTHON_BIN" = /* && -x "$PYTHON_BIN" ]] || { echo "Invalid PYTHON_BIN" >&2; exit 1; }
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root under /data}"
test -f scripts/validate_examples.py
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export MPLCONFIGDIR="$PWD/outputs/.matplotlib"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$MPLCONFIGDIR"
echo "job_id=$SLURM_JOB_ID cpus=$OMP_NUM_THREADS"
"$PYTHON_BIN" --version
"$PYTHON_BIN" scripts/validate_examples.py
"$PYTHON_BIN" scripts/compute_contacts.py
"$PYTHON_BIN" scripts/plot_results.py --out-dir outputs/plots
