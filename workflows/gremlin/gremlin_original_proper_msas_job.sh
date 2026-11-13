#!/bin/bash

#SBATCH --job-name=gremlin_original_proper_msas_job
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --partition=norm
#SBATCH --ntasks=1

# Move to data directory & load required software module
cd /data/schwartzlc/gremlin_test/scripts

# Load required software
set -euo pipefail

PYTHON_MODULE="${PYTHON_MODULE:-python/3.12}"

module purge
module load "$PYTHON_MODULE"

python - <<'PY'
import sys
print("python", sys.version.replace("\n", " "))

import numpy
print("numpy", numpy.__version__)

import scipy
print("scipy", scipy.__version__)

import pandas
print("pandas", pandas.__version__)

import matplotlib
print("matplotlib", matplotlib.__version__)

import tensorflow.compat.v1 as tf
print("tensorflow", tf.__version__)
PY

# Set environmental variables
export GREMLIN_OUTDIR="/data/schwartzlc/gremlin_test/results"

# Loop through .fas files in a specific folder & run the script
for file in /data/AF_interpretation/msa_swaps/original_proper_msas/*/*.fas; do \
    export GREMLIN_FASTA="$file"; \
    echo "Processing $file..."; \
    python gremlin_biowulf_env.py "$file"; \
done

