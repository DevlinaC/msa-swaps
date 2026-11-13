# GREMLIN inference source

`gremlin_biowulf_env.py` is the locally supplied TensorFlow-v1-compatible runner. It reads settings from environment variables, preserves original column indices and writes full pair scores. The original launchers and `READ_ME.md` are retained as historical provenance; they use fixed cluster paths and should not be submitted unchanged.

Within a suitable tested environment and, on Biowulf, a compute allocation:

```bash
export GREMLIN_FASTA=/absolute/path/input.fas
export GREMLIN_OUTDIR=/absolute/path/new-results
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS
python workflows/gremlin/gremlin_biowulf_env.py
```

Use `GREMLIN_DRY_RUN=1` for preprocessing-only inspection. The original default generation gap cutoff is 0.5, distinct from the downstream 0.25 analysis mask. The inference dependency includes `tensorflow.compat.v1`; it is deliberately absent from the lightweight plotting requirements. Exact tested inference versions and upstream license/source attribution still need final reconciliation. Do not equate availability of this runner with proof that every historical matrix used this exact version.

The top-100 cap was in downstream `gremlin_jaccard_v2.py`, not a limit on the full GREMLIN matrix. S2 now uses 100/100/150; S3/S4 use their own documented top-L definitions.
