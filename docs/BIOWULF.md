# NIH Biowulf execution

Use a Slurm compute allocation for scientific computation, substantial tests and environment installation. Login-node work is limited to light inspection/editing and submission. Transfer through Helix or Globus. Keep persistent inputs, code and results under `/data`; request `/lscratch` for substantial temporary I/O. Use pinned modules/environments and bound threads to allocated CPUs.

The portable plotting command can run locally. `scripts/analysis.sbatch` is a small validation/smoke-test template, not a production sizing recommendation. Original workflow launchers are preserved as execution provenance and may contain historical paths/resource choices. Do not submit them unchanged merely because they are included.

Before any production inference or reanalysis: inventory complete/active/pending work, benchmark a representative allocation, preserve validated checkpoints, and group useful work to approximately 10–15 minutes per job where practical. Increasing walltime or concurrency does not correct short-job overhead. Never pad runtime with sleeps. Keep per-input failures visible and avoid duplicate writers. The prior Pairformer run timing is historical; its job sizing must not be adopted automatically for a new production stage.

Consult current NIH tool-specific documentation and available resources before selecting GPUs, memory and threads. No compute jobs were submitted during this repository packaging update.
