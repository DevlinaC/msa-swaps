# MSA Transformer provenance

The underlying pretrained predictor is ESM's `esm_msa1b_t12_100M_UR50S`, accessed through `esm.pretrained`. `MSA_TR` is a project helper function, not the upstream model name.

The available `historical/msatr_af_sampling.py` wrapper imports `MSA_TR` from `modules.MSA_TR_af_sampling` and passes `nseq=256` by default. That production module has not yet been recovered in the local source files. The local `historical/MSA_TR_test.py` defines a different `MSA_TR`: it takes the first 64 MSA rows, runs on CPU and has no `nseq` parameter. Renaming that test file would not reproduce the historical sampling, and calling it with `nseq` would fail.

Saved matrices and downstream Jaccard/rank-tolerant code can be shared and analyzed independently of this missing inference preprocessing. Do not present a newly written wrapper or the first-64-row test as the historical implementation. Recover the production module, model/package version, row-selection method, selected inputs/indices and seed before claiming end-to-end reproduction.
