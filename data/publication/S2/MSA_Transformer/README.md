# MSA Transformer Jaccard across sequence separations

This extends the manuscript Figure 2B rank-tolerant Jaccard from |i−j| ≥6 to
≥12 and ≥24 using the same saved matrices, center-cropping, top-L/top-⌊1.5L⌋
definition, and cutoff-tie handling. All 362 swaps and both query/donor
references are included (2,172 records).

See `../../../reports/figure2_msatransformer_jaccard/README.md` for the caption,
summary, validation, and discussion of the separate GREMLIN top-100 cutoff.
`provenance.json` records input and output hashes; `COMPLETE.json` authenticates
that provenance. `baseline_validation.csv` checks all 724 original ≥6 values.
`summary.csv` includes pooled and length-class summaries; `comparisons.csv`
contains per-swap scores and selection counts, including strict Jaccard as a
diagnostic. The displayed figure uses the manuscript's rank-tolerant metric.
