# Supplied plotting tables for main Figures 2 and 3A

These six CSVs were imported byte-for-byte from the November 13 collaborator delivery. The shared AF2 and AF3 swapped-score files were identical in the supplied Fig2 and Fig3 folders and are stored only once. [Source hashes and both original locations](../../../provenance/figures_2_3_sources.csv) are recorded. All six files are eligible for Git; this plotting package does not depend on the excluded bulk structural-score tables.

| Table | Records | Use |
| --- | ---: | --- |
| `AF2_TMs.csv` | 9,050 | AF2 swapped prediction scores against query (`AF2_A`) and donor (`AF2_B`); Figures 2 and 3A. |
| `tm_align_full_results.csv` | 18,950 | Saved AF3 reference comparisons; scripts select `Job_Type=junk` and match query/donor in `Job_Label`, yielding 9,050 scores per reference. |
| `jaccard_results_gremlin.csv` | 362 | Supplied historical GREMLIN ordinary Jaccard against query/donor; Figure 2. |
| `jaccard_results_msatransformer.csv` | 362 | Supplied historical MSA Transformer rank-tolerant overlap against query/donor; Figure 2. |
| `AF2_single_TMs.csv` | 850 | AF2 single-sequence query-reference scores; Figure 3A. |
| `singleseq_af3_tmscores.csv` | 850 | AF3 single-sequence scores with target, seed and sample identities; Figure 3A. |

## Reproduce the supplied plots

From the repository root, with Python 3.11 and `requirements.txt` installed:

```bash
python scripts/plot_main_figures.py --figure all
# Or choose --figure 2 or --figure 3.
```

The runner verifies source hashes and plotting populations, then executes the three original plotting scripts unchanged in a temporary directory with bounded threads and a noninteractive backend. Outputs go to `outputs/main_figures/`, with population summaries, output hashes, logs and environment information. Use a fresh `--out-dir` for repeat runs. Frozen reviewed exports are provided under [Figure 2](../../../figures/Figure_2/README.md) and [Figure 3](../../../figures/Figure_3/README.md). No inference, alignment or statistical test is run.

## Interpretation and pending decisions

**Figure 2 reproduces the supplied historical contact selection.** Its GREMLIN table is the original fixed-100 calculation, not the newer 100/100/150 S2 analysis. The MSA Transformer column name `Jaccard` denotes the historical rank-tolerant overlap. The supplied figure labels/style are unchanged: the shaded ≥0.6 region must not be read as a validated shared threshold across TM-score and both contact metrics. A revised main-panel contact selection and labeling remain to be reconciled with the manuscript; do not treat these exports as that revision.

**Figure 3A is a pooled-distribution plot, not the query-weighted statistical test.** The supplied AF2 CSVs lack query/model labels. Their pooled swapped query/donor and single-sequence query distributions match the labelled AF2 raw-log audit at five-decimal precision, checked during integration. Do not assign AF2 labels using AF3 row order: that is the documented error in the recovered statistical script. The `ks_cluster_analysis.py` source remains preserved for review but is never executed by this runner. Final corrected D/p reporting remains pending.

The saved AF3 normalization/execution provenance remains subject to the existing [structural-score limitations](../structure_scores/README.md). This package enables faithful plotting of supplied values without claiming a fresh TM-align calculation. Original source images were visually compared with regenerated panels; small renderer-dependent pixel/dimension differences were observed, not byte-identical image reproduction.
