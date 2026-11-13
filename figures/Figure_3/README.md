# Figure 3A pooled structural-score distributions

These PNG/PDF/SVG exports reproduce the supplied Figure 3A plotting script using the [packaged source tables](../../data/publication/Figures_2_3/README.md). There are 9,050 swapped and 850 single-sequence query-reference scores per AF method.

Run `python scripts/plot_main_figures.py --figure 3` from the repository root; outputs go to ignored `outputs/main_figures/`. This command does not run the query-weighted permutation test. Correcting the documented AF2 grouping error affects the statistical calculation, not these pooled distributions. No significance annotation has been added.
