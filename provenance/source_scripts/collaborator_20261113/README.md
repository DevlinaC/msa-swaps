# Preserved Figure 2 and Figure 3 source scripts

These scripts are unchanged copies from the November 13 collaborator delivery, retained for review and traceability. Their hashes are in [the import manifest](../../collaborator_import_20261113.csv). The three plotting scripts are now run unchanged by `scripts/plot_main_figures.py` using the [packaged Figure 2/3 tables](../../../data/publication/Figures_2_3/README.md). The statistical script remains review-only; it references external paths and has the grouping error described below.

- Figure 2 supplied overlap tables duplicate the old fixed-100 GREMLIN analysis. Reconcile the main panel with the agreed revised selection rule before promoting it.
- Figure 3 `ks_cluster_analysis.py` assigns query labels to unlabelled AF2 vectors using AF3 block order. Independent matching to labelled AF2 records confirms only 30/362 swapped and 2/34 single-sequence blocks are correctly assigned. Preserve the supplied implementation but use explicit AF2 identities in a corrected active workflow.
- The audit reports corrected AF2 D=0.1624434, p=0.00019999, compared with supplied D=0.1542383, p=0.0248988 (20,000 permutations, seed 20261029). This is a pending manuscript correction, not a newly adopted manuscript result. Pooled violin distributions are unaffected by relabelling; AF3 is unchanged.

Original audit evidence remains in the parent project's `reports/supervisor_files_audit_20261113/`. Supplied Figure 2 and Figure 3A panels are now reproducible. Their plotting does not update any manuscript statistical claim or resolve Figure 2 contact-selection differences.
