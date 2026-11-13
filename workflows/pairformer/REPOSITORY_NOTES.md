# Pairformer package in this repository

This is the complete locally staged November 13, 2026 **512-row, query-weighting-disabled** package, excluding transient Python caches. All 396 selected inputs and matrices, scientific tables, original scripts, run provenance and inventory are included. Historical absolute paths remain in provenance. The package README's sibling `../msa_pairformer_full/` refers to the external historical project folder, not a dependency for plotting or matrix analysis in this repository.

From the repository root:

```bash
python workflows/pairformer/scripts/plot_summary.py --root workflows/pairformer
```

Raw matrix reanalysis is documented in `../../docs/FIGURES.md`. For inference, the runner still requires the pinned upstream package/weights and the pilot environment/control artifacts referenced by its `--pilot` argument. Those weights/environments are not redistributed. Do not treat plotting success as validation of a fresh inference environment. Preparation from original unsampled MSAs also requires the external source inventory paths.

The four-panel result is the proposed updated S6 source; the locally available SI contains the earlier three-panel artwork and has not been rewritten by this repository update.
