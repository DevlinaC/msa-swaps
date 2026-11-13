# Data in this draft

This guide describes the **complete local package**. The [quick start](../README.md#quick-start-reproduce-packaged-figures) identifies plotting and collaborator-validation commands supported by the Git package. The full protected-source audit and structural-score reconstruction still require excluded records. Approved archive download links are not yet available; see the [upload policy](../docs/GIT_DATA_POLICY.md).

The [main Figure 2/3 plotting tables](publication/Figures_2_3/README.md) are included in full and can be plotted independently of the excluded bulk data. Their supplied selection and the Figure 3 statistical limitation are documented with the files.

`publication/S1` contains saved TM-score matrices, target metadata and cartoon renders. `publication/S2` contains GREMLIN and MSA Transformer comparison tables. `publication/S3` contains Pearson comparison tables and selected pairs. `publication/S4` contains primary and sensitivity comparisons, selection/mapping/input audits and selected-pair records. S4 plotting filters the primary top-L records explicitly. `publication/Table_S1` contains the current target table, sequences, coordinate files and CATH/ECOD annotation provenance. ECOD annotations are PDB-entry matches, not validated chain/domain assignments.

`workflows/pairformer/` contains all 396 selected input alignments and all 396 contact matrices (NumPy and equivalent CSV), alongside metadata and comparison tables. Matrices contain continuous predicted contact scores, not measured distances or experimental contacts. Use the saved donor mappings rather than comparing differently indexed matrices elementwise.

`workflows/chimera/` contains all 50 prediction summary records and representative model coordinates/renders. Full raw ensemble inputs and the NMR dataset are external.

The [Figure 5 package](publication/Figure_5/README.md) adds 100 AF2 mutant coordinate files, four A3Ms and confidence/provenance CSVs. The [case-study package](../workflows/case_studies/README.md) includes supplied Figure 4 and RfaH/NusG AF3 subsample inputs and reference metric tables. These new packages are eligible for Git in full; their commands do not require the excluded bulk publication inputs.

`targets.csv`, `results/`, `examples/` and `validation/` retain the October compact package unchanged locally to support historical overlap checks. These legacy paths are excluded from Git; use the complete local package to run the historical demo scripts. In particular, the >0.75 summaries and fixed-100 GREMLIN examples are not current manuscript summaries. Use `publication/Table_S1/table_S1.csv` for the updated annotations.

Import hashes and source paths are in `../provenance/import_manifest_20261113.json`. Absolute paths within historical metadata identify original execution locations; they are not required by the portable plotting commands. Full AF ensembles and the complete raw GREMLIN/MSA Transformer datasets will need a durable versioned archive and checksums before publication. No download URL has been invented for that archive.

Publication exports belong in `figures/`; the data folders retain scientific tables and required source assets. Superseded S3/S4 report images are in the ignored local archive. S1 individual structure renders and editable sessions remain here because figure generation uses them as source inputs.
