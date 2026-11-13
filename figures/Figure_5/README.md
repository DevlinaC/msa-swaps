# Figure 5 mutant model package

The [reviewed data package](../../data/publication/Figure_5/README.md) contains 100 AF2 models, their four A3M inputs, per-model and per-residue confidence records, and original run configurations.

| Construct | Models | Residues per model | Confidence values |
| --- | ---: | ---: | ---: |
| GB98 T25V-I17T | 25 | 56 | 1,400 |
| GB98 T25E-K28I | 25 | 56 | 1,400 |
| GB98 T25E-W43E | 25 | 56 | 1,400 |
| RfaH L142K-H152V | 25 | 162 | 4,050 |

Each ensemble covers five AF2 networks and five seeds. [Models and identities](../../data/publication/Figure_5/models.csv), [input identities](../../data/publication/Figure_5/inputs.csv), [confidence values](../../data/publication/Figure_5/residue_plddt.csv) and [configurations](../../data/publication/Figure_5/run_config.csv) are available as CSVs. Original PDBs and A3Ms are retained unchanged in the data package; no duplicate model files are needed here.

The [source review](../../provenance/figure5_source_review_20261113.csv) confirms exact agreement with the original archives for all 100 PDBs, four A3Ms, 8,250 score values and all 32 configuration fields in each of the two archive configurations. Original score JSONs contain only pLDDT arrays, fully retained in the CSV. Archive and model/score checksums are recorded in the data package.

To validate the packaged data from the repository root:

```bash
python scripts/validate_collaborator_import.py
```

**These data do not yet reproduce the final Figure 5 fold percentages.** Model-level fold assignments, analyzed regions, confidence inclusion rules and denominators remain unresolved. The 20 feature-selected ensembles and AF3 counterpart classifications are not supplied by these four AF2 mutant ensembles. The supplied bar-plot script uses hard-coded percentages and is not a model-classification analysis; no percentage figures are presented here as validated results.
