# Collaborator notes and repository status

Updated November 13, 2026. This document contains the preparation and release notes previously shown in the main README. The main README introduces the paper and provides instructions for readers.

## Manuscript version and coverage

The collaboration draft follows the locally available November 12 manuscript and SI, with the corrected November 13 Pairformer run. Some final case-study selections, statistical code and experimental data are awaiting collaborator files. This is not yet a complete publication release.

- [Publication metadata and pending release items](PUBLICATION_METADATA.md)
- [Figure coverage and exact commands](FIGURES.md)
- [Missing sources and decisions](MISSING_DATA.md)
- [Instructions for contributing files](../CONTRIBUTING.md)
- [Upload instructions and temporary exclusions](PUBLISHING.md)
- [November 13 update record](UPDATE_20261113.md)

## README artwork

The README now uses a preview recovered from Figure 1 in the current manuscript, with document and embedded-image checksums. Replace it with the final author-supplied editable/high-resolution export when available. The old `MSA_Swap_story.png` remains historical artwork and is no longer embedded.

## Current analysis scope

The available plotting commands regenerate S1–S4, the corrected four-metric Pairformer analysis intended for S6, and S9 prediction panels A–C. They do not regenerate the unavailable main-figure sources, S5, S7, S8 or S9 NMR panels. Compare regenerated figures against the authors' final artwork before release. The current SI still contains the earlier three-panel Pairformer artwork.

The corrected Pairformer package contains all 396 selected MSAs and matrices, uses at most 512 sequences, and disables query weighting. The older 256-row/bias-enabled results are preserved in the parent project as historical results; they have not been relabeled.

The complete model-labelled swapped-score reconciliation contains 9,050 predictions per AF method and both reference TM-scores. Neither method has a saved donor score ≥0.6. This reconciles saved AF2 raw-log records and AF3 score-to-model links; no new TM-align calculation was performed. Historical AF3 execution/reference provenance still needs confirmation.

## Inference and raw-data limitations

| Workflow | Remaining dependency or source question |
| --- | --- |
| GREMLIN | Reconcile the exact original inference environment and runner version. |
| MSA Transformer | Recover production `modules.MSA_TR_af_sampling` and historical subsets. The first-64-row test implementation is not a substitute. |
| MSA Pairformer | Fresh inference requires external upstream code/weights and the pinned pilot environment/control artifacts; table/matrix analysis does not. |
| Structural contacts | Raw coordinate reanalysis requires the external full project data layout. |
| Chimera | Full ensemble reanalysis requires the external 50-model dataset; NMR source data remain outstanding. |

## Local files versus upload contents

The temporary `.gitignore` excludes local environments, JSON/JSONL/JSON5 files except the [15 reviewed JSON exceptions](../provenance/README.md#json-files-eligible-for-git), full contact-matrix CSVs, compressed CSVs and explicitly listed CSVs ≥1 MiB. Small summary/metadata CSVs remain eligible. Files remain locally available. Git cannot filter by file size, so review the explicit path list as data are added.

These exclusions include scientific inputs, configuration and provenance required by the reproduction commands. The complete local package works; a GitHub upload using these exclusions will need the approved files restored or supplied through a versioned archive. Decide the public sharing scope, add archive links and verify a fresh checkout before calling it a complete reproducibility release.

## Historical package and release decisions

The October demo commands (`validate_examples.py`, `compute_contacts.py`, `plot_results.py`) retain their fixed-100/legacy metrics and >0.75 summaries. They are not the current S2–S4 pipeline; original documentation and the three archived plots are in `historical_20261030/` (local-only; excluded from Git). Regenerating those plots writes to `outputs/historical_20261030/plots/`.

Original project code is now covered by the [MIT License](../LICENSE), with an AI-assistance disclosure in the README. See [licensing scope](LICENSING.md) for the distinction between original code, upstream-derived code and non-software assets. Data/asset release terms, the manuscript DOI and the permanent data archive URL remain to be specified. Preserve upstream attribution and licensing. Suggested initial destination is `DevlinaC/msa-swaps`, private during co-author revision. No remote repository or collaborator invitation has been created, and no email has been sent.
