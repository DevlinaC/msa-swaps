# Pipeline and reproducibility scope

1. Curate 34 targets, sequences, reference chains and structural annotations (`data/publication/Table_S1`). S1 uses within-class structural comparisons; 181 unordered pairs yield 362 directed swaps.
2. Create native, swapped and single-sequence prediction inputs; retain corrected donor-2mm4 cases and resolved 60-residue inputs. Full AF inference datasets remain outside this compact repository.
3. Run AF2/ColabFold and AF3, with 25 predictions per input. Shared targets and model samples do not constitute independent biological replicates.
4. Score structures against query and donor experimental references. Current manuscript summaries use TM-score ≥0.6. Historical >0.75 tables in `data/results/` are labeled legacy; do not relabel them. The saved-score reconciliation command now supplies the complete model-labelled swapped table (9,050 predictions per method), both reference scores and explicit ≥0.6 summaries. Exact final Figure 2/3 plotting selections remain unconfirmed. Historical AF3 TM-align execution provenance still needs confirmation.
5. Infer GREMLIN, MSA Transformer and Pairformer contact patterns. GREMLIN's runner is provided. MSA Transformer production preprocessing is incompletely recovered. The corrected Pairformer package includes selected inputs, matrices, configuration and validation.
6. Calculate contact-set overlap, continuous-score correlation and MSA-to-structure agreement. See `contact_methods.md` for the differing definitions; S4 uses rank-1 predicted structures rather than experimental references.
7. Regenerate the available figures with `scripts/reproduce_figures.py`. The figure index separates table-to-plot reproduction, raw-data reanalysis and inference.

## Scientific source freeze

`provenance/import_manifest_20261113.json` identifies the manuscript/SI by SHA-256 and every imported source file by original project-relative path. Imported numerical tables are unchanged. New adapters relocate figure inputs/outputs. Active chimera scripts change figure numbering from S7 to S9; preserved original source snapshots remain available. The original manuscript and SI are not copied into this source-code repository.

The preserved October README/methods are under `docs/historical_20261030/`. Planning documents dated November 12 are historical inventories; consult `MISSING_DATA.md` for updates, particularly the now-completed 512-row Pairformer run.

## Outstanding reproduction limits

Main Figure 3's permutation test, S5's ensemble comparison, final case-study selections/calculations, GA95-L45Y outputs and S9 NMR sources are not yet reproduced by this package. The old AF3 TM-align parser has a known syntax defect and is not promoted as a runnable replacement for its historical execution. Exact upstream environments and original MSA Transformer sampling require further provenance. See the missing-source list for collaborator additions.
