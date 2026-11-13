# Panel C: ensemble results

`Figure_S9C.pdf`, `.svg` and 350-dpi `.png` contain both results plots. `CTD_structural_state.svg/.png` and `CTD_confidence.svg/.png` provide the two plots separately for independent placement.

`data/models.json` contains all 50 per-model records, including CTD confidence, secondary-structure assignments, ranking and coordinate-source hashes. `data/summary.json` contains method-level summaries. These are identical copies of the canonical parent `data/` files at packaging time.

All 25 predictions per method have beta-sheet CTDs. Mean CTD pLDDT is >=70 in 18/25 AF2 and 25/25 AF3 predictions. Confidence points retain all models, with median lines and a threshold line at 70.

Reproduce displays using `../scripts/export_panels.py`. The current drawing preserves the original state counts and asserts the observed strand/helix criteria; it is intended for this audited dataset, not arbitrary replacement models. Recalculating raw-model measurements uses `../scripts/analyze_chimera_s7.py` and requires the full prediction ensemble.
