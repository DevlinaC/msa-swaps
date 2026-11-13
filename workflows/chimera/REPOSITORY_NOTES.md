# Chimera numbering and reproducibility

This package is **current Figure S9 A–C**. It was copied from the original local `Figure_S7` folder. Original file names and source snapshots are retained for traceability; active plotting/export scripts now produce `Figure_S9.*` and individual S9A/B/C files. Script `plot_chimera_s7.py` retains its historical filename.

Run `python scripts/reproduce_figures.py --figure S9` from the repository root. This recreates prediction panels from saved model summaries and renders. It does not create NMR panels D/E. Full ensemble reanalysis requires the external raw chimera dataset via `LYDIA_PROJECT_ROOT`, as documented by the original `_paths.py`. No AF inference was rerun here. The import manifest records original source hashes; active numbering edits intentionally differ from those source hashes.
