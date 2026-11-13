# Panel B: representative structures

`Figure_S9B.pdf`, `.svg` and 350-dpi `.png` contain the two labeled structure views with whole-model and CTD pLDDT values. Individual structure images are in `renders/` (full-model and CTD-only views for AF2 and AF3). `editable/` contains the PyMOL session and saved views; the molecular images embedded in SVG/PDF are raster renders, not vector molecular drawings.

`data/representative_models.json` contains the two rank-one model records. Coordinates remain in `../structures/`. Representatives were selected using the original rankings, independently of CTD state.

Reproduce images with `../scripts/render_chimera_s7.py`, then exports with `../scripts/export_panels.py`. Panel-specific renders/session/data are snapshots of the canonical files in the parent folder.
