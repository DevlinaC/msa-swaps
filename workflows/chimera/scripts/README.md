# Figure S9 scripts

`original/` contains byte-identical copies of the four scripts used for the original combined figure. They retain their historical directory assumptions; use the portable scripts in this folder for routine work.

## Commands for each panel and the combined figure

Run these commands locally from the project root:

```bash
cd /Users/chakravartyd2/Documents/AF_paper/Lydia_summer_student

# Panel A: construct schematic (PDF, PNG, SVG)
/opt/anaconda3/bin/python3 Manuscript/supporting_data/Figure_S9/scripts/export_panels.py --panel A

# Panel B: representative structures (PDF, PNG, SVG)
/opt/anaconda3/bin/python3 Manuscript/supporting_data/Figure_S9/scripts/export_panels.py --panel B

# Panel C: ensemble results (PDF, PNG, SVG), plus separate state/confidence plots (PNG, SVG)
/opt/anaconda3/bin/python3 Manuscript/supporting_data/Figure_S9/scripts/export_panels.py --panel C

# Assemble panels A–C into the combined Figure S9 (PDF, PNG, SVG)
/opt/anaconda3/bin/python3 Manuscript/supporting_data/Figure_S9/scripts/plot_chimera_s7.py
```

Individual panels are saved in `Panel_A_construct/`, `Panel_B_structures/`, and `Panel_C_results/`. The combined outputs are `Figure_S9.pdf`, `Figure_S9.png`, and `Figure_S9.svg` in the Figure S9 folder.

These commands read the canonical saved JSON data and structure renders. The combined command builds the complete layout directly from those sources; it does not stitch the exported panel files or incorporate manual edits made to those files. No analysis or structure rendering is rerun. To regenerate structure images first, use the optional rendering command below. On another computer, replace the project path and Python executable with your equivalents.

## Export the display items

```bash
python /path/to/Figure_S9/scripts/export_panels.py
```

Requires NumPy and Matplotlib. Reads the saved JSONs in `data/` and structure images in `renders/`; no raw predictions or PyMOL needed. Produces panel A, B and C PDF/PNG/SVG exports and the two separate C plots as PNG/SVG. It preserves the combined figure and analysis data.

`plot_chimera_s7.py` provides `build_figure()` for reuse. Running it directly regenerates the combined figure. This is a refactoring of the original drawing code, not a new analysis.

## Regenerate the structure images

```bash
python /path/to/Figure_S9/scripts/render_chimera_s7.py
python /path/to/Figure_S9/scripts/export_panels.py
```

Requires PyMOL (`pymol2`). Uses the two representative coordinate files in `structures/`. Writes the source images to `renders/` and the session/views to `editable/`. Copies inside `Panel_B_structures/` are packaging snapshots; refresh those copies if images or views are changed.

## Recalculate and validate the underlying analysis

```bash
export LYDIA_PROJECT_ROOT=/path/to/Lydia_summer_student
python /path/to/Figure_S9/scripts/analyze_chimera_s7.py
python /path/to/Figure_S9/scripts/validate_chimera_s7.py
```

Requires NumPy, Biopython, tmtools and PyMOL. The full project's raw inputs and all 50 prediction/confidence files must be present under their recorded paths, particularly `new_af3_runs/chimera_20261106/` and `fasta_files/RfaH.fasta`. `LYDIA_PROJECT_ROOT` is optional when this folder remains inside the full project. Those raw ensembles are not duplicated in the panel folders.

Analysis recalculation writes the canonical `data/`, representative structures and provenance; panel-specific data copies are snapshots and must be refreshed afterward. Validation checks raw coordinate hashes and independently saved confidence values. It does not perform visual inspection. No new AF2/AF3 inference is performed by these scripts.

Original analysis environment: Python with Biopython 1.78, NumPy 2.2.6, tmtools 0.1.1, Matplotlib 3.11.0 and PyMOL 2; exact captured versions are in `../provenance/validation.json`. Analysis/render scripts bound numerical/PyMOL threads to one. These commands are intended for local work; use an allocated compute node if recalculating on Biowulf.
