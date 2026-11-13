"""Resolve figure data locally and raw prediction inputs from the project root."""
from pathlib import Path
import os
FIGURE_DIR=Path(__file__).resolve().parents[1]
def project_root():
    if os.environ.get('LYDIA_PROJECT_ROOT'):
        return Path(os.environ['LYDIA_PROJECT_ROOT']).expanduser().resolve()
    for p in FIGURE_DIR.parents:
        if (p/'new_af3_runs/chimera_20261106').is_dir():
            return p
    raise FileNotFoundError('Raw ensemble not found. Set LYDIA_PROJECT_ROOT to the full project folder. Plotting/exporting panels does not require it.')
