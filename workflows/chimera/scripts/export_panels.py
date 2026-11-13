"""Export individual Figure S9 panels without changing combined outputs or data.

Run from any directory: python /path/to/Figure_S9/scripts/export_panels.py
Requires numpy and matplotlib; no raw ensemble or PyMOL required.
"""
import argparse
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/lydia_report_mpl")
from pathlib import Path
import json,hashlib,shutil
from datetime import datetime, timezone
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--panel", choices=["A", "B", "C", "all"], default="all",
                    help="Panel to export (default: all); C also exports its two individual plots.")
args = parser.parse_args()
from matplotlib.transforms import Bbox
from plot_chimera_s7 import build_figure,plt
from _paths import FIGURE_DIR as O
fig=build_figure();w,h=fig.get_size_inches()
parts=[
 ('Panel_A_construct/Figure_S9A',(0.025,.825,.97,1.0),['pdf','png','svg']),
 ('Panel_B_structures/Figure_S9B',(.025,.375,.99,.815),['pdf','png','svg']),
 ('Panel_C_results/Figure_S9C',(.025,.005,.99,.375),['pdf','png','svg']),
]
outputs=[]
for stem,(x0,y0,x1,y1),formats in parts:
 if args.panel != "all" and not stem.endswith(args.panel):
  continue
 for ext in formats:
  p=O/f'{stem}.{ext}';fig.savefig(p,dpi=350,bbox_inches=Bbox.from_extents(x0*w,y0*h,x1*w,y1*h));outputs.append(p)
# Split the two C plots for independent placement; no duplicate panel letter.
for t in fig.texts:
 if t.get_text()=='C':t.set_visible(False)
for name,bounds in [('CTD_structural_state',(.025,.005,.50,.375)),('CTD_confidence',(.53,.005,.995,.375))]:
 if args.panel not in ['C','all']:
  continue
 for ext in ['png','svg']:
  p=O/'Panel_C_results'/f'{name}.{ext}';x0,y0,x1,y1=bounds;fig.savefig(p,dpi=350,bbox_inches=Bbox.from_extents(x0*w,y0*h,x1*w,y1*h));outputs.append(p)
plt.close(fig)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
inputs=[O/'data/models.json',O/'data/summary.json',O/'renders/AF2_full.png',O/'renders/AF3_full.png']
manifest = 'panel_exports.json' if args.panel == 'all' else f'panel_{args.panel}_exports.json'
(O/'provenance').mkdir(parents=True, exist_ok=True)
(O/'provenance'/manifest).write_text(json.dumps({'created':datetime.now(timezone.utc).isoformat(),'inputs':{str(p.relative_to(O)):sha(p) for p in inputs},'outputs':{str(p.relative_to(O)):sha(p) for p in outputs},'scripts':{p.name:sha(p) for p in Path(__file__).parent.glob('*.py')},'note':'Same layout, values, colors and structure images as original plot; export bounds isolate each display. Original combined files are unchanged.'},indent=2)+'\n')
print(f'Exported {len(outputs)} files for panel selection: {args.panel}.')
