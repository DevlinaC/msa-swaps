"""Render actual rank-one chimera coordinates with shared CTD alignment/view."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pymol2,json
from _paths import FIGURE_DIR
O=FIGURE_DIR
(O/'editable').mkdir(exist_ok=True)
(O/'renders').mkdir(exist_ok=True)
with pymol2.PyMOL() as pm:
 c=pm.cmd;c.viewport(1200,1050);c.set('max_threads',1);c.bg_color('white');c.set('orthoscopic',1);c.set('ray_shadows',0);c.set('specular',.20);c.set('antialias',2);c.set('ray_opaque_background',1)
 c.set_color('burgundy',[172/255,0,53/255]);c.set_color('ntd_grey',[.68,.68,.68])
 for method in ['AF2','AF3']:
  p=next((O/'structures').glob(method+'_rank1.*'));c.load(str(p),method);c.dss(method)
  c.hide('everything',method);c.show('cartoon',method);c.color('ntd_grey',method);c.color('burgundy',f'{method} and resi 131-183')
 c.align('AF3 and resi 131-183 and name CA','AF2 and resi 131-183 and name CA',cycles=0)
 c.disable('all');c.enable('AF2');c.orient('AF2');c.zoom('AF2',buffer=9);view=c.get_view()
 for method in ['AF2','AF3']:
  c.disable('all');c.enable(method);c.set_view(view);c.zoom(method,buffer=9)
  c.png(str(O/'renders'/f'{method}_full.png'),width=1200,height=1050,dpi=300,ray=1)
  c.scene(method+'_full','store')
 c.disable('all');c.enable('AF2');c.hide('cartoon','resi 1-130');c.orient('AF2 and resi 131-183');c.zoom('AF2 and resi 131-183',buffer=3);ctdview=c.get_view()
 for method in ['AF2','AF3']:
  c.disable('all');c.enable(method);c.set_view(ctdview)
  c.png(str(O/'renders'/f'{method}_CTD.png'),width=900,height=750,dpi=300,ray=1);c.scene(method+'_CTD','store')
 c.show('cartoon','all');c.disable('all');c.enable('AF2');c.set_view(view)
 c.save(str(O/'editable/Figure_S9_structures.pse'))
 (O/'editable/views.json').write_text(json.dumps({'full':list(view),'ctd':list(ctdview)},indent=2)+'\n')
print('Rendered rank-one structures and CTD views; saved editable session.')
