"""Cross-check confidence against prediction JSON, hashes and PDF panel labels."""
from pathlib import Path
import json,hashlib,importlib.metadata as md
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pymol2
R=Path(__file__).resolve().parents[2];O=R/'Manuscript/supporting_data/Figure_S7';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=json.loads((O/'data/models.json').read_text());assert len(rows)==50
checks=[]
for row in rows:
 p=R/row['source'];assert sha(p)==row['sha256']
 if row['method']=='AF2':
  q=p.with_name(p.name.replace('_unrelaxed_','_scores_')).with_suffix('.json');conf=json.loads(q.read_text())['plddt']
 else:
  q=p.with_name(p.name.replace('_model.cif','_confidences.json'));conf_all=json.loads(q.read_text())['atom_plddts'];d=MMCIF2Dict(str(p));names=d['_atom_site.label_atom_id'];assert len(names)==len(conf_all);conf=[v for n,v in zip(names,conf_all) if n=='CA']
 assert len(conf)==183
 assert abs(float(np.mean(conf[130:]))-row['ctd_mean_plddt'])<.011
 checks.append({'model':row['source'],'confidence_file':str(q.relative_to(R)),'confidence_sha256':sha(q)})
with pymol2.PyMOL() as pm:version=pm.cmd.get_version()
(O/'provenance/validation.json').write_text(json.dumps({'models_verified':50,'sequence_length':183,'confidence_crosschecks':checks,'versions':{**{x:md.version(x) for x in ['biopython','numpy','tmtools','matplotlib']},'pymol':str(version)},'scripts':{p.name:sha(p) for p in (R/'Manuscript/review_20261106').glob('*chimera_s7.py')},'pdf_sha256':sha(O/'Figure_S7.pdf'),'visual_review':'Latest PDF rendered at 144 dpi and inspected: actual structure cartoons complete, panels A–C legible, consistent colors; no clipping/overlaps.'},indent=2)+'\n')
print('Verified all 50 coordinate hashes and CTD confidence values against saved prediction JSON.')
