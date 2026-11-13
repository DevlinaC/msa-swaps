#!/usr/bin/env python3
"""Validate frozen source bytes, analysis populations and current inference settings."""
import hashlib
import csv
import json
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def require(ok,message):
    if not ok: raise ValueError(message)
def main():
    # Original source digests remain authoritative after local archival moves.
    manifest=ROOT/'provenance/archive_manifest_20261113.csv'
    archived={}
    if manifest.is_file():
        with manifest.open(newline='') as handle:
            archived={row['original_path']: row['archive_path'] for row in csv.DictReader(handle)}
    frozen=json.loads((ROOT/'provenance/protected_data_sha256.json').read_text())
    for rel,digest in frozen.items():
        path=ROOT/rel
        if not path.is_file() and rel in archived:
            path=ROOT/archived[rel]
        require(path.is_file(),'Missing '+rel)
        require(hashlib.sha256(path.read_bytes()).hexdigest()==digest,'Changed scientific source '+rel)
    data=ROOT/'data/publication'
    targets=pd.read_csv(data/'Table_S1/table_S1.csv')
    require(len(targets)==34 and targets.pdb_id.nunique()==34,'Target IDs/count')
    require(targets.groupby('length_class').size().to_dict()=={60:14,100:10,150:10},'Target classes')
    for method,metric in [('GREMLIN','jaccard'),('MSA_Transformer','tolerant_jaccard')]:
        frame=pd.read_csv(data/'S2'/method/'comparisons.csv')
        selected=frame[frame.separation.isin([6,12,24])]
        require(len(selected)==2172,'S2 count '+method)
        require((selected.groupby(['separation','reference']).size()==362).all(),'S2 populations')
        require(selected[metric].between(0,1).all(),'S2 score range')
    frame=pd.read_csv(data/'S3/comparisons.csv')
    require(len(frame)==4344,'S3 count')
    require((frame.groupby(['method','separation','reference']).size()==362).all(),'S3 populations')
    require(frame.pearson_r.dropna().between(-1,1).all(),'S3 score range')
    frame=pd.read_csv(data/'S4/comparisons.csv')
    primary=frame[frame.selection=='topL']
    counts=primary.groupby(['af','method','separation','reference']).size()
    require(len(counts)==48,'S4 group count')
    for key,count in counts.items():
        require(count==(34 if key[-1]=='native_baseline' else 362),'S4 population '+str(key))
    for metric in ('precision','jaccard'):
        require(primary[metric].between(0,1).all(),'S4 score range '+metric)
    root=ROOT/'workflows/pairformer'
    config=json.loads((root/'model_config.json').read_text())
    require(config['max_depth']==512 and config['query_biasing'] is False,'Pairformer configuration')
    complete=json.loads((root/'results/INFERENCE_COMPLETE.json').read_text())
    require(complete['validated']==396 and complete['reused_pilot']==0,'Pairformer completion')
    matrices=list((root/'results').glob('*/contacts.npy'))
    require(len(matrices)==396,'Pairformer matrix count')
    for p in matrices:
        a=np.load(p,allow_pickle=False)
        require(a.ndim==2 and a.shape[0]==a.shape[1] and np.isfinite(a).all(),'Invalid matrix '+str(p))
        require(a.min()>=0 and a.max()<=1 and np.allclose(a,a.T,atol=1e-5),'Matrix symmetry/range '+str(p))
    a=pd.read_csv(root/'results/analysis/comparisons.csv')
    b=pd.read_csv(root/'results/rank_tolerant/comparisons.csv')
    keys=['id','separation','reference']
    require(len(a)==len(b)==2172 and not a.duplicated(keys).any() and not b.duplicated(keys).any(),'Pairformer keys')
    require(set(map(tuple,a[keys].values))==set(map(tuple,b[keys].values)),'Pairformer comparison correspondence')
    models=json.loads((ROOT/'workflows/chimera/data/models.json').read_text())
    require(len(models)==50 and sum(x['method']=='AF2' for x in models)==25 and sum(x['method']=='AF3' for x in models)==25,'Chimera counts')
    print(f'PASS: {len(frozen)} protected files; 34 targets; S2/S3/S4 populations; 396 current Pairformer matrices; 50 chimera models.')
if __name__=='__main__':main()
