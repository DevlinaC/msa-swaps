#!/usr/bin/env python3
"""Paired comparison of preserved 256/biased and manuscript 512/unweighted runs."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,csv,io,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent/'pilot_lib'))
from common import sha256,require,atomic_bytes,write_json

def table(p):
    with p.open() as f:return list(csv.DictReader(f))

def save(p,rows):
    b=io.StringIO();w=csv.DictWriter(b,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows);atomic_bytes(p,b.getvalue().encode())

def compare(root,old):
    out=root/'results/configuration_comparison';out.mkdir(exist_ok=True)
    a=json.loads((old/'inputs/manifest.json').read_text());b=json.loads((root/'inputs/manifest.json').read_text())
    prior={s['id']:s for s in a['inputs']};current={s['id']:s for s in b['inputs']}
    require(prior.keys()==current.keys(),'Input identities differ')
    source_hashes={}
    for project in (old,root):
        for rel in ('inputs/manifest.json','results/analysis/VALIDATED.json','results/RUN.json','results/analysis/comparisons_with_rank_tolerance.csv'):
            source_hashes[str(project/rel)]=sha256(project/rel)
    maps=[]
    for name,s in current.items():
        t=prior[name]
        require(s['source_sha256']==t['source_sha256'] and s['donor_mapping']==t['donor_mapping'],'Source or mapping changed')
        require(set(t['selected_row_indices']).issubset(s['selected_row_indices']),'Old selected rows not nested')
        x=np.load(old/'results'/name/'contacts.npy',allow_pickle=False);y=np.load(root/'results'/name/'contacts.npy',allow_pickle=False)
        maps.append(dict(id=name,kind=s['kind'],query=s['query'],donor=s['donor'],old_depth=t['selected_depth'],new_depth=s['selected_depth'],selection_changed=t['selected_sha256']!=s['selected_sha256'],old_npy_sha256=sha256(old/'results'/name/'contacts.npy'),new_npy_sha256=sha256(root/'results'/name/'contacts.npy'),max_abs_change=float(abs(x-y).max()),mean_abs_change=float(abs(x-y).mean())))
    save(out/'matrix_changes.csv',maps)
    key=lambda r:(r['id'],r['separation'],r['reference'])
    oldrows={key(r):r for r in table(old/'results/analysis/comparisons_with_rank_tolerance.csv')}
    newrows=table(root/'results/analysis/comparisons_with_rank_tolerance.csv')
    require(len(oldrows)==len(newrows)==2172,'Incomplete comparison table')
    metrics=('ordinary_jaccard','tolerant_jaccard','pearson_all_eligible','pearson_union_top_L');changes=[];summary=[]
    for r in newrows:
        t=oldrows[key(r)];require(r['L']==t['L'],'Mapped length changed')
        for metric in metrics:
            x,y=float(t[metric]),float(r[metric]);changes.append(dict(id=r['id'],separation=r['separation'],reference=r['reference'],metric=metric,old=x,new=y,delta=y-x))
    save(out/'metric_changes.csv',changes)
    for sep in ('6','12','24'):
        for ref in ('query','donor'):
            for metric in metrics:
                group=[r for r in changes if r['separation']==sep and r['reference']==ref and r['metric']==metric]
                summary.append(dict(separation=sep,reference=ref,metric=metric,n=len(group),old_mean=float(np.mean([r['old'] for r in group])),new_mean=float(np.mean([r['new'] for r in group])),mean_delta=float(np.mean([r['delta'] for r in group]))))
    save(out/'summary.csv',summary)
    write_json(out/'VALIDATED.json',dict(status='PASS',old_configuration='max256/query_bias_on',new_configuration='max512/query_bias_off',matrices=396,comparison_rows=2172,source_and_mapping_checks=396,nested_selection_checks=396,depth_increased=sum(r['new_depth']>r['old_depth'] for r in maps),matrices_changed=sum(r['max_abs_change']>0 for r in maps),source_hashes=source_hashes,code_sha256=sha256(Path(__file__)),caution='Two settings changed together; differences do not separately identify effects of depth and query weighting.',summary=summary))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent);p.add_argument('--old',type=Path,required=True);args=p.parse_args();compare(args.root.resolve(),args.old.resolve())
