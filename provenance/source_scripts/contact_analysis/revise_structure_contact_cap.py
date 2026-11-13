"""Item 3: raise only class-150 secondary GREMLIN cap, retaining z-score >=1.

Run on a compute node. Reuse validated distances for the already exported
top-150 pairs; preserve every primary result and avoid duplicate structure work.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import csv
import gzip
import hashlib
import json
import platform
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();src,out=a.source,a.out
    assert os.getenv('SLURM_JOB_ID'), 'Run the revision on a Slurm compute node'
    assert not out.exists()
    pro=json.loads((src/'provenance.json').read_text())
    assert sha(src/'provenance.json')==json.loads((src/'COMPLETE.json').read_text())['provenance_sha256']
    for name,digest in pro['outputs'].items():assert sha(src/name)==digest,name
    frame=pd.read_csv(src/'comparisons.csv');original=frame.copy(deep=True)
    assert len(frame)==20160
    cols=['name','af','method','reference','separation']
    primary=frame[(frame.method=='GREMLIN')&(frame.length_class==150)&(frame.selection=='topL')]
    secondary=frame[(frame.method=='GREMLIN')&(frame.length_class==150)&(frame.selection=='top100_z1')]
    assert len(primary)==len(secondary)==1680
    assert primary.L.eq(150).all() and primary.n_selected.eq(150).all()
    primary_ids=set(primary.record_id);members=defaultdict(list)
    with gzip.open(src/'selected_pairs.csv.gz','rt') as f:
        reader=csv.DictReader(f);fields=reader.fieldnames
        for row in reader:
            if int(row['record_id']) in primary_ids:members[int(row['record_id'])].append(row)
    lookup={tuple(getattr(r,k) for k in cols):r.record_id for r in primary.itertuples(index=False)}
    replacements={};checks=[]
    for row in secondary.itertuples(index=False):
        pid=lookup[tuple(getattr(row,k) for k in cols)]
        candidates=members[pid];assert len(candidates)==150
        ranked=sorted([r for r in candidates if float(r['score'])>=1],
                      key=lambda r:(-float(r['score']),int(r['i']),int(r['j'])))
        # Independent regression: the previously selected top 100 must be
        # recoverable from these same top-150 pairs before extending the cap.
        previous=ranked[:100]
        hits=sum(float(r['distance_A'])<=8 for r in previous)
        assert len(previous)==row.n_selected and hits==row.intersection,(row.record_id,'baseline')
        assert abs(min(float(r['score']) for r in previous)-row.cutoff)<1e-12
        chosen=ranked[:150]
        # Export in coordinate order, matching the original pair-list convention.
        chosen=sorted(chosen,key=lambda r:(int(r['i']),int(r['j'])))
        replacements[row.record_id]=[dict(r,record_id=str(row.record_id)) for r in chosen]
        n=len(chosen);tp=sum(float(r['distance_A'])<=8 for r in chosen)
        union=n+row.n_structural-tp
        updates=dict(n_selected=n,intersection=tp,union=union,
                     cutoff=min(float(r['score']) for r in chosen),precision=tp/n,
                     recall=tp/row.n_structural if row.n_structural else np.nan,jaccard=tp/union)
        index=frame.index[frame.record_id==row.record_id];assert len(index)==1
        for name,value in updates.items():frame.loc[index,name]=value
        checks.append(dict(record_id=row.record_id,old_selected=row.n_selected,new_selected=n,
                           old_precision=row.precision,new_precision=tp/n,baseline_reproduced=True))
    # Every result outside this specific secondary subgroup must be untouched.
    unchanged=~frame.record_id.isin(replacements)
    pd.testing.assert_frame_equal(frame.loc[unchanged],original.loc[unchanged])
    frame.loc[frame.selection=='top100_z1','selection']='class_cap_z1'
    out.mkdir(parents=True)
    frame.to_csv(out/'comparisons.csv',index=False)
    pd.DataFrame(checks).to_csv(out/'cap_revision_audit.csv',index=False)
    with gzip.open(src/'selected_pairs.csv.gz','rt') as fin,gzip.open(out/'selected_pairs.csv.gz','wt',newline='') as fout:
        reader=csv.DictReader(fin);writer=csv.DictWriter(fout,fieldnames=fields);writer.writeheader();emitted=set()
        for row in reader:
            rid=int(row['record_id'])
            if rid in replacements:
                if rid not in emitted:writer.writerows(replacements[rid]);emitted.add(rid)
            else:writer.writerow(row)
        assert emitted==set(replacements)
    for name in ('structure_audit.csv','prediction_input_audit.csv','mapping_audit.csv'):
        shutil.copy2(src/name,out/name)
    summaries=[]
    for keys,sub in frame.groupby(['af','method','selection','reference','separation']):
        r=dict(zip(['af','method','selection','reference','separation'],keys))
        r.update(n=len(sub),n_defined=int(sub.precision.notna().sum()))
        for metric in ('precision','recall','jaccard'):
            r[metric+'_mean']=sub[metric].mean();r[metric+'_median']=sub[metric].median()
        summaries.append(r)
    pd.DataFrame(summaries).to_csv(out/'summary.csv',index=False)
    provenance=dict(created_utc=datetime.now(timezone.utc).isoformat(),host=platform.node(),job_id=os.getenv('SLURM_JOB_ID'),
        source_directory=str(src),source_provenance_sha256=sha(src/'provenance.json'),
        script_sha256=sha(Path(__file__)),parameters=dict(primary='unchanged top L with cutoff ties',
            secondary='class_cap_z1: cap 150 in class 150; cap 100 in classes 60 and 100',secondary_min_zscore=1.0,
            secondary_ties='deterministic descending score then i then j, unchanged'),
        change_note='Original supplied GREMLIN Jaccard scripts set GREMLIN_TOP_N=100. '
                    'Only the secondary class-150 cap is raised here; primary results already selected 150 contacts.',
        validation=dict(secondary_records_revised=1680,old_top100_records_reproduced=1680,
            unchanged_result_records=int(unchanged.sum()),new_selected_min=min(x['new_selected'] for x in checks),
            new_selected_max=max(x['new_selected'] for x in checks)),
        outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()})
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (out/'COMPLETE.json').write_text(json.dumps(dict(provenance_sha256=sha(out/'provenance.json')),indent=2)+'\n')
    print(json.dumps(provenance['validation']),flush=True)


if __name__=='__main__':main()
