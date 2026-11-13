#!/usr/bin/env python3
"""Item 1 revision: cap 150 for the 150-residue class; cap 100 otherwise.

Preserve z-score >=1, gap filtering, historical sorting/ties and offsets.
Only the cap for the 150-residue class changes.
The original scripts and outputs are preserved and independently replayed.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/lydia_report_mpl')
import gzip
import csv
import json
import platform
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import analyze_gremlin_separation as g

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'for_github_repo/analysis/gremlin_cap150_separation'
SEPS=(6,12,24)


def select(frame,mask,sep,cap):
    ij=frame[['i','j']].to_numpy(dtype=int)
    eligible=frame.loc[mask[ij[:,0]] & mask[ij[:,1]] & ((ij[:,1]-ij[:,0])>=sep) & (frame.zscore>=1)].copy()
    selected=eligible.sort_values('zscore',ascending=False,kind='quicksort').head(cap)
    pairs=set(map(tuple,selected[['i','j']].to_numpy(dtype=int)))
    assert len(pairs)==min(cap,len(eligible))
    assert selected.zscore.ge(1).all()
    return pairs,selected,dict(contact_cap=cap,n_eligible=len(eligible),n_selected=len(pairs),
        status='ok',cutoff=float(selected.zscore.min()) if len(selected) else np.nan)


def main():
    assert not OUT.exists(),f'Choose a fresh output directory: {OUT}'
    olddir=ROOT/'for_github_repo/analysis/gremlin_separation'
    oldpro=json.loads((olddir/'provenance.json').read_text())
    assert g.sha256(olddir/'provenance.json')==json.loads((olddir/'COMPLETE.json').read_text())['provenance_sha256']
    for name,digest in oldpro['output_sha256'].items():assert g.sha256(olddir/name)==digest,name
    sources={}
    def track(path,check_old=False):
        key=str(path.relative_to(ROOT));sources[key]=dict(sha256=g.sha256(path),bytes=path.stat().st_size)
        if check_old:assert sources[key]['sha256']==oldpro['source_files'][key]['sha256'],key
        return path
    saved=pd.read_csv(track(ROOT/'for_github_repo/TM_Jaccard_codes/jaccard_results_gremlin.csv',True))
    targets=pd.read_csv(track(ROOT/'for_github_repo/msa-swaps-github/data/targets.csv',True))
    classes=dict(zip(targets.pdb_id.str.lower(),targets.length_class))
    expected={(a,b) for a in classes for b in classes if a!=b and classes[a]==classes[b]}
    assert len(saved)==362 and set(zip(saved.Protein_A,saved.Protein_B))==expected
    OUT.mkdir(parents=True)
    audit=[];results=[];baseline=[]
    with gzip.open(OUT/'selected_contacts.csv.gz','wt',newline='') as f:
        writer=csv.writer(f);writer.writerow(['identifier','condition','separation','i','j','zscore'])
        def load(identifier,condition,score_path,alignment_path):
            mask,nseq=g.gap_mask(track(alignment_path,True))
            score=g.read_scores(track(score_path,True),mask)
            selections={}
            historical=g.select_contacts(score,mask,5)[0]
            for sep in SEPS:
                contacts,selected,info=select(score,mask,sep,150 if classes[identifier.split('_')[0]]==150 else 100)
                selections[sep]=(contacts,info)
                audit.append(dict(identifier=identifier,condition=condition,separation=sep,
                    length_class=classes[identifier.split('_')[0]],alignment_columns=len(mask),
                    retained_columns=int(mask.sum()),n_sequences=nseq,score_rows=len(score),**info))
                for r in selected.itertuples(index=False):
                    writer.writerow([identifier,condition,sep,int(r.i),int(r.j),float(r.zscore)])
            return selections,historical
        native={a:load(a,'native',ROOT/f'gremlin/results_gremlin_proper/full_{a}.csv',
                     ROOT/f'original_proper_msas_af_sampling/{a}/{a}.fas') for a in sorted(classes)}
        for number,row in enumerate(saved.sort_values(['Protein_A','Protein_B']).itertuples(index=False),1):
            a,b=row.Protein_A,row.Protein_B;pair=f'{a}_{b}_aligned';fixed=b=='2mm4' and a in g.FIXED
            alignment=ROOT/(f'2mm4_fixes/{pair}/{pair}.fas' if fixed else f'original_junk_msas/{a}/{pair}/{pair}.fas')
            swapped,oldswap=load(pair,'swapped',ROOT/f'gremlin/results_gremlin_junk/full_{pair}.csv',alignment)
            for ref,protein in (('A',a),('B',b)):
                old=g.compare(oldswap,native[protein][1],ref=='B')
                previous=getattr(row,f'Jaccard_vs_Protein{ref}');oldoff=int(row.Best_Offset_ProteinB) if ref=='B' else 0
                error=abs(old['jaccard']-previous)
                assert error<1e-10 and old['offset']==oldoff
                baseline.append(dict(pair=pair,reference=ref,saved=previous,reproduced=old['jaccard'],absolute_error=error))
                for sep in SEPS:
                    S,si=swapped[sep];T,ti=native[protein][0][sep]
                    result=g.compare(S,T,ref=='B')
                    if si['status']!='ok' or ti['status']!='ok':
                        result.update(jaccard=np.nan,status='insufficient_eligible_for_L')
                    results.append(dict(pair=pair,protein_a=a,protein_b=b,length_class=classes[a],
                        corrected=fixed,separation=sep,reference=ref,cap_swapped=si['contact_cap'],cap_native=ti['contact_cap'],**result))
            if number%60==0 or number==362:print(f'Processed {number}/362 swaps',flush=True)
    frame=pd.DataFrame(results);aud=pd.DataFrame(audit)
    assert len(frame)==2172 and not frame.duplicated(['pair','reference','separation']).any()
    # Exact selection count check for the class at issue, not just configuration.
    a150=aud[aud.length_class==150]
    assert len(a150)==300 and a150.contact_cap.eq(150).all()
    assert (a150.n_selected==np.minimum(a150.n_eligible,150)).all()
    dirs=[]
    for sep in SEPS:
        wide=frame[frame.separation==sep].pivot(index='pair',columns='reference',values='jaccard').dropna()
        dirs.append(dict(separation=sep,n_both_defined=len(wide),n_B_gt_A=int((wide.B>wide.A).sum()),
                         n_equal=int((wide.B==wide.A).sum()),n_A_gt_B=int((wide.A>wide.B).sum())))
    old=pd.read_csv(olddir/'comparisons.csv');old=old[old.separation.isin(SEPS)]
    changes=frame.merge(old,on=['pair','reference','separation'],suffixes=('_revised','_top100'),validate='one_to_one')
    changes['delta_jaccard']=changes.jaccard_revised-changes.jaccard_top100
    changes['offset_changed']=changes.offset_revised!=changes.offset_top100
    unchanged=changes[changes.length_class_revised!=150]
    assert np.allclose(unchanged.jaccard_revised,unchanged.jaccard_top100,atol=1e-14,rtol=0)
    assert (unchanged.offset_revised==unchanged.offset_top100).all()
    for name,table in [('comparisons.csv',frame),('summary.csv',g.summarize(frame)),
                       ('contact_selection_audit.csv',aud),('baseline_validation.csv',pd.DataFrame(baseline)),
                       ('donor_query_direction.csv',pd.DataFrame(dirs)),('changes_from_top100.csv',changes)]:
        table.to_csv(OUT/name,index=False)
    for path in (Path(__file__),Path(g.__file__),olddir/'provenance.json'):track(path)
    pro=dict(created_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        parameters=dict(separations=SEPS,selection='highest-scoring pairs capped at 150 for class 150; 100 otherwise',
            caps_by_length_class={60:100,100:100,150:150},min_zscore=1.0,gap_fraction_strictly_below=.25,
            offsets_A=[0],offsets_B=list(range(-5,6)),insufficient_pairs='retain all eligible pairs when fewer than the cap remain'),
        change_note='Original supplied gremlin_jaccard.py and gremlin_jaccard_v2.py used GREMLIN_TOP_N=100. '
                    'Raise the cap to 150 only for class 150; retain z-score >=1 and all other settings. '
                    'Raw GREMLIN inference was not capped at 100.',
        validation=dict(baseline_comparisons=724,maximum_baseline_error=max(x['absolute_error'] for x in baseline),
            cap_150=150,n_150_input_separation_records=len(a150),
            class150_records_at_cap=int(a150.n_selected.eq(150).sum()),
            unchanged_other_class_comparisons=len(unchanged),
            undefined_comparisons=int(frame.jaccard.isna().sum()),offset_changes=int(changes.offset_changed.sum())),
        source_files=sources,output_sha256={p.name:g.sha256(p) for p in OUT.iterdir() if p.is_file()})
    (OUT/'provenance.json').write_text(json.dumps(pro,indent=2)+'\n')
    (OUT/'COMPLETE.json').write_text(json.dumps(dict(provenance_sha256=g.sha256(OUT/'provenance.json')),indent=2)+'\n')
    print(g.summarize(frame).query("length_class=='all'").to_string(index=False))
    print(json.dumps(pro['validation']))


if __name__=='__main__':main()
