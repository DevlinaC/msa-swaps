#!/usr/bin/env python3
"""Build model-labelled paired score records from saved AF2/AF3 audits; no new alignments."""
from pathlib import Path
import json
import hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'data/publication/structure_scores'
def main():
    records=[]
    af2=pd.read_csv(D/'af2_tmalign_records.csv')
    af2=af2[af2.condition=='swapped']
    if len(af2)!=18100 or not af2.matches_archived_prediction.all():raise ValueError('AF2 population/link mismatch')
    for row in af2.itertuples():
        records.append(dict(af='AF2',query=row.query,donor=row.donor,
            prediction_id=row.prediction,model=str(row.model),seed=row.seed,sample='',
            reference=row.reference_role,tm_score=row.tm_chain2,other_normalization_score=row.tm_chain1,
            normalization='Chain_2; matches historical extracted score',source_record=row.source_member))
    af3=pd.read_csv(D/'af3_saved_scores_with_links.csv')
    af3=af3[af3.Job_Type=='junk']
    if len(af3)!=18100:raise ValueError('AF3 population mismatch')
    for row in af3.itertuples():
        query,donor,*_=row.Job_Label.split('_')
        if row.Compared_Protein not in (query,donor):raise ValueError('Unexpected AF3 reference')
        records.append(dict(af='AF3',query=query,donor=donor,
            prediction_id=f'{row.Job_Label}/seed-{row.Seed_Number}_sample-{row.Sample_Number}',
            model='',seed=row.Seed_Number,sample=row.Sample_Number,
            reference='query' if row.Compared_Protein==query else 'donor',tm_score=row.TM_Score,
            other_normalization_score=None,normalization='Saved score; supplied parser intends reference/Chain_1; historical execution unconfirmed',
            source_record=row.candidate_model_paths))
    long=pd.DataFrame(records)
    key=['af','query','donor','prediction_id']
    if long.duplicated(key+['reference']).any():raise ValueError('Duplicate model/reference')
    pivot=long.pivot(index=key,columns='reference',values='tm_score').rename(columns={'query':'tm_query','donor':'tm_donor'}).reset_index()
    if len(pivot)!=18100 or pivot[['tm_query','tm_donor']].isna().any().any():raise ValueError('Incomplete paired scores')
    metadata=long[long.reference=='query'].drop(columns=['reference','tm_score','other_normalization_score'])
    wide=pivot.merge(metadata,on=key,validate='one_to_one')
    if wide.groupby('af').size().to_dict()!={'AF2':9050,'AF3':9050}:raise ValueError('AF method counts')
    if not (wide.groupby(['af','query','donor']).size()==25).all():raise ValueError('Per-pair sample counts')
    if not wide[['tm_query','tm_donor']].ge(0).all().all() or not wide[['tm_query','tm_donor']].le(1).all().all():raise ValueError('Score range')
    wide['query_ge_0_6']=wide.tm_query>=.6;wide['donor_ge_0_6']=wide.tm_donor>=.6
    out=D/'derived';out.mkdir(parents=True,exist_ok=True)
    wide.to_csv(out/'swapped_predictions.csv',index=False)
    long.to_csv(out/'swapped_comparisons.csv',index=False)
    summary=[]
    for af,group in wide.groupby('af'):
        summary.append(dict(af=af,predictions=len(group),directed_pairs=group[['query','donor']].drop_duplicates().shape[0],threshold=.6,operator='>=',query_ge_0_6=int(group.query_ge_0_6.sum()),donor_ge_0_6=int(group.donor_ge_0_6.sum()),maximum_donor_tm=float(group.tm_donor.max())))
    pd.DataFrame(summary).to_csv(out/'summary.csv',index=False)
    sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in D.glob('*.csv')}
    (out/'provenance.json').write_text(json.dumps(dict(sources=sources,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limitation='Reconciliation of saved scores only; no new TM-align. AF3 historical execution normalization and reference provenance remain to be confirmed.'),indent=2)+'\n')
    print(pd.DataFrame(summary).to_string(index=False))
if __name__=='__main__':main()
