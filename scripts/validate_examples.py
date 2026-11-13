#!/usr/bin/env python3
"""Check portable examples against saved scientific results and edge cases."""
from pathlib import Path
import numpy as np
import pandas as pd
from compute_contacts import calculate,standard_scores,ROOT
import msatr_contacts as m


def main():
    # Nonempty-union standard metric behavior, including one empty set.
    for a,b,expected in [({(1,7)},{(1,7)},1.),({(1,7)},{(2,8)},0.),
                         ({(1,7),(2,8)},{(1,7),(3,9)},1/3),(set(),{(1,7)},0.)]:
        result=standard_scores(a,b)
        assert abs(result['standard_jaccard']-expected)<1e-12
    try:
        standard_scores(set(),set())
    except ValueError:
        pass
    else:
        raise AssertionError('Empty union must fail explicitly')
    # Show a real methodological difference: swapped ranks within wider sets
    # can produce tolerant overlap=1 while strict top-L contact sets differ.
    n=12; eligible=list(zip(*np.triu_indices(n,k=6)))
    x=np.zeros((n,n));y=x.copy()
    order=list(range(len(eligible)))
    changed=order.copy(); changed[8:16]=order[12:16]+order[8:12]
    for rank,index in enumerate(order):x[eligible[index]]=len(order)-rank
    for rank,index in enumerate(changed):y[eligible[index]]=len(order)-rank
    mask=np.triu(np.ones((n,n),bool),k=6)
    a,b=[set(map(tuple,np.argwhere(m.top(z,mask,n)))) for z in [x,y]]
    assert m.tolerant_jaccard(x,y)==1.
    assert standard_scores(a,b)['standard_jaccard']<1.
    actual=calculate(ROOT/'data/examples')
    for method,suffix in [('GREMLIN','gremlin'),('MSA Transformer','msatransformer')]:
        saved=pd.read_csv(ROOT/f'data/results/jaccard_results_{suffix}.csv')
        saved['pair']=saved.Protein_A+'_'+saved.Protein_B+'_aligned'
        saved=saved.set_index('pair')
        for r in actual[actual.method==method].itertuples():
            assert abs(r.reported_score-saved.loc[r.pair,f'Jaccard_vs_Protein{r.reference}'])<1e-10
            if method=='GREMLIN' and r.reference=='B':
                assert r.offset==saved.loc[r.pair,'Best_Offset_ProteinB']
    for suffix in ['gremlin','msatransformer']:
        table=pd.read_csv(ROOT/f'data/results/jaccard_results_{suffix}.csv')
        assert len(table)==362 and not table.duplicated(['Protein_A','Protein_B']).any()
        scores=table[['Jaccard_vs_ProteinA','Jaccard_vs_ProteinB']].to_numpy()
        assert np.isfinite(scores).all() and ((scores>=0)&(scores<=1)).all()
    table=pd.read_csv(ROOT/'data/results/singleseq_structure_confidence.csv')
    assert len(table)==1700 and set(table.Experiment)=={'singleseq'}
    assert table.groupby(['Model','Target']).size().eq(25).all()
    # Reconcile the retained majority classifications with the summary.
    detail=pd.read_csv(ROOT/'data/results/high_tmscore_reproducible_pairs.csv')
    summary=pd.read_csv(ROOT/'data/results/high_tmscore_summary.csv')
    for row in summary.to_dict('records'):
        block=detail[(detail.Model==row['Model'])&(detail.Reference==row['Reference'])]
        assert len(block)==row['N_Pairs_Total']==362
        assert block.N_Total.sum()==row['N_Predictions_Total']
        assert block['N_Above_0.75'].sum()==row['N_Predictions_Above_0.75']
        assert (block['N_Above_0.75']>block.N_Total/2).sum()==row['N_Pairs_Reproducible']
    print(f'PASS: {len(actual)} example scores, offsets, library agreement, edge cases and table counts')

if __name__=='__main__':
    main()
