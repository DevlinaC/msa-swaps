#!/usr/bin/env python3
"""Regenerate historical October demo plots from bundled result tables."""
# Historical October 2026 demonstration; see docs/FIGURES.md for manuscript workflows.
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def save(fig,path):
    fig.tight_layout()
    fig.savefig(path,dpi=180,bbox_inches='tight')
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-dir',type=Path,default=ROOT/'data/results')
    ap.add_argument('--out-dir',type=Path,default=ROOT/'outputs/historical_20261030/plots')
    args=ap.parse_args(); args.out_dir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    summary=pd.read_csv(args.data_dir/'high_tmscore_summary.csv')
    fig,axes=plt.subplots(1,2,figsize=(10,4.4),sharey=True)
    settings=[('N_Predictions_Above_0.75','N_Predictions_Total','Individual predictions'),
              ('N_Pairs_Reproducible','N_Pairs_Total','Pairs with ≥13/25 high scores')]
    for ax,(num,den,title) in zip(axes,settings):
        for pos,ref in enumerate(['A','B']):
            for shift,model,color in [(-.18,'AF2','#4C72B0'),(.18,'AF3','#DD8452')]:
                r=summary[(summary.Model==model)&(summary.Reference==ref)].iloc[0]
                pct=100*r[num]/r[den]
                ax.bar(pos+shift,pct,width=.32,color=color,label=model if pos==0 else None)
                ax.text(pos+shift,pct+.7,f'{pct:.2f}%',ha='center',fontsize=9)
        ax.set_xticks([0,1],['A: sequence source','B: MSA source'])
        ax.set_title(title);ax.set_ylim(0,29);ax.legend(frameon=False)
    axes[0].set_ylabel('Percent with TM-score > 0.75')
    save(fig,args.out_dir/'swap_recovery.png')
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),sharey=True)
    for ax,suffix,title in zip(axes,['gremlin','msatransformer'],
            ['GREMLIN: standard Jaccard','MSA Transformer: tolerant contact overlap']):
        table=pd.read_csv(args.data_dir/f'jaccard_results_{suffix}.csv')
        series=[table[f'Jaccard_vs_Protein{r}'].to_numpy() for r in ['A','B']]
        ax.violinplot(series,showmedians=True,showextrema=True)
        ax.set_xticks([1,2],['A: sequence source','B: MSA source'])
        ax.set_title(title,fontsize=11);ax.set_ylim(-.02,1.06)
        for i,values in enumerate(series,1):
            ax.text(i,1.015,f'n={len(values)}',ha='center',fontsize=9)
    axes[0].set_ylabel('Contact overlap (different definitions)')
    save(fig,args.out_dir/'contact_overlap.png')
    table=pd.read_csv(args.data_dir/'singleseq_structure_confidence.csv')
    fig,ax=plt.subplots(figsize=(6,4.5))
    values=[pd.to_numeric(table.loc[table.Model==model,'TM_Score'],errors='raise').to_numpy()
            for model in ['AF2','AF3']]
    ax.violinplot(values,showmedians=True,showextrema=True)
    ax.set_xticks([1,2],[f'{model}\n(n={len(v)})' for model,v in zip(['AF2','AF3'],values)])
    ax.set_ylabel('TM-score against native target');ax.set_ylim(-.02,1.02)
    ax.set_title('Single-sequence predictions without an MSA')
    save(fig,args.out_dir/'singleseq_tmscores.png')
    print(f'Wrote three figures to {args.out_dir}')

if __name__=='__main__':
    main()
