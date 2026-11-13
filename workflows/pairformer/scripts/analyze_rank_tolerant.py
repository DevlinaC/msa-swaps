#!/usr/bin/env python3
"""Rank-tolerance sensitivity analysis of saved, registered Pairformer maps."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/msa_pairformer_plot_mpl')
import argparse, ast, csv, hashlib, io, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0,str(Path(__file__).resolve().parent/'pilot_lib'))
from common import require, sha256, atomic_bytes, write_json


def write_csv(path, rows):
    buf=io.StringIO(); writer=csv.DictWriter(buf,fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows); atomic_bytes(path,buf.getvalue().encode())


def main(root, legacy):
    dest=root/'results/rank_tolerant'; dest.mkdir(exist_ok=True)
    sources=[root/'inputs/manifest.json',root/'results/analysis/comparisons.csv',root/'results/analysis/VALIDATED.json',legacy,root/'figure_palette.json']
    hashes={str(p):sha256(p) for p in sources}
    manifest=json.loads(sources[0].read_text()); validation=json.loads(sources[2].read_text())
    # Execute the unmodified compare function from the historical analysis,
    # without executing that script's main or importing its plotting dependencies.
    tree=ast.parse(legacy.read_text()); node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='compare')
    scope={'np':np}; exec(compile(ast.Module(body=[node],type_ignores=[]),str(legacy),'exec'),scope)
    historical_compare=scope['compare']
    matrices={}
    for item in manifest['inputs']:
        p=root/'results'/item['id']/'contacts.npy'
        require(sha256(p)==validation['input_hashes'][item['id']],f'Matrix hash mismatch: {p}')
        matrices[item['id']]=np.load(p,allow_pickle=False)
    with sources[1].open() as f: prior={(r['id'],int(r['separation']),r['reference']):r for r in csv.DictReader(f)}
    records=[]
    for item in manifest['inputs']:
        if item['kind']!='swap': continue
        mapping=item['donor_mapping']; require(mapping['status']=='ok','Mapping not validated')
        idx=np.arange(mapping['start'],mapping['stop']); didx=idx+mapping['offset']; L=len(idx)
        a=matrices[item['id']][np.ix_(idx,idx)]
        refs={'query':matrices['proper_'+item['query']][np.ix_(idx,idx)],'donor':matrices['proper_'+item['donor']][np.ix_(didx,didx)]}
        for sep in (6,12,24):
            ii,jj=np.triu_indices(L,k=sep); x=a[ii,jj]; wide=int(1.5*L)
            require(len(x)>=wide,'Insufficient eligible pairs')
            for ref,b in refs.items():
                y=b[ii,jj]
                masks=[v>=np.partition(v,-k)[-k] for v,k in ((x,L),(y,L),(x,wide),(y,wide))]
                A,B,AW,BW=masks; union=A|B; credited=(A&BW)|(B&AW)
                ordinary=float(np.sum(A&B)/sum(union)); tolerant=float(sum(credited)/sum(union))
                old=prior[item['id'],sep,ref]
                require(abs(ordinary-float(old['ordinary_jaccard']))<1e-12,'Ordinary score changed')
                # Historical implementation uses full sorting and explicit pair-set replay.
                historical=historical_compare(a,b,sep)
                require(abs(tolerant-historical['tolerant_jaccard'])<1e-12,'Historical metric mismatch')
                require(0<=ordinary<=tolerant<=1,'Invalid metric ordering')
                records.append(dict(id=item['id'],query=item['query'],donor=item['donor'],reference=ref,separation=sep,L=L,wide_k=wide,source_depth=item['depth'],selected_depth=item['selected_depth'],ordinary_jaccard=ordinary,tolerant_jaccard=tolerant,delta=tolerant-ordinary,top_swap=int(A.sum()),top_reference=int(B.sum()),wide_swap=int(AW.sum()),wide_reference=int(BW.sum()),credited=int(credited.sum()),union=int(union.sum())))
    require(len(records)==2172,'Incomplete comparison table')
    write_csv(dest/'comparisons.csv',records)
    summary=[]; direction=[]
    for sep in (6,12,24):
        for ref in ('query','donor'):
            group=[r for r in records if r['separation']==sep and r['reference']==ref]
            for metric in ('ordinary_jaccard','tolerant_jaccard'):
                vals=[r[metric] for r in group]
                summary.append(dict(separation=sep,reference=ref,metric=metric,n=len(vals),mean=float(np.mean(vals)),median=float(np.median(vals)),minimum=min(vals),maximum=max(vals)))
        group={ref:{r['id']:r['tolerant_jaccard'] for r in records if r['separation']==sep and r['reference']==ref} for ref in ('query','donor')}
        direction.append(dict(separation=sep,donor_greater=sum(group['donor'][k]>v for k,v in group['query'].items()),equal=sum(group['donor'][k]==v for k,v in group['query'].items()),query_greater=sum(group['donor'][k]<v for k,v in group['query'].items())))
    write_csv(dest/'summary.csv',summary); write_csv(dest/'paired_direction.csv',direction)
    palette=json.loads((root/'figure_palette.json').read_text())['palette']
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(9,4.5),constrained_layout=True)
    for ax,metric,title in zip(axes,('ordinary_jaccard','tolerant_jaccard'),('Ordinary top-L Jaccard','Rank-tolerant top-L Jaccard')):
        for si,sep in enumerate((6,12,24)):
            for ri,ref in enumerate(('query','donor')):
                group=[r for r in records if r['separation']==sep and r['reference']==ref]
                vals=[r[metric] for r in group]; pos=1+3*si+ri
                xx=[pos+.56*(int(hashlib.sha256(r['id'].encode()).hexdigest()[:8],16)/0xffffffff-.5) for r in group]
                ax.scatter(xx,vals,s=10,alpha=.42,color=palette[ref],edgecolors='none')
                ax.plot([pos-.31,pos+.31],[np.median(vals)]*2,color='black',lw=1.7)
        ax.set(xticks=[1.5,4.5,7.5],xticklabels=['≥6','≥12','≥24'],xlabel='Minimum sequence separation',title=title,ylim=(-.025,1.025),xlim=(.4,8.6))
        ax.grid(axis='y',alpha=.2); ax.set_axisbelow(True)
    axes[0].legend(handles=[Patch(facecolor=palette[r],label=l) for r,l in [('query','Query A'),('donor','Donor B')]],loc='upper center',bbox_to_anchor=(.5,1.23),ncol=2,borderaxespad=0)
    for suffix in ('png','pdf','svg'):
        b=io.BytesIO(); fig.savefig(b,format=suffix,dpi=300); atomic_bytes(dest/('jaccard_comparison.'+suffix),b.getvalue())
    plt.close(fig)
    require(all(sha256(Path(p))==h for p,h in hashes.items()),'Source file changed')
    write_json(dest/'VALIDATED.json',dict(status='PASS',comparison_records=len(records),matrices_checked=396,historical_function_replays=len(records),ordinary_score_replays=len(records),source_hashes=hashes,code_sha256=sha256(Path(__file__)),mapping='Existing homolog registration retained; no center recropping.',definition='|(A intersect Bwide) union (B intersect Awide)| / |A union B|; A/B top L, wide top floor(1.5L), cutoff ties retained.',figure='All 362 points per group; horizontal jitter only; black median bars; manuscript colors.',summary=summary,paired_direction=direction,output_hashes={p.name:sha256(p) for p in dest.iterdir() if p.is_file() and p.name!='VALIDATED.json'}))
    for row in summary: print(row)
    print('PAIRED',direction)
    for sep in (6,12,24):
        group=[r for r in records if r['separation']==sep and r['reference']=='donor' and r['donor']=='2mm4']
        print('2mm4',sep,'n',len(group),'ordinary',np.mean([r['ordinary_jaccard'] for r in group]),'tolerant',np.mean([r['tolerant_jaccard'] for r in group]))

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent)
    p.add_argument('--legacy',type=Path,default=Path(__file__).resolve().parent/'reference/analyze_msatransformer_separation.py')
    args=p.parse_args(); main(args.root.resolve(),args.legacy.resolve())
