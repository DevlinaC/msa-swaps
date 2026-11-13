"""Item 1: MSA contact agreement with confidence-ranked AF structures.

All computation is intended for a Slurm compute node. Shared inputs are read
only. Donor registration is inferred from shared MSA rows, never optimized on
contact overlap. All top-L cutoff ties are retained.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import csv
import gzip
import hashlib
import json
import platform
import io
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import Bio
from Bio.PDB import PDBParser, MMCIFParser
from Bio.SeqUtils import seq1
import analyze_gremlin_separation as g

SEPARATIONS = (6,12,24)


def fasta(path):
    records = []; name = None; parts = []
    with (path if hasattr(path,'read') else open(path)) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if line.startswith('>'):
                if name is not None:
                    records.append((name, ''.join(parts)))
                name = line[1:]; parts = []
            else:
                # A3M lowercase insertions and dots do not occupy query columns.
                parts.append(''.join(c for c in line if not c.islower() and c != '.'))
    if name is not None:
        records.append((name, ''.join(parts)))
    return records


def infer_offset(swapped, native):
    """Infer positional correspondence from homolog rows shared by two MSAs."""
    ns = {h:s for h,s in native[1:]}
    shared = [(s,ns[h]) for h,s in swapped[1:] if h in ns]
    if not shared:
        return dict(status='no_shared_homolog_rows', offset=None, shared_rows=0)
    sl, nl = len(swapped[0][1]), len(native[0][1])
    g.require(all(len(s)==sl and len(n)==nl for s,n in shared),'Malformed homolog alignment')
    sm=np.frombuffer(''.join(s for s,n in shared).encode('ascii'),dtype=np.uint8).reshape(len(shared),sl)
    nm=np.frombuffer(''.join(n for s,n in shared).encode('ascii'),dtype=np.uint8).reshape(len(shared),nl)
    scores = []
    for offset in range(-12,13):
        lo, hi = max(0,-offset), min(sl,nl-offset)
        if hi-lo < .8*min(sl,nl):
            continue
        x,y=sm[:,lo:hi],nm[:,lo+offset:hi+offset]
        informative=(x!=ord('-'))|(y!=ord('-'))
        equal=x==y
        matches=int((equal&informative).sum()); total=int(informative.sum()); exact=int(equal.all(axis=1).sum())
        scores.append((matches/total if total else 0, exact, offset,total))
    scores.sort(reverse=True)
    best=scores[0]
    unique=sum(abs(s[0]-best[0])<1e-12 for s in scores)==1
    status='ok' if unique and best[0]>=.99 and len(shared)>=3 else 'ambiguous_or_low_agreement'
    return dict(status=status,offset=best[2],shared_rows=len(shared),identity=best[0],
                exact_rows=best[1], informative_characters=best[3], runner_up_identity=scores[1][0])


def structure(path, sequence):
    parser=MMCIFParser(QUIET=True) if path.suffix=='.cif' else PDBParser(QUIET=True)
    model=next(parser.get_structure('prediction',str(path)).get_models())
    matches=[]
    for chain in model:
        residues=[r for r in chain if r.id[0]==' ']
        observed=''.join(seq1(r.resname,custom_map={'MSE':'M'}) for r in residues)
        if observed==sequence:
            matches.append((chain.id,residues))
    g.require(len(matches)==1, f'Expected one exact query-sequence chain in {path}')
    chain,residues=matches[0]
    xyz=np.full((len(sequence),3),np.nan)
    for i,r in enumerate(residues):
        atom='CA' if sequence[i]=='G' else 'CB'
        if atom in r:
            xyz[i]=r[atom].coord
    distance=np.linalg.norm(xyz[:,None,:]-xyz[None,:,:],axis=2)
    return distance,dict(chain=chain,sequence_length=len(sequence),missing_contact_atoms=int((~np.isfinite(xyz[:,0])).sum()))


def calculate(ij,scores,distances,offset,separation,length,selection):
    """Compare on measurable, registered pairs; never fill missing scores."""
    mapped=ij+offset
    valid=(mapped[:,0]>=0)&(mapped[:,1]<len(distances))&(ij[:,1]-ij[:,0]>=separation)
    ij, scores, mapped=ij[valid],scores[valid],mapped[valid]
    ds=distances[mapped[:,0],mapped[:,1]]
    valid=np.isfinite(ds)
    ij,scores,ds=ij[valid],scores[valid],ds[valid]
    result=dict(L=length,n_eligible=len(scores),offset=offset)
    if selection=='topL':
        if len(scores)<length:
            return dict(result,status='insufficient_eligible_for_L',precision=np.nan,recall=np.nan,jaccard=np.nan),None
        cutoff=np.sort(scores)[-length]; selected=scores>=cutoff
    else:
        selected=np.zeros(len(scores),dtype=bool)
        eligible=np.flatnonzero(scores>=1)
        # Deterministic tie break by i then j for the secondary top-100 definition.
        order=sorted(eligible,key=lambda k:(-scores[k],ij[k,0],ij[k,1]))[:100]
        selected[order]=True; cutoff=float(min(scores[selected])) if selected.any() else np.nan
    truth=ds<=8.0
    tp=int((selected&truth).sum()); predicted=int(selected.sum()); structural=int(truth.sum())
    union=predicted+structural-tp
    result.update(status='ok' if predicted else 'empty_selected', n_selected=predicted,
        n_structural=structural,intersection=tp,union=union,cutoff=float(cutoff),
        precision=tp/predicted if predicted else np.nan,
        recall=tp/structural if structural else np.nan,jaccard=tp/union if union else np.nan)
    return result,(ij, scores,ds,selected,truth)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--shared',type=Path,required=True)
    p.add_argument('--overlay',type=Path)
    p.add_argument('--targets',type=Path,required=True)
    p.add_argument('--expected',type=Path,required=True)
    p.add_argument('--audit',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--limit',type=int)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=False)
    expected=json.loads(a.expected.read_text()); hashes={}
    def source(rel):
        f=a.overlay/rel if a.overlay and (a.overlay/rel).exists() else a.shared/rel
        if rel not in hashes:
            digest=g.sha256(f)
            g.require(digest==expected[rel],f'Source differs from items 2/3: {f}')
            hashes[rel]=dict(path=str(f),sha256=digest)
        return f
    targets=pd.read_csv(a.targets)
    sequences=dict(zip(targets.pdb_id.str.lower(),targets.sequence))
    classes=dict(zip(targets.pdb_id.str.lower(),targets.length_class))
    audit=json.loads((a.audit/'structures.json').read_text())
    models={}; structure_audit=[]; chosen_models={}
    for entry in audit:
        cand=[c for c in entry['candidates'] if not Path(c['path']).name.startswith('.')]
        if not cand:
            structure_audit.append(dict(method=entry['method'],name=entry['name'],status='missing'))
            continue
        # Same deterministic timestamp policy for reruns, recorded explicitly.
        chosen=sorted(cand,key=lambda c:c['run'])[-1]
        f=Path(chosen['path'])
        g.require(g.sha256(f)==chosen['sha256'],f'Structure changed: {f}')
        try:
            dist,details=structure(f,sequences[entry['query']])
        except ValueError as e:
            structure_audit.append(dict(method=entry['method'],name=entry['name'],status='sequence_mismatch',error=str(e)))
            continue
        if entry['method']=='AF2':
            g.require(chosen['recycles']==3 and chosen['rank_by']=='plddt',f'Unexpected AF2 settings: {f}')
            g.require(len(cand)==1, f'Ambiguous AF2 rank 1: {f}')
        else:
            g.require(chosen['n_ranked']==25 and chosen['ranking_ties']==1,f'Incomplete/ambiguous AF3 ranking: {f}')
        alternate_difference=0.0
        for alternate in cand:
            if alternate['path']==chosen['path']:
                continue
            alternate_dist,_=structure(Path(alternate['path']),sequences[entry['query']])
            alternate_difference=max(alternate_difference,float(np.nanmax(abs(alternate_dist-dist))))
        models[(entry['method'],entry['name'])]=dist
        chosen_models[(entry['method'],entry['name'])]=chosen
        record=dict(method=entry['method'],name=entry['name'],status='ok',n_candidates=len(cand),**chosen,**details)
        record['alternate_max_distance_difference_A']=alternate_difference
        record['ranking_sha256']=g.sha256(Path(chosen['ranking']))
        structure_audit.append(record)
    pd.DataFrame(structure_audit).to_csv(a.out/'structure_audit.csv',index=False)
    native_rows={}; native_scores={}
    def read_scores(name,query,donor=None):
        native=donor is None
        fixed=donor=='2mm4' and query in g.FIXED
        relfas=(f'original_proper_msas_af_sampling/{query}/{query}.fas' if native else
                f'2mm4_fixes/{name}/{name}.fas' if fixed else f'original_junk_msas/{query}/{name}/{name}.fas')
        rows=fasta(source(relfas)); g.require(rows[0][1]==sequences[query],f'Query columns do not match target: {name}')
        mask,_=g.gap_mask(source(relfas))
        grel=f"gremlin/results_gremlin_{'proper' if native else 'junk'}/full_{name}.csv"
        gf=g.read_scores(source(grel),mask)
        gij=gf[['i','j']].to_numpy(dtype=int)
        keep=mask[gij[:,0]]&mask[gij[:,1]]
        mrel=(f'original_proper_msas_af_sampling/{query}/full_{query}.csv' if native else
              f'2mm4_fixes/{name}/full_{name}.csv' if fixed else f'original_junk_msas_af_sampling/{query}/{name}/full_{name}.csv')
        matrix=np.loadtxt(source(mrel),delimiter=',')
        g.require(matrix.shape==(len(mask),len(mask)) and np.isfinite(matrix).all(),f'Invalid MSAT shape: {name}')
        i,j=np.triu_indices(len(mask),1)
        return rows,{'GREMLIN':(gij[keep],gf.zscore.to_numpy()[keep]),'MSA Transformer':(np.column_stack((i,j)),matrix[i,j])}
    for q in sequences:
        native_rows[q],native_scores[q]=read_scores(q,q)
    comparisons=[]; mappings=[]; input_audit=[]; pair_count=0
    def audit_prediction_msa(name,rows):
        contact_set=set(s for h,s in rows[1:])
        for af in ('AF2','AF3'):
            c=chosen_models.get((af,name))
            if not c:
                continue
            if af=='AF2':
                f=Path(c['input_msa'])
                predicted=fasta(f) if f.exists() else []
            else:
                inputs=list(Path(c['run']).glob('*_data.json'))
                f=inputs[0] if len(inputs)==1 else None
                data=json.loads(f.read_text()) if f else {}
                proteins=[s['protein'] for s in data.get('sequences',[]) if 'protein' in s]
                predicted=fasta(io.StringIO(proteins[0].get('unpairedMsa',''))) if len(proteins)==1 else []
            pred_set=set(s for h,s in predicted[1:])
            n_shared=len(contact_set&pred_set); denom=min(len(contact_set),len(pred_set))
            fraction=n_shared/denom if denom else np.nan
            input_audit.append(dict(af=af,name=name,path=str(f),sha256=g.sha256(f) if f and f.exists() else None,
                query_match=bool(predicted and predicted[0][1]==rows[0][1]),contact_unique_rows=len(contact_set),
                prediction_unique_rows=len(pred_set),shared_unique_rows=n_shared,shared_fraction_of_smaller=fraction,
                status='compatible_rows' if predicted and predicted[0][1]==rows[0][1] and fraction>=.99 else 'requires_review'))
    with gzip.open(a.out/'selected_pairs.csv.gz','wt',newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['record_id','i','j','structure_i','structure_j','score','distance_A','contact'])
        def run(name,q,donor,score_data,reference,structure_name,offset,shared_length=None,mapping_status='native'):
            # L is the shared physical coordinate span, before gap filtering.
            target_length=len(sequences[donor]) if reference=='native_B' else len(sequences[q])
            length=shared_length or min(len(sequences[q])+offset,target_length)-max(0,offset)
            for af in ('AF2','AF3'):
                if (af,structure_name) not in models:
                    continue
                for method,(ij,scores) in score_data.items():
                    for separation in SEPARATIONS:
                        for selection in (('topL','top100_z1') if method=='GREMLIN' else ('topL',)):
                            result,members=calculate(ij,scores,models[(af,structure_name)],offset,separation,length,selection)
                            record=dict(record_id=len(comparisons),name=name,query=q,donor=donor,af=af,method=method,
                                reference=reference,separation=separation,selection=selection,length_class=classes[q],mapping_status=mapping_status,**result)
                            comparisons.append(record)
                            if members:
                                pairs,values,ds,sel,truth=members
                                for k in np.flatnonzero(sel):
                                    i,j=map(int,pairs[k]); writer.writerow([record['record_id'],i,j,i+offset,j+offset,float(values[k]),float(ds[k]),int(truth[k])])
        for q in sequences:
            audit_prediction_msa(q,native_rows[q])
            run(q,q,None,native_scores[q],'native_baseline',q,0)
        for q in sequences:
            for donor in sequences:
                if q==donor or classes[q]!=classes[donor]:
                    continue
                if a.limit and pair_count>=a.limit:
                    continue
                name=f'{q}_{donor}_aligned'; rows,data=read_scores(name,q,donor)
                audit_prediction_msa(name,rows)
                mapping=infer_offset(rows,native_rows[donor]); mappings.append(dict(name=name,query=q,donor=donor,**mapping))
                span=len(sequences[q])
                if mapping['status']=='ok':
                    off=mapping['offset']; lo=max(0,-off); hi=min(len(sequences[q]),len(sequences[donor])-off)
                    span=hi-lo
                    # Use the same candidate pairs and L for all three swapped comparisons.
                    data={m:(ij[(ij[:,0]>=lo)&(ij[:,1]<hi)],s[(ij[:,0]>=lo)&(ij[:,1]<hi)]) for m,(ij,s) in data.items()}
                    run(name,q,donor,data,'native_B',donor,off,span,mapping['status'])
                run(name,q,donor,data,'native_A',q,0,span,mapping['status'])
                run(name,q,donor,data,'swapped_output',name,0,span,mapping['status'])
                pair_count+=1
                if pair_count%50==0:
                    print(f'Processed {pair_count} swaps',flush=True)
    frame=pd.DataFrame(comparisons)
    frame.to_csv(a.out/'comparisons.csv',index=False)
    pd.DataFrame(mappings).to_csv(a.out/'mapping_audit.csv',index=False)
    pd.DataFrame(input_audit).to_csv(a.out/'prediction_input_audit.csv',index=False)
    summaries=[]
    for key,sub in frame.groupby(['af','method','selection','reference','separation']):
        row=dict(zip(['af','method','selection','reference','separation'],key))
        row.update(n=len(sub),n_defined=int(sub.precision.notna().sum()))
        for metric in ('precision','recall','jaccard'):
            row.update({f'{metric}_{stat}':getattr(sub[metric],stat)() for stat in ('mean','median','min','max')})
        summaries.append(row)
    pd.DataFrame(summaries).to_csv(a.out/'summary.csv',index=False)
    primary=frame[frame.selection=='topL']
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharey=True)
    refs=['native_baseline','native_A','native_B','swapped_output']
    colors=['#777777','#31688e','#b45a2a','#439970']
    for row,af in enumerate(('AF2','AF3')):
        for col,method in enumerate(('GREMLIN','MSA Transformer')):
            ax=axes[row,col]
            for r,(ref,color) in enumerate(zip(refs,colors)):
                values=[primary.loc[(primary.af==af)&(primary.method==method)&(primary.reference==ref)&(primary.separation==sep),'precision'].dropna().values for sep in SEPARATIONS]
                bp=ax.boxplot(values,positions=np.arange(3)*5+r+1,widths=.7,patch_artist=True,
                    medianprops=dict(color='black'),flierprops=dict(marker='.',markersize=2,alpha=.25))
                for box in bp['boxes']: box.set_facecolor(color); box.set_alpha(.65)
            ax.set_xticks([2.5,7.5,12.5]); ax.set_xticklabels(['≥6','≥12','≥24'])
            ax.set_title(f'{af} · {method}'); ax.set_ylim(-.03,1.03)
            ax.set_xlabel('Minimum sequence separation'); ax.set_ylabel('Fraction of top-L pairs in structural contact')
            ax.spines[['top','right']].set_visible(False)
    from matplotlib.patches import Patch
    fig.legend([Patch(facecolor=c,alpha=.65) for c in colors],['Native MSA → native fold','Swapped MSA → native A','Swapped MSA → native B','Swapped MSA → swapped output'],loc='lower center',ncol=2)
    fig.suptitle('MSA-derived contacts compared with rank-1 AlphaFold structures')
    fig.tight_layout(rect=[0,.09,1,.95])
    fig.savefig(a.out/'structure_contact_precision.png',dpi=250)
    fig.savefig(a.out/'structure_contact_precision.pdf'); plt.close(fig)
    provenance=dict(created_utc=datetime.now(timezone.utc).isoformat(),host=platform.node(),job_id=os.getenv('SLURM_JOB_ID'),
        python=platform.python_version(),versions=dict(numpy=np.__version__,pandas=pd.__version__,biopython=Bio.__version__,matplotlib=matplotlib.__version__),
        parameters=dict(distance_A=8,atom='CB; CA for glycine',separations=SEPARATIONS,primary='top L retaining cutoff ties',
            universe='scored pairs passing the input GREMLIN gap mask where applicable, within the common A/B span for swaps, mapped to finite structure distances',
            donor_registration='shared homolog rows; unique offset in -12..12; identity >=0.99; >=3 rows; no contact optimization',
            reruns='lexically latest run directory (timestamp suffix); all candidates in audit',
            selection='AF2 rank_001 by pLDDT, AF3 maximum ranking_score among 25'),
        sources=hashes,scripts={Path(__file__).name:g.sha256(Path(__file__)),Path(g.__file__).name:g.sha256(Path(g.__file__))},
        n_swaps=pair_count,smoke_test=bool(a.limit),n_records=len(frame))
    provenance['outputs']={f.name:g.sha256(f) for f in a.out.iterdir() if f.is_file()}
    (a.out/'provenance.json').write_text(json.dumps(provenance,indent=2))
    (a.out/'COMPLETE.json').write_text(json.dumps(dict(provenance_sha256=g.sha256(a.out/'provenance.json'))))
    print(f'Complete: {len(frame)} comparisons, {pair_count} swaps',flush=True)


if __name__=='__main__':
    main()
