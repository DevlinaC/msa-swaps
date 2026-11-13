#!/usr/bin/env python3
"""Registered contact-map comparisons for all 362 swaps, retaining undefined cases."""
import argparse
import csv
import io
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent/'pilot_lib'))
from common import require,sha256,atomic_bytes,write_json
from run_inference import validated_checkpoint
from analyze_pilot import compare


def write_csv(path,records):
    stream=io.StringIO(); keys=list(dict.fromkeys(k for row in records for k in row))
    writer=csv.DictWriter(stream,fieldnames=keys); writer.writeheader(); writer.writerows(records)
    atomic_bytes(path,stream.getvalue().encode())


def replay(a,b,L,separation,result):
    """Independent explicit set reconstruction of top-L selection and score pairing."""
    pairs=[(i,j) for i in range(L) for j in range(i+separation,L)]
    x=np.array([a[i,j] for i,j in pairs]); y=np.array([b[i,j] for i,j in pairs])
    xcut=sorted(x,reverse=True)[L-1]; ycut=sorted(y,reverse=True)[L-1]
    sx={p for p,v in zip(pairs,x) if v>=xcut}
    sy={p for p,v in zip(pairs,y) if v>=ycut}
    expected=len(sx&sy)/len(sx|sy)
    require(abs(expected-result['ordinary_jaccard'])<1e-12,'Jaccard replay mismatch')
    require(len(sx)==result['top_swap'] and len(sy)==result['top_reference'],'Selection count mismatch')
    # Constant vectors legitimately yield undefined correlations.
    union=sx|sy; idx=[k for k,p in enumerate(pairs) if p in union]
    for key,xx,yy in [('pearson_all_eligible',x,y),('pearson_union_top_L',x[idx],y[idx])]:
        expected=None if np.std(xx)==0 or np.std(yy)==0 else float(np.corrcoef(xx,yy)[0,1])
        require((expected is None and result[key] is None) or
                (expected is not None and result[key] is not None and abs(expected-result[key])<1e-12),
                'Correlation replay mismatch')


def analyze(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    manifest=json.loads((root/'inputs/manifest.json').read_text())
    out=root/'results'; run=json.loads((out/'RUN.json').read_text())
    matrices={}; hashes={}
    for item in manifest['inputs']:
        folder=out/item['id']
        require(validated_checkpoint(folder,run['run_fingerprint'],item['length']), 'Unvalidated inference input')
        matrices[item['id']]=np.load(folder/'contacts.npy',allow_pickle=False)
        hashes[item['id']]=sha256(folder/'contacts.npy')
    records=[]; map_records=[]; count=0
    for item in manifest['inputs']:
        if item['kind']!='swap': continue
        donor_map=item['donor_mapping']; map_records.append(dict(id=item['id'],**donor_map))
        if donor_map['status']=='ok':
            idx=np.arange(donor_map['start'],donor_map['stop']); didx=idx+donor_map['offset']
            a=matrices[item['id']][np.ix_(idx,idx)]
            references={'query':matrices['proper_'+item['query']][np.ix_(idx,idx)],
                        'donor':matrices['proper_'+item['donor']][np.ix_(didx,didx)]}
        for separation in (6,12,24):
            for ref in ('query','donor'):
                base=dict(id=item['id'],query=item['query'],donor=item['donor'],reference=ref,
                          separation=separation,source_depth=item['depth'],selected_depth=item['selected_depth'])
                if donor_map['status']!='ok':
                    records.append(dict(base,status='undefined_mapping',reason=json.dumps(donor_map))); continue
                if sum(max(0,len(idx)-i-separation) for i in range(len(idx)))<len(idx):
                    records.append(dict(base,status='insufficient_eligible_pairs')); continue
                comparison=compare(a,references[ref],len(idx),separation)
                replay(a,references[ref],len(idx),separation,comparison)
                records.append(dict(base,status='ok',**comparison)); count+=1
    require(len(records)==362*3*2,'Incomplete comparison table')
    dest=out/'analysis'; dest.mkdir(exist_ok=True)
    write_csv(dest/'comparisons.csv',records); write_csv(dest/'mapping_audit.csv',map_records)
    metrics=('ordinary_jaccard','pearson_all_eligible','pearson_union_top_L')
    summary=[]; directions=[]
    for sep in (6,12,24):
        for metric in metrics:
            by_ref={}
            for ref in ('query','donor'):
                by_ref[ref]={r['id']:r[metric] for r in records if r['reference']==ref and r['separation']==sep
                             and r['status']=='ok' and r.get(metric) is not None}
                values=list(by_ref[ref].values())
                summary.append(dict(separation=sep,metric=metric,reference=ref,n=len(values),
                              mean=float(np.mean(values)) if values else None,median=float(np.median(values)) if values else None))
            common=by_ref['query'].keys()&by_ref['donor'].keys()
            directions.append(dict(separation=sep,metric=metric,defined=len(common),
                donor_greater=sum(by_ref['donor'][k]>by_ref['query'][k] for k in common),
                query_greater=sum(by_ref['donor'][k]<by_ref['query'][k] for k in common),
                equal=sum(by_ref['donor'][k]==by_ref['query'][k] for k in common)))
    write_csv(dest/'summary.csv',summary); write_csv(dest/'paired_direction.csv',directions)
    fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for ax,metric,title in zip(axes,metrics,['Ordinary top-L Jaccard','Pearson: all eligible pairs','Pearson: union of top-L pairs']):
        groups=[]
        for sep in (6,12,24):
            for ref in ('query','donor'):
                groups.append([r[metric] for r in records if r['separation']==sep and r['reference']==ref
                               and r['status']=='ok' and r.get(metric) is not None])
        boxes=ax.boxplot(groups,positions=[1,2,4,5,7,8],patch_artist=True,showfliers=False,
                         medianprops={'color':'black'})
        for k,box in enumerate(boxes['boxes']): box.set_facecolor('#d48080' if k%2==0 else '#79bfa8')
        ax.set_xticks([1.5,4.5,7.5],['≥6','≥12','≥24']); ax.set_xlabel('Minimum sequence separation')
        ax.set_title(title); ax.set_ylim((0,1) if metric=='ordinary_jaccard' else (-1,1)); ax.grid(axis='y',alpha=.2)
    from matplotlib.patches import Patch
    axes[0].legend(handles=[Patch(facecolor='#d48080',label='Query A'),Patch(facecolor='#79bfa8',label='Donor B')],loc='upper right')
    for suffix in ('png','pdf'):
        b=io.BytesIO(); fig.savefig(b,format=suffix,dpi=200); atomic_bytes(dest/('contact_comparison.'+suffix),b.getvalue())
    plt.close(fig)
    write_json(dest/'VALIDATED.json',dict(status='PASS',matrices=396,comparison_records=len(records),
               independently_replayed=count,undefined=len(records)-count,
               mapping_pass=sum(m['status']=='ok' for m in map_records),input_hashes=hashes,
               run_fingerprint=run['run_fingerprint'],analysis_sha256=sha256(Path(__file__)),
               caution='Repeated proteins are not independent replicates. Top-L union selection can induce negative correlation. No structural accuracy benchmark.',
               selection=manifest['selection'],summary=summary,paired_direction=directions))
    for row in summary: print(row,flush=True)
    print(f'Validated {count}/{len(records)} comparisons; {len(records)-count} undefined',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--root',type=Path,required=True)
    analyze(p.parse_args().root.resolve())
