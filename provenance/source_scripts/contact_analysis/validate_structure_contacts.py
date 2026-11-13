"""Independent tabular replay and provenance sensitivity, on a compute node."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args();d=a.directory
    provenance=json.loads((d/'provenance.json').read_text())
    for name,digest in provenance['outputs'].items():
        assert hashlib.sha256((d/name).read_bytes()).hexdigest()==digest,name
    frame=pd.read_csv(d/'comparisons.csv'); structures=pd.read_csv(d/'structure_audit.csv')
    inputs=pd.read_csv(d/'prediction_input_audit.csv'); mappings=pd.read_csv(d/'mapping_audit.csv')
    assert len(structures)==792 and set(structures.status)=={'ok'}
    assert structures.missing_contact_atoms.eq(0).all()
    assert structures.alternate_max_distance_difference_A.eq(0).all()
    assert len(mappings)==362
    count=defaultdict(int);hits=defaultdict(int);pairs=defaultdict(list)
    with gzip.open(d/'selected_pairs.csv.gz','rt') as f:
        for r in csv.DictReader(f):
            k=int(r['record_id']); count[k]+=1
            truth=float(r['distance_A'])<=8
            assert truth==bool(int(r['contact']))
            hits[k]+=truth; pairs[k].append((int(r['i']),int(r['j'])))
    max_error=0.
    for r in frame.itertuples(index=False):
        if math.isnan(r.precision):
            assert count[r.record_id]==0
            continue
        assert count[r.record_id]==r.n_selected
        assert hits[r.record_id]==r.intersection
        manual=hits[r.record_id]/count[r.record_id]
        max_error=max(max_error,abs(manual-r.precision))
        assert abs(manual-r.precision)<1e-12
        assert abs(r.intersection/(r.n_selected+r.n_structural-r.intersection)-r.jaccard)<1e-12
    # Matching correspondence cannot change which MSA pairs were selected.
    paired_checks=0
    swap=frame[frame.reference!='native_baseline']
    for key,sub in swap.groupby(['af','method','selection','name','separation']):
        if set(sub.reference)=={'native_A','native_B','swapped_output'}:
            lists=[pairs[r.record_id] for r in sub.itertuples()]
            assert lists[0]==lists[1]==lists[2],str(key)
            paired_checks+=1
    ok={(r.af,r.name):r.status=='compatible_rows' for r in inputs.itertuples()}
    def compatible(r):
        if r['reference']=='native_baseline':
            return ok.get((r['af'],r['query']),False)
        return all(ok.get((r['af'],n),False) for n in (r['query'],r['donor'],r['name']))
    frame['input_compatible']=frame.apply(compatible,axis=1)
    # Every AF3 run passed the input audit; no AF3 comparison may be dropped.
    if inputs[inputs.af=='AF3'].status.eq('compatible_rows').all():
        assert frame.loc[frame.af=='AF3','input_compatible'].all()
    frame.to_csv(d/'comparisons_with_input_flags.csv',index=False)
    summary=[];directions=[]
    for label,subset in [('all',frame),('compatible_inputs',frame[frame.input_compatible])]:
        for keys,sub in subset.groupby(['af','method','selection','reference','separation']):
            row=dict(zip(['af','method','selection','reference','separation'],keys));row.update(subset=label,n=len(sub),n_defined=int(sub.precision.notna().sum()))
            for metric in ('precision','recall','jaccard'):
                row[f'{metric}_mean']=sub[metric].mean();row[f'{metric}_median']=sub[metric].median()
            summary.append(row)
        sw=subset[subset.reference!='native_baseline']
        for keys,sub in sw.groupby(['af','method','selection','separation']):
            for metric in ('precision','jaccard'):
                pivot=sub.pivot(index='name',columns='reference',values=metric)
                for other in ('native_A','swapped_output'):
                    paired=pivot[['native_B',other]].dropna()
                    row=dict(zip(['af','method','selection','separation'],keys))
                    row.update(subset=label,metric=metric,comparison=f'native_B_vs_{other}',n=len(paired),
                        donor_greater=int((paired.native_B>paired[other]).sum()),equal=int((paired.native_B==paired[other]).sum()),
                        donor_lower=int((paired.native_B<paired[other]).sum()),mean_difference=(paired.native_B-paired[other]).mean())
                    directions.append(row)
    pd.DataFrame(summary).to_csv(d/'sensitivity_summary.csv',index=False)
    pd.DataFrame(directions).to_csv(d/'paired_direction.csv',index=False)
    plt.rcParams.update({'font.size':11,'pdf.fonttype':42})
    primary=frame[frame.selection=='topL']
    refs=['native_baseline','native_A','native_B','swapped_output']
    colors=['#777777','#31688e','#b45a2a','#439970']
    for metric,label in [('precision','Contact precision'),('jaccard','Jaccard overlap')]:
        fig,axes=plt.subplots(2,2,figsize=(9,7),sharey=True)
        for row,af in enumerate(('AF2','AF3')):
            for col,method in enumerate(('GREMLIN','MSA Transformer')):
                ax=axes[row,col]
                for r,(ref,color) in enumerate(zip(refs,colors)):
                    values=[primary.loc[(primary.af==af)&(primary.method==method)&(primary.reference==ref)&(primary.separation==sep),metric].dropna().values for sep in (6,12,24)]
                    bp=ax.boxplot(values,positions=np.arange(3)*5+r+1,widths=.7,patch_artist=True,
                        medianprops=dict(color='black'),flierprops=dict(marker='.',markersize=2,alpha=.3))
                    for box in bp['boxes']:box.set_facecolor(color);box.set_alpha(.65)
                ax.set_xticks([2.5,7.5,12.5]);ax.set_xticklabels(['≥6','≥12','≥24'])
                ax.set_title(f'{af} · {method}',fontsize=12);ax.set_ylim(-.03,1.03)
                ax.set_xlabel('Minimum sequence separation')
                if col==0:ax.set_ylabel(label)
                ax.spines[['top','right']].set_visible(False)
        from matplotlib.patches import Patch
        fig.legend([Patch(facecolor=c,alpha=.65) for c in colors],['Native MSA → native fold','Swapped MSA → native A','Swapped MSA → native B','Swapped MSA → swapped output'],loc='lower center',ncol=2,fontsize=10)
        fig.tight_layout(rect=[0,.12,1,1])
        fig.savefig(d/f'structure_contact_{metric}_report.png',dpi=300)
        fig.savefig(d/f'structure_contact_{metric}_report.pdf');plt.close(fig)
    result=dict(job_id=os.getenv('SLURM_JOB_ID'),records=len(frame),defined_precision=int(frame.precision.notna().sum()),
        selected_pairs=sum(count.values()),manual_precision_max_error=max_error,matched_selection_checks=paired_checks,
        mapping_status=mappings.status.value_counts().to_dict(),input_status=inputs.groupby(['af','status']).size().to_string(),
        missing_atoms=int(structures.missing_contact_atoms.sum()),alternative_distance_difference=float(structures.alternate_max_distance_difference_A.max()))
    (d/'validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
    artifacts=['validation.json','sensitivity_summary.csv','paired_direction.csv','comparisons_with_input_flags.csv',
        'structure_contact_precision_report.png','structure_contact_precision_report.pdf',
        'structure_contact_jaccard_report.png','structure_contact_jaccard_report.pdf']
    (d/'VALIDATED.json').write_text(json.dumps(dict(analysis_provenance_sha256=hashlib.sha256((d/'provenance.json').read_bytes()).hexdigest(),
        validation_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        outputs={n:hashlib.sha256((d/n).read_bytes()).hexdigest() for n in artifacts}),indent=2))


if __name__=='__main__':main()
