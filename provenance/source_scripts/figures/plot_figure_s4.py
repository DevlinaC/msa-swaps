"""Format primary precision and Jaccard as panels A and B; no reanalysis."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/lydia_report_mpl')
from pathlib import Path
import json,hashlib
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
R=Path(__file__).resolve().parents[2];O=R/'Manuscript/supporting_data/Figure_S4'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p=O/'data/comparisons.csv';source_hash=sha(p)
df=pd.read_csv(p);df=df[df.selection=='topL']
refs=['native_baseline','native_A','native_B','swapped_output']
colors=['#BFBFBF','#AC0035','#DAA520','#737373'];alphas=[1,1,.85,1]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
counts=[];outputs=[]
for stem,indices in [('Figure_S4_combined',[0,1]),('Figure_S4A',[0]),('Figure_S4B',[1])]:
 combined=len(indices)==2
 fig,axs=plt.subplots(4 if combined else 2,2,figsize=(9,13.6 if combined else 7.2))
 fig.subplots_adjust(left=.10,right=.98,bottom=.17 if combined else .25,top=.94 if combined else .88,hspace=.60 if combined else .45,wspace=.27)
 for slot,mi in enumerate(indices):
  metric,label=[('precision','Contact precision'),('jaccard','Jaccard overlap')][mi]
  fig.text(.10,(.982 if mi==0 else .557) if combined else .97,f'{chr(65+mi)}  {label}',fontsize=15,fontweight='bold',va='top')
  for ai,af in enumerate(['AF2','AF3']):
   for col,method in enumerate(['GREMLIN','MSA Transformer']):
    ax=axs[slot*2+ai,col];centers=np.arange(1,4)
    for ri,(ref,color,alpha) in enumerate(zip(refs,colors,alphas)):
     values=[]
     for sep in [6,12,24]:
      sub=df[(df.af==af)&(df.method==method)&(df.reference==ref)&(df.separation==sep)]
      assert len(sub)==(34 if ref=='native_baseline' else 362) and sub[metric].notna().all()
      values.append(sub[metric].to_numpy())
      if combined:counts.append(dict(metric=metric,af=af,method=method,reference=ref,separation=sep,n=len(sub)))
     bp=ax.boxplot(values,positions=centers+(ri-1.5)*.20,widths=.17,patch_artist=True,whis=1.5,medianprops=dict(color='black',lw=1),flierprops=dict(marker='o',markersize=1.8,markerfacecolor=color,markeredgecolor=color,alpha=.35))
     for box in bp['boxes']:box.set_facecolor(color);box.set_alpha(alpha)
    ax.set_title(f'{af} · {method}',fontsize=11,fontweight='bold',pad=8)
    ax.set_ylim(-.03,1.03);ax.set_xlim(.45,3.55);ax.set_yticks([0,.25,.5,.75,1])
    ax.set_xticks([1,2,3],['≥6','≥12','≥24'])
    if col==0:ax.set_ylabel(label,fontsize=11,fontweight='bold')
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',color='#e4e4e4',lw=.5);ax.set_axisbelow(True)
    for tick in ax.get_xticklabels()+ax.get_yticklabels():tick.set_fontweight('bold');tick.set_fontsize(10)
 fig.text(.54,.12 if combined else .16,'Minimum sequence separation |i − j|',ha='center',fontsize=12,fontweight='bold')
 labels=['Native MSA → native fold','Swapped MSA → native query (A)','Swapped MSA → native donor (B)','Swapped MSA → swapped output']
 fig.legend(handles=[Patch(facecolor=c,alpha=a,edgecolor='black',label=l) for c,a,l in zip(colors,alphas,labels)],loc='lower center',bbox_to_anchor=(.54,.025),ncol=2,frameon=False,prop={'weight':'bold','size':10},columnspacing=1.5)
 for ext in ['pdf','png','svg']:
  path=O/f'{stem}.{ext}';fig.savefig(path,dpi=350);outputs.append(path)
 plt.close(fig);assert sha(p)==source_hash
pro=dict(script_sha256=sha(Path(__file__)),source_tables={str(p.relative_to(R)):source_hash},selection='topL, unchanged primary analysis',panels={'A':'Contact precision','B':'Jaccard overlap'},palette=dict(zip(refs,colors)),box_alpha=dict(zip(refs,alphas)),counts=counts,outputs={p.name:sha(p) for p in outputs})
(O/'figure_provenance.json').write_text(json.dumps(pro,indent=2)+'\n')
print('Exported Figure S4 precision and Jaccard in PDF, PNG and SVG.')
