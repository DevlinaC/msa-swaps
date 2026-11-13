"""Plot Figures S2 and S3 from authenticated saved tables, without reanalysis."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/lydia_report_mpl')
from pathlib import Path
import pandas as pd,numpy as np,hashlib,json,argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'Manuscript/supporting_data'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def load(d,key):
 p=R/'for_github_repo/analysis'/d
 pro=json.loads((p/'provenance.json').read_text())
 assert sha(p/'provenance.json')==json.loads((p/'COMPLETE.json').read_text())['provenance_sha256']
 assert sha(p/'comparisons.csv')==pro[key]['comparisons.csv']
 return pd.read_csv(p/'comparisons.csv'),p/'comparisons.csv'
g,gp=load('gremlin_cap150_separation','output_sha256');m,mp=load('msatransformer_separation','outputs');p,pp=load('pearson_topl','output_sha256')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none'})
parser=argparse.ArgumentParser()
parser.add_argument('--figure',type=int,choices=[2,3])
args=parser.parse_args()
for number in ([args.figure] if args.figure else [2,3]):
 if number==2:
  from matplotlib.patches import Patch
  fig,axs=plt.subplots(1,2,figsize=(9,4.7))
  fig.subplots_adjust(left=.09,right=.98,bottom=.28,top=.88,wspace=.30)
  for row,(method,data,metric) in enumerate([('GREMLIN',g,'jaccard'),('MSA Transformer',m,'tolerant_jaccard')]):
   ax=axs[row]
   for ref,color,offset,alpha in [('A','#AC0035',-.18,1),('B','#DAA520',.18,.85)]:
    values=[]
    for sep in [6,12,24]:
     sub=data[(data.reference==ref)&(data.separation==sep)]
     assert len(sub)==362 and sub[metric].notna().all()
     values.append(sub[metric].to_numpy())
    bp=ax.boxplot(values,positions=np.arange(1,4)+offset,widths=.30,patch_artist=True,whis=1.5,
      medianprops=dict(color='black',lw=1),flierprops=dict(marker='o',markersize=2,markerfacecolor=color,markeredgecolor=color,alpha=1 if ref=='A' else .45))
    for box in bp['boxes']:box.set_facecolor(color);box.set_alpha(alpha)
   ax.set_title(f'{chr(65+row)}  {method}',loc='left',fontsize=13,fontweight='bold',pad=12)
   ax.set_ylim(-.025,1.12);ax.set_xlim(.5,3.5);ax.set_yticks([0,.25,.5,.75,1])
   ax.text(2,1.045,'n = 362 per distribution',ha='center',fontsize=9,fontweight='bold')
   ax.set_xticks([1,2,3],['≥6','≥12','≥24'])
   ax.set_ylabel('Standard Jaccard' if row==0 else 'Rank-tolerant Jaccard',fontweight='bold')
   ax.spines[['top','right']].set_visible(False)
   ax.grid(axis='y',color='#e4e4e4',lw=.5);ax.set_axisbelow(True)
   for tick in ax.get_xticklabels()+ax.get_yticklabels():tick.set_fontweight('bold')
  fig.text(.535,.165,'Minimum sequence separation |i − j|',ha='center',fontsize=11,fontweight='bold')
  fig.legend(handles=[Patch(facecolor='#AC0035',edgecolor='black',label='Query'),Patch(facecolor='#DAA520',alpha=.85,edgecolor='black',label='Donor')],loc='lower center',bbox_to_anchor=(.535,.025),ncol=2,frameon=False,prop={'weight':'bold','size':11})
  dest=O/'Figure_S2'
  for ext in ['png','pdf','svg']:fig.savefig(dest/f'Figure_S2_contact_overlap.{ext}',dpi=350)
  plt.close(fig)
  pro=dict(script_sha256=sha(Path(__file__)),source_tables={str(q.relative_to(R)):sha(q) for q in [gp,mp]},outputs={q.name:sha(q) for q in dest.glob('Figure_S2_contact_overlap.*')},query_color='#AC0035',donor_color='#DAA520',query_opacity=1,donor_box_opacity=.85,color_reference='Manuscript/supporting_data/main_figures_embedded/Figure_2_embedded.png',color_reference_sha256=sha(O/'main_figures_embedded/Figure_2_embedded.png'),formatting='Two panels: A GREMLIN, B MSA Transformer; paired query/donor boxes at each separation; shared color legend; no overall title; bold panel text.')
  (dest/'figure_provenance.json').write_text(json.dumps(pro,indent=2)+'\n')
  continue
 from matplotlib.patches import Patch
 fig,axs=plt.subplots(2,1,figsize=(8,8))
 fig.subplots_adjust(left=.14,right=.97,bottom=.18,top=.94,hspace=.40)
 counts={}
 for row,method in enumerate(['GREMLIN','MSA Transformer']):
  ax=axs[row];data=p[p.method==method]
  for ref,color,offset,alpha in [('A','#AC0035',-.20,1),('B','#DAA520',.20,.85)]:
   values=[]
   for sep in [6,12,24]:
    sub=data[(data.reference==ref)&(data.separation==sep)];assert len(sub)==362
    values.append(sub.pearson_r.dropna().to_numpy())
   counts[f'{method}_{ref}']=[len(v) for v in values]
   positions=np.arange(1,4)+offset
   bp=ax.boxplot(values,positions=positions,widths=.32,patch_artist=True,whis=1.5,
     medianprops=dict(color='black',lw=1.2),flierprops=dict(marker='o',markersize=2.5,markerfacecolor=color,markeredgecolor=color,alpha=1 if ref=='A' else .45))
   for box in bp['boxes']:box.set_facecolor(color);box.set_alpha(alpha)
   for x,v in zip(positions,values):ax.text(x,1.10,f'n={len(v)}',ha='center',fontsize=9,fontweight='bold')
  ax.set_title(f'{chr(65+row)}  {method}',loc='left',fontsize=14,fontweight='bold',pad=12)
  ax.set_ylim(-1.08,1.28);ax.set_xlim(.5,3.5);ax.set_yticks([-1,-.5,0,.5,1])
  ax.set_xticks([1,2,3],['≥6','≥12','≥24'])
  ax.set_ylabel('Pearson R',fontsize=12,fontweight='bold')
  ax.spines[['top','right']].set_visible(False)
  ax.grid(axis='y',color='#e4e4e4',lw=.5);ax.set_axisbelow(True)
  for tick in ax.get_xticklabels()+ax.get_yticklabels():tick.set_fontweight('bold');tick.set_fontsize(12)
 fig.text(.555,.105,'Minimum sequence separation |i − j|',ha='center',fontsize=12,fontweight='bold')
 fig.legend(handles=[Patch(facecolor='#AC0035',edgecolor='black',label='Query'),Patch(facecolor='#DAA520',alpha=.85,edgecolor='black',label='Donor')],loc='lower center',bbox_to_anchor=(.555,.025),ncol=2,frameon=False,prop={'weight':'bold','size':12})
 dest=O/'Figure_S3'
 for ext in ['png','pdf','svg']:fig.savefig(dest/f'Figure_S3_topL_Pearson.{ext}',dpi=350)
 plt.close(fig)
 pro=dict(script_sha256=sha(Path(__file__)),source_tables={str(pp.relative_to(R)):sha(pp)},outputs={q.name:sha(q) for q in dest.glob('Figure_S3_topL_Pearson.*')},query_color='#AC0035',donor_color='#DAA520',query_opacity=1,donor_box_opacity=.85,color_reference='Manuscript/supporting_data/main_figures_embedded/Figure_2_embedded.png',color_reference_sha256=sha(O/'main_figures_embedded/Figure_2_embedded.png'),defined_counts=counts,formatting='Two stacked panels: A GREMLIN, B MSA Transformer; paired query/donor boxes; no overall title; bold text. Undefined correlations excluded exactly as before.')
 (dest/'figure_provenance.json').write_text(json.dumps(pro,indent=2)+'\n')
print('Requested figures exported')
