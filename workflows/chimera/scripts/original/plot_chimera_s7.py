"""Assemble Figure S7 from audited chimera models and actual PyMOL renders."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/lydia_report_mpl')
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Patch
R=Path(__file__).resolve().parents[2];O=R/'Manuscript/supporting_data/Figure_S7'
rows=json.loads((O/'data/models.json').read_text());summary=json.loads((O/'data/summary.json').read_text())
assert all(x['ctd_strand_residues']>=34 and x['ctd_helix_residues']<=3 for x in rows)
red='#AC0035';gold='#DAA520';grey='#ADADAD'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'font.weight':'bold','axes.labelweight':'bold','axes.titleweight':'bold','pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(8.8,9.4),facecolor='white')
def txt(x,y,t,**kw):fig.text(x,y,t,va='top',**kw)
txt(.04,.98,'A',fontsize=17);txt(.10,.976,'NusG–RfaH chimera',fontsize=14)
ax=fig.add_axes([.11,.86,.79,.08]);ax.set_xlim(0,183);ax.set_ylim(0,1);ax.axis('off')
ax.add_patch(Rectangle((0,.30),130,.36,color=grey));ax.add_patch(Rectangle((130,.30),53,.36,color=red))
ax.text(65,.48,'N-terminal region',ha='center',va='center',fontsize=11)
ax.text(156.5,.48,'RfaH CTD',ha='center',va='center',color='white',fontsize=11)
for x,label,ha in [(0,'1','left'),(130,'131','center'),(183,'183','right')]:ax.text(x,.20,label,ha=ha,va='top',fontsize=10)
txt(.50,.85,'RfaH residues 110–162 → chimera residues 131–183',ha='center',fontsize=10)
txt(.04,.80,'B',fontsize=17)
for i,method in enumerate(['AF2','AF3']):
 x=.09+i*.46;cx=x+.185
 txt(cx,.795,'AF2 (ColabFold)' if method=='AF2' else 'AF3',ha='center',fontsize=14)
 ax=fig.add_axes([x,.425,.37,.35]);im=plt.imread(O/'renders'/f'{method}_full.png');ys,xs=np.where(np.any(im[...,:3]<.97,axis=2));im=im[max(0,ys.min()-25):ys.max()+26,max(0,xs.min()-25):xs.max()+26];ax.imshow(im);ax.axis('off')
 rr=next(r for r in rows if r['method']==method and r['rank']==1)
 txt(cx,.428,'β-sheet CTD',ha='center',color=red,fontsize=12)
 txt(cx,.405,f"Model pLDDT = {rr['mean_plddt']:.1f}  |  CTD pLDDT = {rr['ctd_mean_plddt']:.1f}",ha='center',fontsize=9)
txt(.04,.36,'C',fontsize=17)
a=fig.add_axes([.12,.105,.34,.21]);a.set_title('CTD structural state',fontsize=12,pad=16)
for i,method in enumerate(['AF2','AF3']):
 a.bar(i-.16,0,width=.29,color=gold);a.bar(i+.16,100,width=.29,color=red)
 a.text(i-.16,3,'0/25',ha='center',fontsize=9);a.text(i+.16,103,'25/25',ha='center',fontsize=9)
a.set_ylim(0,118);a.set_yticks([0,25,50,75,100]);a.set_ylabel('Predictions (%)',fontsize=11);a.set_xticks([0,1],['AF2','AF3']);a.spines[['top','right']].set_visible(False)
a.legend(handles=[Patch(facecolor=gold,label='α-helical hairpin'),Patch(facecolor=red,label='β-sheet')],loc='upper center',bbox_to_anchor=(.5,-.18),frameon=False,ncol=1,fontsize=9)
b=fig.add_axes([.61,.105,.32,.21]);b.set_title('CTD confidence',fontsize=12,pad=16)
for i,method in enumerate(['AF2','AF3']):
 vals=np.array([r['ctd_mean_plddt'] for r in rows if r['method']==method]);jit=np.linspace(-.13,.13,len(vals));b.scatter(i+jit,vals,s=19,color=red,alpha=.8,edgecolors='white',linewidths=.3,zorder=3);b.plot([i-.20,i+.20],[np.median(vals)]*2,color='black',lw=2,zorder=4)
b.axhline(70,color='#777777',ls='--',lw=1);b.set_ylim(50,100);b.set_yticks([50,60,70,80,90,100]);b.set_xlim(-.5,1.5);b.set_ylabel('Mean CTD pLDDT',fontsize=11);b.set_xticks([0,1],['AF2','AF3']);b.spines[['top','right']].set_visible(False)
txt(.77,.063,'CTD pLDDT ≥70: 18/25 AF2; 25/25 AF3',ha='center',fontsize=9)
outputs=[]
for ext in ['pdf','png','svg']:
 p=O/f'Figure_S7.{ext}';fig.savefig(p,dpi=350);outputs.append(p)
plt.close(fig)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(O/'provenance/figure.json').write_text(json.dumps(dict(script_sha256=sha(Path(__file__)),inputs={str(p.relative_to(R)):sha(p) for p in [O/'data/models.json',O/'data/summary.json',O/'renders/AF2_full.png',O/'renders/AF3_full.png']},outputs={p.name:sha(p) for p in outputs},palette={'beta_CTD':red,'alpha_hairpin':gold,'N_terminal_region':grey},representative_selection='Independent rank 1 for each prediction method; no selection by fold outcome.'),indent=2)+'\n')
print('Exported Figure S7 PDF, PNG and SVG.')
