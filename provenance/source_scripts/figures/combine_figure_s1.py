"""Lay out existing PyMOL renderings and unchanged TM-score matrices as S1A/B."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/lydia_s1_mpl')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/lydia_s1_cache')
from pathlib import Path
import csv,json,hashlib,datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

R=Path(__file__).resolve().parents[2];O=R/'Manuscript/supporting_data/Figure_S1'
W=8.6;AH=7.92;BH=3.80;GAP=.12
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.5})
sources={};groups={};seen=[]
for size,count in [(60,14),(100,10),(150,10)]:
    p=O/f'data/tm_max_subgroup_{size}.csv';sources[str(p.relative_to(R))]=sha(p)
    with p.open() as f:
        rows=list(csv.reader(f));ids=rows[0][1:]
    assert len(ids)==count and [x[0] for x in rows[1:]]==ids
    values=np.array([[float(v) if v!='NA' else np.nan for v in row[1:]] for row in rows[1:]])
    assert np.allclose(values,values.T,equal_nan=True)
    assert np.isnan(np.diag(values)).all()
    assert np.nanmax(values)<.6
    groups[size]=(ids,values);seen+=ids
assert len(seen)==len(set(seen))==34
images={}
for size,(ids,_) in groups.items():
    for pid in ids:
        p=R/f'for_supporting_data/outputs/pymol/panels_{size}/{pid}.png'
        sources[str(p.relative_to(R))]=sha(p)
        pixels=plt.imread(p);ys,xs=np.where(np.any(pixels[...,:3]<.97,axis=2));assert len(xs)
        # Same framing-only crop and square padding as the original gallery script.
        crop=pixels[max(0,ys.min()-8):ys.max()+9,max(0,xs.min()-8):xs.max()+9]
        side=max(crop.shape[:2])+24;canvas=np.ones((side,side,pixels.shape[2]),dtype=pixels.dtype)
        y,x=(side-crop.shape[0])//2,(side-crop.shape[1])//2
        canvas[y:y+crop.shape[0],x:x+crop.shape[1]]=crop
        images[pid]=canvas

def text(fig,x,y,s,**kw):
    h=fig.get_figheight();fig.text(x/W,1-y/h,s,va='top',**kw)
def axes(fig,x,y,w,h):
    return fig.add_axes([x/W,1-(y+h)/fig.get_figheight(),w/W,h/fig.get_figheight()])
def draw_a(fig,top=0):
    text(fig,.10,top+.06,'A',fontsize=18,fontweight='bold')
    text(fig,.46,top+.11,'Protein folds',fontsize=14,fontweight='bold')
    for gi,(size,(ids,_)) in enumerate(groups.items()):
        y=top+.50+gi*2.36
        text(fig,W/2,y,f'~{size} residues  |  {len(ids)} proteins',fontsize=12,fontweight='bold',ha='center')
        cols=7 if size==60 else 5
        pitch=1.17;left=(W-cols*pitch)/2
        for i,pid in enumerate(ids):
            row,col=divmod(i,cols);cx=left+(col+.5)*pitch;yy=y+.28+row*1.02
            text(fig,cx,yy,pid,fontsize=11,fontweight='bold',ha='center')
            ax=axes(fig,cx-.385,yy+.19,.77,.77);ax.imshow(images[pid]);ax.set_axis_off()
    # One shared secondary-structure legend for all 34 proteins.
    for x,label,color in [(2.43,'Helix','#E67A0D'),(3.76,'Strand','#0073B3'),(5.26,'Loop','#7A7A7A')]:
        yy=top+7.66
        fig.patches.append(Rectangle((x/W,1-(yy+.11)/fig.get_figheight()),.25/W,.11/fig.get_figheight(),transform=fig.transFigure,color=color,clip_on=False))
        text(fig,x+.34,yy-.01,label,fontsize=11,fontweight='bold')
def draw_b(fig,top=0):
    text(fig,.10,top+.06,'B',fontsize=18,fontweight='bold')
    text(fig,.46,top+.11,'Within-class structural similarity',fontsize=14,fontweight='bold')
    cmap=plt.get_cmap('cividis').copy();cmap.set_bad('#D9D9D9')
    for i,(size,(ids,values)) in enumerate(groups.items()):
        ax=axes(fig,.52+i*2.65,top+.82,2.10,2.10)
        im=ax.imshow(np.ma.masked_invalid(values),vmin=0,vmax=1,cmap=cmap,interpolation='nearest')
        ax.set_xticks(range(len(ids)),ids,rotation=90,fontsize=7.5,fontweight='normal')
        ax.set_yticks(range(len(ids)),ids,fontsize=7.5,fontweight='normal')
        ax.tick_params(length=2,width=.4,pad=3)
        ax.set_title(f'~{size} residues ({len(ids)} proteins)',fontsize=10,pad=12,fontweight='normal')
    cax=axes(fig,8.04,top+.82,.10,2.10)
    cb=fig.colorbar(im,cax=cax,ticks=np.linspace(0,1,6));cb.ax.tick_params(labelsize=8,length=2,width=.4,pad=2)
    cb.set_label('Maximum TM-score',fontsize=9,labelpad=5,fontweight='normal')
    text(fig,W/2,top+3.61,'Grey diagonal: self-comparisons excluded.',fontsize=9,ha='center',color='#444444')

outputs=[]
for stem,h,panels in [('Figure_S1A_protein_folds',AH,'A'),('Figure_S1B_TM_score_heatmaps',BH,'B'),('Figure_S1_combined',AH+GAP+BH,'AB')]:
    fig=plt.figure(figsize=(W,h),facecolor='white')
    if 'A' in panels:draw_a(fig)
    if 'B' in panels:draw_b(fig,AH+GAP if 'A' in panels else 0)
    for ext in ['pdf','png','svg']:
        p=O/f'{stem}.{ext}';fig.savefig(p,dpi=350,facecolor='white');outputs.append(p)
    plt.close(fig)
pro=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),script_sha256=sha(Path(__file__)),
    panel_A='34 original cartoon renderings, grouped by nominal length; one label and one shared secondary-structure legend',
    panel_B='Three unchanged within-class maximum-normalized TM-score matrices; one label and one shared colorbar',
    source_sha256=sources,protein_ids=seen,counts={'60':14,'100':10,'150':10,'total':34,'unordered_pairs':181},
    changes='Layout and labels only. Same image-framing crop as original montage; no coordinate, orientation, color, TM-score, normalization or threshold changes.',
    dimensions_inches={'A':[W,AH],'B':[W,BH],'combined':[W,AH+GAP+BH]},
    typography_points={'panel_letters':18,'panel_titles':14,'A_group_headings':12,'A_protein_labels':11,'A_legend':11,'B_group_headings':10,'B_protein_labels':7.5,'B_colorbar_ticks':8,'B_colorbar_label':9},
    panel_B_typography_note='Regular-weight protein labels and single-line class headings; wider whitespace between heatmaps.',
    removed_annotation='Structures are independently oriented and zoomed.',
    output_sha256={p.name:sha(p) for p in outputs})
(O/'provenance/combined_figure_provenance.json').write_text(json.dumps(pro,indent=2)+'\n')
print('Exported',*[p.name for p in outputs],sep='\n')
