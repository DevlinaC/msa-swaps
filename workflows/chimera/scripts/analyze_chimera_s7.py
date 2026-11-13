"""Audit all 50 chimera models; map the RfaH CTD by exact sequence identity."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,csv,hashlib,shutil
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import seq1
from tmtools import tm_align
import pymol2
from _paths import FIGURE_DIR,project_root
R=project_root();RUN=R/'new_af3_runs/chimera_20261106';O=FIGURE_DIR
for sub in ['data','structures','provenance']:(O/sub).mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
h=str(next(SeqIO.parse(RUN/'inputs/Rfah_NusG_hybrid.fasta','fasta')).seq)
r=str(next(SeqIO.parse(R/'fasta_files/RfaH.fasta','fasta')).seq)
assert len(h)==183 and len(r)==162 and h[130:]==r[109:]
n=0
for a,b in zip(h[::-1],r[::-1]):
 if a!=b:break
 n+=1
assert n==60
cf=sorted((RUN/'results/colabfold').glob('*unrelaxed_rank_*.pdb'));assert len(cf)==25
base=RUN/'results/af3/Rfah_NusG_hybrid';ranking=list(csv.DictReader((base/'Rfah_NusG_hybrid_ranking_scores.csv').open()))
ranking.sort(key=lambda x:float(x['ranking_score']),reverse=True);assert len(ranking)==25
specs=[('AF2',i+1,p,None) for i,p in enumerate(cf)]
for i,row in enumerate(ranking):
 d=base/f"seed-{row['seed']}_sample-{row['sample']}";specs.append(('AF3',i+1,next(d.glob('*model.cif')),float(row['ranking_score'])))
rows=[];reference={};sources={};representatives={}
with pymol2.PyMOL() as pm:
 c=pm.cmd;c.set('max_threads',1)
 for method,rank,p,score in specs:
  c.load(str(p),'m');c.dss('m');atoms=c.get_model('m and name CA').atom
  assert len(atoms)==183 and ''.join(seq1(a.resn) for a in atoms)==h,p
  assert [int(a.resi) for a in atoms]==list(range(1,184))
  ca=np.array([a.coord for a in atoms]);ctd=ca[130:];ss=''.join(a.ss or 'L' for a in atoms[130:]);conf=np.array([a.b for a in atoms])
  if rank==1:
   reference[method]=ctd
   dest=O/'structures'/f'{method}_rank1{p.suffix}';shutil.copyfile(p,dest);representatives[method]=os.path.relpath(dest,R)
  result=tm_align(ctd,reference[method],h[130:],h[130:])
  rows.append(dict(method=method,rank=rank,source=str(p.relative_to(R)),sha256=sha(p),ranking_score=score,mean_plddt=float(conf.mean()),ctd_mean_plddt=float(conf[130:].mean()),ctd_min_plddt=float(conf[130:].min()),ctd_helix_residues=ss.count('H'),ctd_strand_residues=ss.count('S'),ctd_loop_residues=ss.count('L'),ctd_ss=ss,ctd_tm_to_rank1=float(result.tm_norm_chain1)))
  sources[str(p.relative_to(R))]=sha(p);c.delete('all')
summary={}
for method in ['AF2','AF3']:
 rr=[x for x in rows if x['method']==method]
 summary[method]={key:dict(min=min(x[key] for x in rr),median=float(np.median([x[key] for x in rr])),max=max(x[key] for x in rr)) for key in ['mean_plddt','ctd_mean_plddt','ctd_helix_residues','ctd_strand_residues','ctd_tm_to_rank1']}
 summary[method]['n_models']=len(rr);summary[method]['ctd_plddt_ge70']=sum(x['ctd_mean_plddt']>=70 for x in rr)
 summary[method]['strand_majority']=sum(x['ctd_strand_residues']>26.5 for x in rr)
for p in [RUN/'inputs/Rfah_NusG_hybrid.fasta',R/'fasta_files/RfaH.fasta',base/'Rfah_NusG_hybrid_ranking_scores.csv',RUN/'results/colabfold/config.json']:
 sources[str(p.relative_to(R))]=sha(p)
mapping=dict(chimera_length=183,reference_rfah_length=162,ctd_chimera=[131,183],ctd_rfah=[110,162],ctd_identity='53/53',exact_matching_suffix_chimera=[124,183],exact_matching_suffix_rfah=[103,162],note='CTD interval follows manuscript residues 110–162 in wild-type RfaH, mapped exactly to chimera 131–183. Longest exact C-terminal match is 60 residues. Exact engineered splice boundary is not inferred from overlapping sequence identity.')
(O/'data/models.json').write_text(json.dumps(rows,indent=2)+'\n');(O/'data/summary.json').write_text(json.dumps(summary,indent=2)+'\n');(O/'data/domain_mapping.json').write_text(json.dumps(mapping,indent=2)+'\n')
(O/'provenance/analysis.json').write_text(json.dumps(dict(script_sha256=sha(Path(__file__)),sources=sources,representatives=representatives,secondary_structure='PyMOL DSS; CA residue assignments H helix, S strand, L loop. Geometric annotation, not a learned fold classifier.',confidence='Unweighted mean of CA pLDDT over chimera residues 131–183; all models retained.',structural_comparison='tmtools TM-align, CTD normalized by 53 residues, compared with independently rank-selected model within each prediction method.'),indent=2)+'\n')
print(json.dumps(summary,indent=2))
