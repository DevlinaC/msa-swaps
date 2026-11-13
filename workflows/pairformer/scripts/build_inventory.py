#!/usr/bin/env python3
"""Create shareable path indexes and validate saved matrices without inference."""
import argparse, csv, hashlib, io, json, os
from datetime import datetime, timezone
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def tsv(path, rows):
    stream=io.StringIO(); w=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t')
    w.writeheader(); w.writerows(rows); path.write_text(stream.getvalue())


def category(p):
    s=str(p)
    if s.startswith('results/swap_'): return 'swapped prediction or checkpoint'
    if s.startswith('results/proper_'): return 'proper reference prediction or checkpoint'
    if '/archive_' in s: return 'superseded figure archive'
    if s.startswith('results/rank_tolerant/'): return 'rank-tolerant sensitivity analysis'
    if s.startswith('results/configuration_comparison/'): return 'paired 256-biased versus 512-unweighted comparison'
    if s.startswith('results/analysis/'): return 'ordinary Jaccard and Pearson analysis'
    if s.startswith('inputs/'): return 'selected MSA or input metadata'
    if s.startswith('logs/'): return 'Slurm execution log'
    if p.suffix=='.md': return 'documentation'
    if p.suffix in ('.py','.sbatch'): return 'code or batch script'
    return 'configuration or provenance'


def main(root):
    out=root/'inventory'; out.mkdir(exist_ok=True)
    manifest=json.loads((root/'inputs/manifest.json').read_text())
    valid=json.loads((root/'results/analysis/VALIDATED.json').read_text())
    run=json.loads((root/'results/RUN.json').read_text())
    assert manifest['selection']['cap']==run['max_depth']==512 and run['query_biasing'] is False
    rows=[]
    for item in sorted(manifest['inputs'],key=lambda i:i['id']):
        ident=item['id']; folder=Path('results')/ident; npy=folder/'contacts.npy'; csv_path=folder/'full_contacts.csv'
        a=np.load(root/npy,allow_pickle=False); b=np.loadtxt(root/csv_path,delimiter=',')
        assert a.shape==(item['length'],)*2 and np.isfinite(a).all() and np.min(a)>=0 and np.max(a)<=1,ident
        assert np.max(abs(a-a.T))<=1e-6 and np.array_equal(a,b),ident
        h=sha(root/npy); complete=json.loads((root/folder/'COMPLETE.json').read_text())
        assert complete['query_biasing'] is False and complete['max_depth']==512 and complete['sequence_weighting']=='uniform',ident
        assert h==valid['input_hashes'][ident]==complete['contacts_sha256'],ident
        msa=Path(item['frozen_path']); metadata=Path('inputs')/(ident+'.json')
        assert sha(root/msa)==item['selected_sha256'],ident
        assert json.loads((root/metadata).read_text())==item,ident
        mapping=item.get('donor_mapping') or {}
        rows.append(dict(id=ident,kind=item['kind'],query=item['query'],donor=item['donor'],matrix_rows=a.shape[0],matrix_columns=a.shape[1],dtype=str(a.dtype),source_depth=item['depth'],selected_depth=item['selected_depth'],corrected_input=item.get('repaired',False),npy_path=str(npy),csv_path=str(csv_path),checkpoint_path=str(folder/'COMPLETE.json'),selected_msa_path=str(msa),metadata_path=str(metadata),proper_query_npy=f"results/proper_{item['query']}/contacts.npy",proper_donor_npy=f"results/proper_{item['donor']}/contacts.npy",mapping_start=mapping.get('start',''),mapping_stop_exclusive=mapping.get('stop',''),donor_offset=mapping.get('offset',''),npy_bytes=(root/npy).stat().st_size,csv_bytes=(root/csv_path).stat().st_size,npy_sha256=h,csv_sha256=sha(root/csv_path),selected_msa_sha256=item['selected_sha256']))
    for row in rows:
        row.update(max_depth=512,query_biasing=False,sequence_weighting='uniform')
    swapped=[r for r in rows if r['kind']=='swap']; assert len(rows)==396 and len(swapped)==362
    tsv(out/'all_matrices.tsv',rows); tsv(out/'swapped_matrices.tsv',swapped)
    for fmt in ('csv','npy'):
        (out/f'swapped_{fmt}_files.txt').write_text(''.join(r[f'{fmt}_path']+'\n' for r in swapped))
    files=[]; directories=[]
    for p in sorted(root.rglob('*')):
        rel=p.relative_to(root)
        if any(x.startswith('.') or x in ('__pycache__','inventory') for x in rel.parts): continue
        if p.is_dir(): directories.append(str(rel)); continue
        if not p.is_file(): continue
        files.append(dict(path=str(rel),bytes=p.stat().st_size,sha256=sha(p),purpose=category(rel)))
    tsv(out/'files.tsv',files)
    (out/'filesystem.txt').write_text('Paths relative to msa_pairformer_512_unweighted/; caches, hidden files, and inventory/ excluded.\n\nDIRECTORIES\n'+''.join(p+'/\n' for p in directories)+'\nFILES\n'+''.join(f['path']+'\n' for f in files))
    (out/'SHA256SUMS').write_text(''.join(f"{f['sha256']}  {f['path']}\n" for f in files))
    report=dict(status='PASS',created_utc=datetime.now(timezone.utc).isoformat(),root_convention='All paths relative to project root',matrices=396,swapped=362,proper=34,csv_npy_exact_matches=396,original_matrix_hash_matches=396,selected_msa_hash_matches=396,indexed_files=len(files),indexed_directories=len(directories),indexed_bytes=sum(f['bytes'] for f in files),swapped_csv_bytes=sum(r['csv_bytes'] for r in swapped),swapped_npy_bytes=sum(r['npy_bytes'] for r in swapped),excluded=['inventory/ (generated indexes listed separately below)','hidden files/directories','__pycache__/'],inventory_outputs={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='VALIDATED.json'})
    (out/'VALIDATED.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='inventory_outputs'},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent)
    main(parser.parse_args().root.resolve())
