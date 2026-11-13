"""Bounded inventory of item 1 structures; run under Slurm on Biowulf."""
import argparse
import csv
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--shared', type=Path, required=True)
    p.add_argument('--targets', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--expected', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    targets = list(csv.DictReader(a.targets.open()))
    records = []
    for t in targets:
        query = t['pdb_id'].lower()
        for donor in [None] + [x['pdb_id'].lower() for x in targets if x['length_class'] == t['length_class'] and x != t]:
            condition = 'native' if donor is None else 'swapped'
            name = query if donor is None else f'{query}_{donor}_aligned'
            for method in ('AF2', 'AF3'):
                candidates = []
                if method == 'AF2':
                    folder = a.shared/f"targets_{t['length_class']}"/query/f'{name}_models_3'
                    for f in sorted(folder.glob('*unrelaxed_rank_001_*.pdb')):
                        if f.name.startswith('.'):
                            continue
                        config = folder/'config.json'
                        settings = json.loads(config.read_text()) if config.exists() else {}
                        candidates.append(dict(path=str(f), ranking=str(f).replace('_unrelaxed_', '_scores_').replace('.pdb','.json'),
                            config=str(config), input_msa=str(folder/f'{name}.a3m'),
                            rank_by=settings.get('rank_by'), recycles=settings.get('num_recycles'),
                            run=str(folder), sha256=sha(f)))
                else:
                    base = a.shared/'af3_results_full'/('proper' if donor is None else 'junk')/query
                    if donor is not None:
                        base = base/name
                    for rank in sorted(base.glob('*/ranking_scores.csv')):
                        rows = list(csv.DictReader(rank.open()))
                        ordered = sorted(rows, key=lambda r: -float(r['ranking_score']))
                        if not ordered:
                            continue
                        best = ordered[0]
                        f = rank.parent/f"seed-{best['seed']}_sample-{best['sample']}"/'model.cif'
                        if f.exists():
                            candidates.append(dict(path=str(f), ranking=str(rank), run=str(rank.parent),
                                seed=best['seed'], sample=best['sample'], score=float(best['ranking_score']),
                                n_ranked=len(rows), ranking_ties=sum(float(r['ranking_score']) == float(best['ranking_score']) for r in rows),
                                sha256=sha(f)))
                records.append(dict(method=method, condition=condition, query=query, donor=donor,
                                    name=name, candidates=candidates, n_candidates=len(candidates)))
    (a.out/'structures.json').write_text(json.dumps(records, indent=2))
    sources = json.loads(a.expected.read_text())
    checks = []
    for rel, expected in sources.items():
        f = a.shared/rel
        actual = sha(f) if f.is_file() else None
        checks.append(dict(relative_path=rel, expected=expected, actual=actual,
                           status='ok' if actual == expected else 'missing' if actual is None else 'different'))
    (a.out/'source_checks.json').write_text(json.dumps(checks, indent=2))
    versions = {}
    for pkg in ('numpy', 'pandas', 'scipy', 'matplotlib', 'Bio'):
        try:
            m = importlib.import_module(pkg)
            versions[pkg] = m.__version__
        except Exception as e:
            versions[pkg] = repr(e)
    summary = dict(host=platform.node(), job_id=os.getenv('SLURM_JOB_ID'), python=platform.python_version(), versions=versions,
        coverage={f'{m}_{c}': {str(n): sum(r['method']==m and r['condition']==c and r['n_candidates']==n for r in records)
                  for n in sorted(set(r['n_candidates'] for r in records))} for m in ('AF2','AF3') for c in ('native','swapped')},
        sources={s:sum(x['status']==s for x in checks) for s in ('ok','missing','different')})
    (a.out/'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
