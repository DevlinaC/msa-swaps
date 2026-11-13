#!/usr/bin/env python3
"""Reproduce a fixed random sample without changing source results.

Use existing contact outputs; no model inference. Requires numpy, pandas,
scipy and scikit-learn. Native GREMLIN alignments may be supplied explicitly;
the default is the af_sampling copy, whose provenance must be checked.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import random
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.spatial.distance import jaccard
from sklearn.metrics import jaccard_score

REPO = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-root', type=Path, required=True)
    ap.add_argument('--native-alignments', type=Path)
    ap.add_argument('--out-dir', type=Path, required=True)
    ap.add_argument('--seed', type=int, default=20261030)
    ap.add_argument('--per-size', type=int, default=10)
    args = ap.parse_args()
    root = args.data_root.resolve()
    native = args.native_alignments or root / 'original_proper_msas_af_sampling'
    sources = {}

    def track(path):
        path = Path(path).resolve()
        if str(path) not in sources:
            sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    gp = track(REPO / 'TM_Jaccard_codes/Jaccard_code/gremlin_jaccard_v2.py')
    mp = track(REPO / 'TM_Jaccard_codes/Jaccard_code/msa_transformer_jaccard_v3.py')
    g, m = load_module('gremlin_audit', gp), load_module('msatr_audit', mp)
    saved = {}
    for method, suffix in [('GREMLIN', 'gremlin'), ('MSAT', 'msatransformer')]:
        p = track(REPO / f'TM_Jaccard_codes/jaccard_results_{suffix}.csv')
        with p.open() as f:
            rows = list(csv.DictReader(f))
        saved[method] = {(r['Protein_A'], r['Protein_B']): r for r in rows}
        assert len(saved[method]) == len(rows) == 362
    assert saved['GREMLIN'].keys() == saved['MSAT'].keys()

    @lru_cache(None)
    def matrix(path):
        a = np.loadtxt(track(path), delimiter=',')
        assert a.ndim == 2 and a.shape[0] == a.shape[1] and np.isfinite(a).all()
        return a

    @lru_cache(None)
    def contacts(path, fas):
        path, fas = track(path), track(fas)
        mask = g.compute_gap_mask(str(fas))
        assert mask is not None, f'Invalid alignment: {fas}'
        df = g.parse_gremlin(str(path))
        assert np.isfinite(df[['i', 'j', 'score']].to_numpy()).all()
        assert df[['i', 'j']].min().min() >= 0
        assert df[['i', 'j']].max().max() < len(mask), path
        return g.get_gremlin_contact_set(str(path), gap_mask=mask)

    bins = {60: [], 100: [], 150: []}
    for a, b in sorted(saved['MSAT']):
        # Only read the first row to identify nominal native target size.
        p = root / f'original_proper_msas_af_sampling/{a}/full_{a}.csv'
        with p.open() as f:
            n = len(next(csv.reader(f)))
        size = 60 if n <= 60 else n
        bins[size].append((a, b))
    rng = random.Random(args.seed)
    sample = [(pair, 'random', size) for size, pairs in bins.items()
              for pair in rng.sample(pairs, args.per_size)]
    chosen = {pair for pair, _, _ in sample}
    fixes = [(a, '2mm4') for a in sorted(m.FIXED_2MM4_TARGETS)
             if (a, '2mm4') not in chosen]
    sample += [(pair, 'repair_supplement', 60) for pair in rng.sample(fixes, min(2, len(fixes)))]
    results = []

    def compare_sets(S, T):
        union = sorted(S | T)
        assert union, 'Empty union: requires an explicit undefined-score policy'
        u = np.array([p in S for p in union], dtype=bool)
        v = np.array([p in T for p in union], dtype=bool)
        return (len(S & T) / len(union), float(jaccard_score(u, v, zero_division=0)),
                float(1 - jaccard(u, v)))

    for (a, b), selection, size in sample:
        pair = f'{a}_{b}_aligned'
        fixed = b == '2mm4' and a in m.FIXED_2MM4_TARGETS
        swapped_fas = (root / f'2mm4_fixes/{pair}/{pair}.fas' if fixed else
                       root / f'original_junk_msas/{a}/{pair}/{pair}.fas')
        S = contacts(root / f'gremlin/results_gremlin_junk/full_{pair}.csv', swapped_fas)
        for ref, protein in [('A', a), ('B', b)]:
            T = contacts(root / f'gremlin/results_gremlin_proper/full_{protein}.csv',
                         native / protein / f'{protein}.fas')
            shifted, best, offset = S, 0.0, 0
            if ref == 'B':
                for off in range(-5, 6):
                    trial = {(i + off, j + off) for i, j in S}
                    score = g.calculate_jaccard(trial, T)
                    if score > best:
                        shifted, best, offset = trial, score, off
            manual, sk, sp = compare_sets(shifted, T)
            actual = g.calculate_jaccard(shifted, T)
            previous = saved['GREMLIN'][(a, b)]
            results.append(dict(method='GREMLIN', pair=pair, reference=ref,
                selection=selection, size=size, corrected=fixed,
                saved=float(previous[f'Jaccard_vs_Protein{ref}']), reproduced=actual,
                standard=manual, sklearn=sk, scipy_similarity=sp,
                contacts_swapped=len(shifted), contacts_native=len(T), intersection=len(shifted&T),
                union=len(shifted|T), offset=offset,
                saved_offset=int(previous['Best_Offset_ProteinB']) if ref=='B' else 0))
        swapped_path = (root / f'2mm4_fixes/{pair}/full_{pair}.csv' if fixed else
                        root / f'original_junk_msas_af_sampling/{a}/{pair}/full_{pair}.csv')
        raw = matrix(swapped_path)
        for ref, protein in [('A', a), ('B', b)]:
            other = matrix(root / f'original_proper_msas_af_sampling/{protein}/full_{protein}.csv')
            n = min(len(raw), len(other))
            x, y = m.crop_center(raw, n), m.crop_center(other, n)
            mask = np.triu(np.ones((n,n), bool), k=6)
            A, B = m.top(x, mask, n), m.top(y, mask, n)
            S, T = set(map(tuple,np.argwhere(A))), set(map(tuple,np.argwhere(B)))
            manual, sk, sp = compare_sets(S, T)
            # Independently select contacts by ranking and expand threshold ties.
            def ranked(z, k):
                eligible = [(float(z[i,j]), (int(i),int(j))) for i,j in np.argwhere(mask)]
                threshold = sorted(score for score,_ in eligible)[-k]
                return {ij for score,ij in eligible if score >= threshold}
            assert S == ranked(x,n) and T == ranked(y,n)
            Sw, Tw = ranked(x,int(1.5*n)), ranked(y,int(1.5*n))
            tolerant = len((S&Tw)|(T&Sw)) / len(S|T)
            actual = float(m.tolerant_jaccard(raw,other))
            assert abs(tolerant-actual) < 1e-12
            results.append(dict(method='MSAT', pair=pair, reference=ref,
                selection=selection, size=size, corrected=fixed,
                saved=float(saved['MSAT'][(a,b)][f'Jaccard_vs_Protein{ref}']), reproduced=actual,
                standard=manual, sklearn=sk, scipy_similarity=sp,
                contacts_swapped=len(S),contacts_native=len(T),intersection=len(S&T),union=len(S|T),
                offset='',saved_offset=''))
    df = pd.DataFrame(results)
    df['saved_abs_error'] = abs(df.saved-df.reproduced)
    df['library_max_error'] = np.maximum(abs(df.standard-df.sklearn),abs(df.standard-df.scipy_similarity))
    df['tolerant_minus_standard'] = df.reproduced-df.standard
    assert (df.library_max_error < 1e-12).all()
    summary=[]
    for (method,ref,selection), block in df.groupby(['method','reference','selection']):
        summary.append(dict(method=method,reference=ref,selection=selection,n=len(block),
            max_saved_error=float(block.saved_abs_error.max()),
            saved_mismatches=int((block.saved_abs_error>1e-10).sum()),
            max_library_error=float(block.library_max_error.max()),
            mean_standard=float(block.standard.mean()),mean_current=float(block.reproduced.mean()),
            mean_difference=float(block.tolerant_minus_standard.mean()),
            max_difference=float(block.tolerant_minus_standard.max()),
            n_difference_ge_005=int((block.tolerant_minus_standard>=0.05).sum())))
    out=args.out_dir; out.mkdir(parents=True,exist_ok=True)
    df.to_csv(out/'comparisons.csv',index=False)
    pd.DataFrame(summary).to_csv(out/'summary.csv',index=False)
    metadata=dict(seed=args.seed,per_size=args.per_size,python=sys.version,
        versions={p.__name__:p.__version__ for p in [np,pd,scipy,sklearn]},
        data_root=str(root),native_alignment_candidate=str(native),
        tolerance=1e-10,screening_difference=0.05,source_sha256=sources)
    (out/'provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(pd.DataFrame(summary).to_string(index=False))
    print('Offset mismatches:',sum(df[df.method=='GREMLIN'].offset != df[df.method=='GREMLIN'].saved_offset))
    print('Pairs checked:',len(sample))

if __name__ == '__main__':
    main()
