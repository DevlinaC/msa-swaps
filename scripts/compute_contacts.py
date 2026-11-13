#!/usr/bin/env python3
"""Recompute contact-overlap scores from a portable input manifest.

Inputs are existing GREMLIN tables and MSA Transformer matrices, not model
inference inputs. Output contains both the historical and standard measures.
"""
# Historical October 2026 demonstration; see docs/FIGURES.md for manuscript workflows.
import argparse
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import jaccard as distance_jaccard
from sklearn.metrics import jaccard_score
import gremlin_contacts as g
import msatr_contacts as m

ROOT = Path(__file__).resolve().parents[1]


def standard_scores(a, b):
    """Same Boolean contact indicators for independent standard calculations."""
    union = sorted(a | b)
    if not union:
        raise ValueError('No selected contacts in either input; similarity is undefined')
    u = np.array([pair in a for pair in union], dtype=bool)
    v = np.array([pair in b for pair in union], dtype=bool)
    manual = len(a & b) / len(union)
    sk = float(jaccard_score(u, v, zero_division=0))
    sp = float(1 - distance_jaccard(u, v))
    if not np.allclose([manual, sk, sp], manual, atol=1e-12, rtol=0):
        raise ValueError('Standard Jaccard implementations disagree')
    return dict(standard_jaccard=manual, sklearn_jaccard=sk,
                scipy_similarity=sp, selected_swapped=len(a), selected_native=len(b),
                intersection=len(a & b), union=len(union))


@lru_cache(None)
def gremlin_set(directory):
    path = directory / 'gremlin.csv'
    mask = g.compute_gap_mask(str(directory / 'alignment.fas'))
    if mask is None:
        raise ValueError(f'Alignment cannot supply a gap mask: {directory}')
    frame = g.parse_gremlin(str(path))
    vals = frame[['i', 'j', 'score']].to_numpy()
    idx = frame[['i', 'j']].to_numpy()
    if (not np.isfinite(vals).all() or not np.equal(idx, np.floor(idx)).all()
            or (idx < 0).any() or (idx >= len(mask)).any()):
        raise ValueError(f'Invalid GREMLIN values or alignment indices: {path}')
    return g.get_gremlin_contact_set(str(path), gap_mask=mask)


@lru_cache(None)
def matrix(directory):
    a = np.loadtxt(directory / 'msatr.csv', delimiter=',')
    if a.ndim != 2 or a.shape[0] != a.shape[1] or not np.isfinite(a).all():
        raise ValueError(f'Expected a finite square matrix: {directory}')
    return a


def calculate(examples):
    pairs = pd.read_csv(examples / 'pairs.csv')
    if pairs['pair'].duplicated().any():
        raise ValueError('Duplicate manifest pairs')
    output = []
    for row in pairs.to_dict('records'):
        pair, a, b = row['pair'], row['protein_a'], row['protein_b']
        if pair != f'{a}_{b}_aligned':
            raise ValueError(f'Inconsistent pair identifiers: {pair}')
        for identifier in [pair,a,b]:
            if not identifier or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789_' for c in identifier):
                raise ValueError(f'Invalid path identifier: {identifier}')
        swapped = examples / 'swapped' / pair
        S = gremlin_set(swapped)
        raw = matrix(swapped)
        for reference, protein in [('A',a),('B',b)]:
            native = examples / 'native' / protein
            T = gremlin_set(native)
            selected, offset, best = S, 0, 0.0
            if reference == 'B':
                for trial_offset in range(-5,6):
                    trial = {(i+trial_offset,j+trial_offset) for i,j in S}
                    score = g.calculate_jaccard(trial,T)
                    if score > best:
                        selected,offset,best = trial,trial_offset,score
            output.append(dict(pair=pair,reference=reference,method='GREMLIN',
                reported_metric='standard_jaccard',reported_score=g.calculate_jaccard(selected,T),
                offset=offset,**standard_scores(selected,T)))
            other = matrix(native)
            n = min(len(raw),len(other))
            mask = np.triu(np.ones((n,n),bool),k=6)
            if mask.sum() < int(1.5*n):
                raise ValueError(f'Too few eligible contacts for historical top-1.5L selection: {pair}')
            x,y = m.crop_center(raw,n),m.crop_center(other,n)
            A,B = m.top(x,mask,n),m.top(y,mask,n)
            sets = [set(map(tuple,np.argwhere(z))) for z in [A,B]]
            output.append(dict(pair=pair,reference=reference,method='MSA Transformer',
                reported_metric='tolerant_contact_overlap',reported_score=m.tolerant_jaccard(raw,other),
                offset=0,**standard_scores(*sets)))
    return pd.DataFrame(output)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--examples',type=Path,default=ROOT/'data/examples')
    ap.add_argument('--output',type=Path,default=ROOT/'outputs/example_scores.csv')
    args = ap.parse_args()
    frame = calculate(args.examples.resolve())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    tmp = args.output.with_name(args.output.name+'.tmp')
    frame.to_csv(tmp,index=False)
    tmp.replace(args.output)
    print(f'Wrote {len(frame)} comparisons to {args.output}')

if __name__ == '__main__':
    main()
