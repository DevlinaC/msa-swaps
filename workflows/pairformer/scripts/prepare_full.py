#!/usr/bin/env python3
"""Audit all 396 inputs, freeze deterministic 512-row subsets and donor mappings."""
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent / 'pilot_lib'))
from common import require, read_fasta, sha256, fingerprint, fasta_bytes, atomic_bytes, write_json


def select_indices(rows, cap=512):
    """Greedily maximize mean Hamming distance; preserve query, break ties by row index."""
    n = len(rows)
    if n <= cap:
        return np.arange(n, dtype=int)
    a = np.frombuffer(''.join(s for _, s in rows).encode('ascii'), dtype=np.uint8).reshape(n, -1)
    selected = [0]
    used = np.zeros(n, dtype=bool)
    used[0] = True
    distance_sum = np.zeros(n, dtype=np.int64)
    while len(selected) < cap:
        distance_sum += (a != a[selected[-1]]).sum(axis=1)
        scores = np.where(used, -1, distance_sum)
        index = int(scores.argmax())
        selected.append(index)
        used[index] = True
    return np.sort(selected)


def mapping(swapped, donor):
    reference = {}
    for h, s in donor[1:]:
        require(h not in reference or reference[h] == s, 'Conflicting duplicate donor header')
        reference[h] = s
    shared = {h: (s, reference[h]) for h, s in swapped[1:] if h in reference}
    if len(shared) < 3:
        return dict(status='undefined', reason='fewer than three shared homologs', shared_homologs=len(shared))
    ls, ln = len(swapped[0][1]), len(donor[0][1])
    sm = np.frombuffer(''.join(s for s, _ in shared.values()).encode(), dtype=np.uint8).reshape(-1, ls)
    dm = np.frombuffer(''.join(d for _, d in shared.values()).encode(), dtype=np.uint8).reshape(-1, ln)
    candidates = []
    for offset in range(-12, 13):
        lo, hi = max(0, -offset), min(ls, ln-offset)
        if hi-lo < .8*min(ls, ln):
            continue
        a, b = sm[:, lo:hi], dm[:, lo+offset:hi+offset]
        informative = (a != ord('-')) | (b != ord('-'))
        total = int(informative.sum())
        identity = float(((a == b) & informative).sum())/total if total else 0
        candidates.append((identity, offset, lo, hi))
    candidates.sort(reverse=True)
    best = candidates[0]
    unique = sum(abs(c[0]-best[0]) < 1e-12 for c in candidates) == 1
    return dict(status='ok' if unique and best[0] >= .99 else 'undefined',
                offset=best[1], start=best[2], stop=best[3], identity=best[0],
                shared_homologs=len(shared), unique=unique,
                convention='donor_column = swap_column + offset')


def prepare(root, source):
    inventory = json.loads((root/'source_inventory.json').read_text())
    config = json.loads((root/'model_config.json').read_text())
    require(config['max_depth']==512 and config['query_biasing'] is False, 'Wrong manuscript configuration')
    settings = dict(cap=config['max_depth'], method='greedy maximum mean Hamming distance',
                    tie_break='earliest source row', query_row=0, hhfilter=False,
                    output_order='original source order', crop=False,
                    historical_subsets='not independently recovered; source MSAs match historical run files')
    expected = fingerprint(dict(inventory=inventory, settings=settings, code=sha256(Path(__file__))))
    inputs = root/'inputs'; inputs.mkdir(parents=True, exist_ok=True)
    manifest_path = inputs/'manifest.json'
    if manifest_path.exists():
        prior = json.loads(manifest_path.read_text())
        require(prior['preparation_fingerprint'] == expected, 'Changed preparation; use a new run directory')
        for item in prior['inputs']:
            require(sha256(root/item['frozen_path']) == item['selected_sha256'], 'Corrupt frozen input')
        print('Validated prepared inputs: 396; no repeat preparation', flush=True)
        return
    native = {}; prepared = []; failures = []; began = time.monotonic()
    for spec in inventory['inputs']:
        start = time.monotonic()
        try:
            path = source/spec['path']
            require(sha256(path) == spec['source_sha256'], f'Source differs from local audit: {path}')
            rows = read_fasta(path)
            require(len(rows) == spec['depth'] and len(rows[0][1]) == spec['length'], 'Input shape mismatch')
            require(rows[0][1] == spec['query_sequence'], 'Unexpected query')
            if spec['kind'] == 'proper':
                native[spec['query']] = rows
                donor_map = None
            else:
                require(rows[0][1] == native[spec['query']][0][1], 'Swap/proper query mismatch')
                donor_map = mapping(rows, native[spec['donor']])
            indices = select_indices(rows, cap=settings['cap'])
            selected = [rows[i] for i in indices]
            require(indices[0] == 0 and selected[0] == rows[0], 'Query not preserved')
            dest = inputs/(spec['id']+'.fas')
            data = fasta_bytes(selected)
            if dest.exists():
                require(dest.read_bytes() == data, 'Partial-run input differs; use new run directory')
            else:
                atomic_bytes(dest, data)
            record = dict(spec, selected_depth=len(indices), selected_row_indices=indices.tolist(),
                          selected_sha256=sha256(dest), frozen_path=str(dest.relative_to(root)),
                          donor_mapping=donor_map, preparation_seconds=time.monotonic()-start)
            write_json(inputs/(spec['id']+'.json'), record)
            prepared.append(record)
            print(f"{len(prepared)}/396 {spec['id']}: {len(rows)} -> {len(selected)} rows; {time.monotonic()-start:.2f}s", flush=True)
        except Exception as e:
            failures.append(dict(id=spec['id'], error=repr(e)))
            print(f"FAILED {spec['id']}: {e}", file=sys.stderr, flush=True)
    write_json(root/'preparation_report.json', dict(prepared=len(prepared), failures=failures, seconds=time.monotonic()-began,
               undefined_mappings=[s['id'] for s in prepared if s['kind']=='swap' and s['donor_mapping']['status']!='ok']))
    require(not failures and len(prepared)==396, 'Input audit/preparation failed; GPU job must not start')
    manifest = dict(preparation_fingerprint=expected, selection=settings, inputs=prepared,
                    inventory_sha256=sha256(root/'source_inventory.json'))
    manifest['fingerprint'] = fingerprint(manifest)
    write_json(manifest_path, manifest)
    print(f'Preparation complete: {time.monotonic()-began:.2f}s', flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True); p.add_argument('--source', type=Path, required=True)
    a=p.parse_args(); prepare(a.root.resolve(), a.source.resolve())
