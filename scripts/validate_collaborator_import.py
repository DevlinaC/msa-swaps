#!/usr/bin/env python3
"""Validate imported source bytes and mutant sequence/confidence correspondence."""
from collections import defaultdict
import csv
import hashlib
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/publication/Figure_5'
AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),
              'ARNDCQEGHILKMFPSTWYV'))


def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    imported = read(ROOT / 'provenance/collaborator_import_20261113.csv')
    derived = read(ROOT / 'provenance/collaborator_derived_sha256.csv')
    for row in imported + derived:
        path = ROOT / row['destination']
        require(path.stat().st_size == int(row['bytes']), f'Size changed: {path}')
        require(hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'], f'Hash changed: {path}')
    models = read(DATA / 'models.csv')
    inputs = {r['construct']: r for r in read(DATA / 'inputs.csv')}
    require(len(models) == 100 and len(inputs) == 4, 'Expected 100 models and four inputs')
    require(len({r['model_id'] for r in models}) == 100, 'Duplicate model identifiers')
    confidence = defaultdict(dict)
    for row in read(DATA / 'residue_plddt.csv'):
        index = int(row['residue_index_1based'])
        require(index not in confidence[row['model_id']], 'Duplicate confidence residue')
        confidence[row['model_id']][index] = float(row['plddt'])
    require(set(confidence) == {r['model_id'] for r in models}, 'Confidence/model correspondence')
    for construct, source in inputs.items():
        lines = (ROOT / source['input_path']).read_text().splitlines()
        seq = next(line for line in lines if line and not line.startswith(('#', '>')))
        require(seq == source['query_sequence'] and len(seq) == int(source['length']), 'Input sequence mismatch')
        group = [r for r in models if r['construct'] == construct]
        require(len(group) == 25, 'Expected 25 predictions per construct')
        require({(int(r['model']), int(r['seed'])) for r in group} ==
                {(m, s) for m in range(1, 6) for s in range(5)}, 'Model/seed coverage')
        require({int(r['rank']) for r in group} == set(range(1, 26)), 'Rank coverage')
        for row in group:
            path = ROOT / row['model_path']
            require(hashlib.sha256(path.read_bytes()).hexdigest() == row['model_sha256'], 'Model digest mismatch')
            atoms = [line for line in path.read_text().splitlines()
                     if line.startswith('ATOM  ') and line[12:16].strip() == 'CA']
            require(''.join(AA[line[17:20]] for line in atoms) == seq, f'Coordinate sequence mismatch: {path}')
            values = confidence[row['model_id']]
            require(set(values) == set(range(1, len(seq) + 1)), 'Confidence residue coverage')
            require(all(0 <= x <= 100 for x in values.values()), 'Confidence range')
            require(all(math.isclose(values[i], float(line[60:66]), abs_tol=0.011)
                        for i, line in enumerate(atoms, 1)), f'PDB/confidence mismatch: {path}')
            require(math.isclose(sum(values.values()) / len(seq), float(row['mean_plddt_all_residues']),
                                 abs_tol=1e-10), 'Mean confidence mismatch')
            require(row['fold_assignment'] == 'unassigned' and
                    row['included_in_fold_percentage'] == 'not_evaluated', 'Unreviewed fold classification')
    print(f'PASS: {len(imported)} imported files; 100 models; four sequences; 8,250 residue confidence values.')


if __name__ == '__main__':
    main()
