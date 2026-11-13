#!/usr/bin/env python3
"""Replay collaborator case studies in isolation and verify supplied metric tables."""
import argparse
import csv
import hashlib
import importlib.metadata
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'workflows/case_studies/sources'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare_tables(expected, actual):
    """Require identical schema/labels and numerically equivalent scalar cells."""
    with expected.open(newline='') as f:
        a = list(csv.reader(f))
    with actual.open(newline='') as f:
        b = list(csv.reader(f))
    if not a or len(a) != len(b) or a[0] != b[0]:
        raise ValueError(f'Table schema/row count changed: {expected.name}')
    for row_number, (left, right) in enumerate(zip(a[1:], b[1:]), 2):
        if len(left) != len(right):
            raise ValueError(f'Column count changed: {expected.name}:{row_number}')
        for column, x, y in zip(a[0], left, right):
            if x == y:
                continue
            try:
                nx, ny = float(x), float(y)
            except ValueError:
                equal = False
            else:
                equal = (math.isnan(nx) and math.isnan(ny)) or math.isclose(
                    nx, ny, rel_tol=1e-10, abs_tol=1e-12)
            if not equal:
                raise ValueError(f'Metric mismatch: {expected.name}:{row_number}:{column}: {x!r} != {y!r}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['all', 'figure4', 'rfah_nusg_512'], default='all')
    parser.add_argument('--out-dir', type=Path, default=ROOT / 'outputs/case_studies')
    args = parser.parse_args()
    out = args.out_dir.resolve()
    # Refuse output into the imported source area, or any directory with prior outputs.
    if out == SOURCE or SOURCE in out.parents or out in SOURCE.parents:
        parser.error('Output must be separate from the imported source tree.')
    if out.exists() and any(out.iterdir()):
        parser.error('Output directory is not empty; choose a fresh --out-dir.')
    with (ROOT / 'provenance/collaborator_import_20261113.csv').open() as f:
        records = {row['destination']: row for row in csv.DictReader(f)}
    groups = ['figure4', 'rfah_nusg_512'] if args.case == 'all' else [args.case]
    # Authenticate all files before executing any preserved source code.
    for group in groups:
        for path in (SOURCE / group).iterdir():
            row = records[str(path.relative_to(ROOT))]
            if sha(path) != row['sha256']:
                raise ValueError(f'Imported source changed: {path}')
    out.mkdir(parents=True, exist_ok=True)
    results = []
    env = dict(os.environ, MPLBACKEND='Agg')
    for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        env[key] = '1'
    # tempfile honors TMPDIR, including allocated /lscratch on Biowulf.
    with tempfile.TemporaryDirectory(prefix='msa-case-studies-') as tmp:
        env['MPLCONFIGDIR'] = str(Path(tmp) / 'matplotlib')
        for group in groups:
            source = SOURCE / group
            work = Path(tmp) / group
            shutil.copytree(source, work)
            # Remove copied expected outputs so a script that fails to write cannot pass.
            for path in work.glob('*.csv'):
                path.unlink()
            scripts = sorted(work.glob('*.py'))
            if len(scripts) != (4 if group == 'figure4' else 1):
                raise ValueError(f'Unexpected script set in {group}')
            target = out / group
            target.mkdir()
            for script in scripts:
                print(f'Running {group}/{script.name}', flush=True)
                with (target / (script.stem + '.log')).open('w') as log:
                    subprocess.run([sys.executable, str(script)], cwd=work, env=env,
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
            for expected in sorted(source.glob('*.csv')):
                actual = work / expected.name
                compare_tables(expected, actual)
                results.append({'case': group, 'table': expected.name, 'status': 'PASS',
                                'expected_sha256': sha(expected), 'replayed_sha256': sha(actual)})
            for path in work.iterdir():
                if path.suffix in ('.csv', '.png', '.svg'):
                    shutil.copy2(path, target / path.name)
    if len(results) != sum(4 if group == 'figure4' else 1 for group in groups):
        raise ValueError('Unexpected metric-table count')
    with (out / 'validation.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    with (out / 'environment.txt').open('w') as f:
        f.write(f'python={sys.version}\nrunner_sha256={sha(Path(__file__))}\n')
        f.write(f'import_manifest_sha256={sha(ROOT / "provenance/collaborator_import_20261113.csv")}\n')
        for package in ['numpy', 'scipy', 'matplotlib', 'biopython']:
            f.write(f'{package}={importlib.metadata.version(package)}\n')
    print(f'PASS: {len(results)} supplied metric tables reproduced; outputs: {out}')


if __name__ == '__main__':
    main()
