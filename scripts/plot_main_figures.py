#!/usr/bin/env python3
"""Replot supplied Figure 2 and Figure 3A panels; no statistical tests or inference.

Figure 2 retains the delivered historical contact selection, not revised S2.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/publication/Figures_2_3'
SOURCES = ROOT / 'provenance/source_scripts/collaborator_20261113'
SCRIPTS = {
    '2': ['figure2_query_similarity_violins_bg90pct_transparency',
          'figure2_donor_similarity_violins_bg90pct_transparency'],
    '3': ['figure3A_query_similarity_swapped_vs_single_50pct_text'],
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure', choices=['2', '3', 'all'], default='all')
    parser.add_argument('--out-dir', type=Path, default=ROOT / 'outputs/main_figures')
    args = parser.parse_args()
    out = args.out_dir.resolve()
    if out.exists() and any(out.iterdir()):
        parser.error('Choose an empty --out-dir; existing results will not be overwritten.')
    for source in (DATA, SOURCES):
        if out == source or source in out.parents or out in source.parents:
            parser.error('Outputs must be separate from source directories.')
    manifest = ROOT / 'provenance/figures_2_3_sources.csv'
    with manifest.open() as f:
        records = list(csv.DictReader(f))
    for row in records:
        path = ROOT / row['destination']
        if path.stat().st_size != int(row['bytes']) or sha(path) != row['sha256']:
            raise ValueError(f'Changed supplied input/source: {path}')

    # Set thread bounds before loading numerical libraries.
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[key] = '1'
    import numpy as np
    import pandas as pd

    af2 = pd.read_csv(DATA / 'AF2_TMs.csv')
    af3 = pd.read_csv(DATA / 'tm_align_full_results.csv')
    swapped = af3[af3.Job_Type.str.lower().eq('junk')].copy()
    parts = swapped.Job_Label.str.split('_')
    swapped['query'] = parts.str[0]
    swapped['donor'] = parts.str[1]
    datasets = []
    for role, column in [('query', 'A'), ('donor', 'B')]:
        datasets.append(('2', 'AF2', role, 'TM-score', af2[f'AF2_{column}'].to_numpy(), 9050))
        selected = swapped[swapped.Compared_Protein.eq(swapped[role])]
        if selected.duplicated(['Job_Label', 'Seed_Number', 'Sample_Number']).any():
            raise ValueError('Duplicate AF3 model/reference records')
        if len(selected.groupby('Job_Label')) != 362 or not selected.groupby('Job_Label').size().eq(25).all():
            raise ValueError('AF3 swapped population mismatch')
        datasets.append(('2', 'AF3', role, 'TM-score', selected.TM_Score.to_numpy(), 9050))
        for suffix, method, metric in [('gremlin', 'GREMLIN', 'ordinary Jaccard'),
                                       ('msatransformer', 'MSA Transformer', 'rank-tolerant overlap')]:
            table = pd.read_csv(DATA / f'jaccard_results_{suffix}.csv')
            if table.duplicated(['Protein_A', 'Protein_B']).any():
                raise ValueError('Duplicate contact comparison keys')
            datasets.append(('2', method, role, metric, table[f'Jaccard_vs_Protein{column}'].to_numpy(), 362))
    for method, file, column in [('AF2', 'AF2_single_TMs.csv', 'AF2_A'),
                                 ('AF3', 'singleseq_af3_tmscores.csv', 'TM_Score')]:
        table = pd.read_csv(DATA / file)
        if method == 'AF3':
            if table.duplicated(['Target', 'Seed_Number', 'Sample_Number']).any():
                raise ValueError('Duplicate AF3 single-sequence keys')
            if table.Target.nunique() != 34 or not table.groupby('Target').size().eq(25).all():
                raise ValueError('AF3 single-sequence population mismatch')
        datasets.append(('3', method, 'single_sequence_query', 'TM-score', table[column].to_numpy(), 850))
    for figure, method, ref, metric, values, expected in list(datasets):
        if figure == '2' and method in ('AF2', 'AF3') and ref == 'query':
            datasets.append(('3', method, 'swapped_query', metric, values, expected))
    summaries = []
    for figure, method, ref, metric, values, expected in datasets:
        if len(values) != expected or not np.isfinite(values).all() or not ((values >= 0) & (values <= 1)).all():
            raise ValueError(f'Invalid plotting population: {method} {ref}')
        summaries.append(dict(figure=figure, method=method, reference_or_condition=ref,
                              metric=metric, n=len(values), median=float(np.median(values)), mean=float(np.mean(values))))

    figures = ['2', '3'] if args.figure == 'all' else [args.figure]
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, MPLBACKEND='Agg')
    exports = []
    with tempfile.TemporaryDirectory(prefix='msa-main-figures-') as tmp:
        work = Path(tmp)
        env['MPLCONFIGDIR'] = str(work / 'matplotlib')
        for path in DATA.glob('*.csv'):
            shutil.copy2(path, work / path.name)
        for figure in figures:
            target = out / f'Figure_{figure}'
            target.mkdir()
            if figure == '2':
                print('Figure 2: supplied historical contact tables; revised main-panel selection is still pending.', flush=True)
            for stem in SCRIPTS[figure]:
                script = SOURCES / f'Fig{figure}' / (stem + '.py')
                print(f'Replotting {stem}', flush=True)
                with (target / (stem + '.log')).open('w') as log:
                    subprocess.run([sys.executable, str(script)], cwd=work, env=env,
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
                for ext in (['png', 'pdf'] if figure == '2' else ['png', 'pdf', 'svg']):
                    path = work / f'{stem}.{ext}'
                    if not path.is_file() or not path.stat().st_size:
                        raise ValueError(f'Missing figure output: {path}')
                    shutil.copy2(path, target / path.name)
                    exports.append(dict(path=str((target / path.name).relative_to(out)), sha256=sha(path)))
    for filename, rows in [('plotting_populations.csv', [r for r in summaries if r['figure'] in figures]),
                           ('output_sha256.csv', exports)]:
        with (out / filename).open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    with (out / 'environment.txt').open('w') as f:
        f.write(f'python={sys.version}\nrunner_sha256={sha(Path(__file__))}\nsource_manifest_sha256={sha(manifest)}\n')
        for package in ('numpy', 'pandas', 'matplotlib'):
            f.write(f'{package}={importlib.metadata.version(package)}\n')
    print(f'PASS: verified plotting populations and {len(exports)} figure exports. No statistical tests run. Outputs: {out}')


if __name__ == '__main__':
    main()
