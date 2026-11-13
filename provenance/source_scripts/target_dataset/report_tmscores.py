"""Summarize the new TM-scores and draw publication heatmaps.

Run after create_tmscore_matrix.py: python report_tmscores.py
No alignments are rerun. All statistics exclude self-pairs and count each
unordered pair once. SD is sample SD (ddof=1); pair scores share proteins, so
these are descriptive statistics, not independent observations for inference.

Only pairs.csv and targets.csv from the current calculation are read. Legacy
matrices are not consulted. Each metric produces a standalone all-target figure
and a separate three-panel subgroup figure, in PDF and PNG.
Default: continuous colorblind-friendly cividis palette, large bold type.
Use --cmap Reds for white-to-red, or --font-size to adjust all text together.
Plots use a fixed 0–1 scale and mask the uncomputed diagonal in gray. Mean and
max of the two normalization scores are separate, explicitly labeled figures.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import tempfile

# Avoid plotting-library caches in home or substantial shared /tmp usage on HPC.
if 'MPLCONFIGDIR' not in os.environ:
    os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp(prefix='tmscore_mpl_', dir=os.environ.get('TMPDIR'))
os.environ.setdefault('XDG_CACHE_HOME', os.environ['MPLCONFIGDIR'])
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from create_tmscore_matrix import BASE_DIR, checksum, save_csv, atomic_write, summarize


def read_rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def load_results(directory):
    """Validate full pair coverage before producing publication summaries."""
    targets = read_rows(directory / 'targets.csv')
    pairs = read_rows(directory / 'pairs.csv')
    ids = [r['pdb_id'] for r in targets]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate target IDs')
    seen = set()
    for row in pairs:
        a, b = row['pdb_id_1'], row['pdb_id_2']
        key = tuple(sorted((a, b)))
        if a == b or a not in ids or b not in ids or key in seen:
            raise ValueError(f'Invalid or duplicate pair: {a}, {b}')
        seen.add(key)
        for field in ('length_class_1', 'length_class_2', 'n_ca_1', 'n_ca_2'):
            row[field] = int(row[field])
        for field in ('tm_norm_1', 'tm_norm_2', 'tm_mean', 'tm_max'):
            row[field] = float(row[field])
            if not np.isfinite(row[field]) or not 0 <= row[field] <= 1:
                raise ValueError(f'Invalid score: {key}')
        if not np.isclose(row['tm_mean'], (row['tm_norm_1']+row['tm_norm_2'])/2, atol=1e-12, rtol=0):
            raise ValueError(f'Incorrect pair mean: {key}')
        if row['tm_max'] != max(row['tm_norm_1'], row['tm_norm_2']):
            raise ValueError(f'Incorrect pair maximum: {key}')
    if len(pairs) != len(ids)*(len(ids)-1)//2 or not pairs:
        raise ValueError('Incomplete pair table')
    return targets, pairs


def plot_scores(targets, pairs, metric, output, palette='cividis', font_size=16, dpi=300):
    """Save all-target and subgroup figures with shared 0–1 limits.

    Every text element, including ticks, colorbar, and the diagonal note, uses
    bold type at least font_size points. Target order is preserved from the
    validated current calculation. The diagonal is masked, not assigned zero.
    """
    ids = [r['pdb_id'] for r in targets]
    index = {pid: i for i, pid in enumerate(ids)}
    matrix = np.full((len(ids), len(ids)), np.nan)
    for row in pairs:
        i, j = index[row['pdb_id_1']], index[row['pdb_id_2']]
        matrix[i, j] = matrix[j, i] = row[metric]
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': font_size, 'font.weight': 'bold',
        'axes.labelsize': font_size + 2, 'axes.labelweight': 'bold',
        'axes.titlesize': font_size + 4, 'axes.titleweight': 'bold',
        'xtick.labelsize': font_size, 'ytick.labelsize': font_size,
        'figure.titlesize': font_size + 6, 'figure.titleweight': 'bold',
        'pdf.fonttype': 42, 'svg.fonttype': 'none',
    })
    cmap = plt.get_cmap(palette).copy()
    cmap.set_bad('#D9D9D9')
    label = 'Mean' if metric == 'tm_mean' else 'Maximum'

    def panel(ax, size, letter=None):
        positions = [i for i, r in enumerate(targets)
                     if size is None or int(r['length_class']) == size]
        if not positions:
            raise ValueError(f'No targets in subgroup {size}')
        labels = [ids[i] for i in positions]
        image = ax.imshow(np.ma.masked_invalid(matrix[np.ix_(positions, positions)]),
                          vmin=0, vmax=1, cmap=cmap, interpolation='nearest')
        ax.set_xticks(range(len(labels)), labels, rotation=90,
                      fontsize=font_size, fontweight='bold')
        ax.set_yticks(range(len(labels)), labels, fontsize=font_size, fontweight='bold')
        ax.tick_params(axis='both', length=4, width=1.3, pad=5)
        count = len(labels)
        title = 'All targets' if size is None else f'~{size} aa subgroup'
        prefix = f'{letter}  ' if letter else ''
        ax.set_title(f'{prefix}{title}\n{count} proteins; {count*(count-1)//2} pairs', pad=16)
        if size is None:
            for i in range(1, len(targets)):
                if targets[i]['length_class'] != targets[i-1]['length_class']:
                    ax.axhline(i-.5, color='white', lw=1.8)
                    ax.axvline(i-.5, color='white', lw=1.8)
        return image

    scale = font_size / 16
    for layout in ('all', 'subgroups'):
        if layout == 'all':
            fig, ax = plt.subplots(figsize=(16*scale, 16*scale), layout='constrained')
            axes = [ax]
            image = panel(ax, None)
        else:
            fig, ax_array = plt.subplots(1, 3, figsize=(26*scale, 10*scale), layout='constrained')
            axes = list(ax_array)
            for ax, size, letter in zip(axes, (60, 100, 150), ('A', 'B', 'C')):
                image = panel(ax, size, letter)
        fig.suptitle(f'{label} of the two length-normalized TM-scores')
        colorbar = fig.colorbar(image, ax=axes, shrink=.85, pad=.02,
                               ticks=np.linspace(0, 1, 6), label=f'{label} TM-score')
        colorbar.ax.tick_params(labelsize=font_size, width=1.3)
        for tick in colorbar.ax.get_yticklabels():
            tick.set_fontweight('bold')
        colorbar.ax.yaxis.label.set_fontweight('bold')
        fig.supxlabel('Gray diagonal: self-comparisons excluded. Normalization: observed Cα counts.',
                      fontsize=font_size, fontweight='bold')
        stem = f'heatmap_all_{metric}' if layout == 'all' else f'heatmaps_subgroups_{metric}'
        for extension in ('png', 'pdf'):
            fig.savefig(output / f'{stem}.{extension}', dpi=dpi, facecolor='white')
        plt.close(fig)
    return matrix


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--results-dir', type=Path, default=BASE_DIR/'outputs/tmscores')
    parser.add_argument('--cmap', choices=('cividis', 'Reds'), default='cividis',
                        help='Continuous palette: colorblind-friendly cividis or white-to-red Reds')
    parser.add_argument('--font-size', type=float, default=16, help='Minimum font size in points; all fonts bold')
    parser.add_argument('--dpi', type=int, default=300, help='PNG resolution; PDF keeps vector artwork')
    args = parser.parse_args(argv)
    if args.font_size < 12 or args.dpi < 72:
        parser.error('--font-size must be at least 12 points and --dpi at least 72')
    targets, pairs = load_results(args.results_dir)
    summary = summarize(pairs)
    atomic_write(args.results_dir/'summary.json', lambda f: json.dump(summary,f,indent=2,allow_nan=False))
    flat = []
    for group, stats in summary.items():
        for metric in ('tm_mean','tm_max'):
            data = stats[metric]
            flat.append(dict(group=group, metric=metric, n_pairs=stats['n_pairs'],
                             mean=data['mean'], median=data['median'], std=data['std'], std_ddof=1,
                             min=data['min'], max=data['max'], max_pair=' / '.join(data['max_pair'])))
    save_csv(args.results_dir/'summary_statistics.csv',list(flat[0]),flat)
    ids = [r['pdb_id'] for r in targets]
    for metric in ('tm_mean','tm_max'):
        matrix = plot_scores(targets, pairs, metric, args.results_dir, args.cmap, args.font_size, args.dpi)
        for size in (60,100,150):
            selected = [i for i,r in enumerate(targets) if int(r['length_class'])==size]
            columns = [ids[i] for i in selected]
            rows = [dict(pdb_id=ids[i], **{ids[j]: 'NA' if i==j else matrix[i,j] for j in selected}) for i in selected]
            save_csv(args.results_dir/f'{metric}_subgroup_{size}.csv',['pdb_id',*columns],rows)
    provenance = dict(script_sha256=checksum(Path(__file__)), matplotlib=matplotlib.__version__,
                      inputs={name:checksum(args.results_dir/name) for name in ('pairs.csv','targets.csv')},
                      std_ddof=1, heatmap_limits=[0,1], diagonal='masked; excluded from statistics',
                      palette=args.cmap, minimum_font_size=args.font_size, font_weight='bold',
                      png_dpi=args.dpi, layouts=['all', 'subgroups'],
                      summary_script_sha256=checksum(BASE_DIR/'create_tmscore_matrix.py'))
    atomic_write(args.results_dir/'report_provenance.json',lambda f:json.dump(provenance,f,indent=2))
    print(f'Saved subgroup matrices, summary statistics and separate all-target/subgroup heatmaps to {args.results_dir}')


if __name__ == '__main__':
    main()
