import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

AF2_SWAPPED = Path('AF2_TMs.csv')
AF2_SINGLE = Path('AF2_single_TMs.csv')
AF3_SWAPPED = Path('tm_align_full_results.csv')
AF3_SINGLE = Path('singleseq_af3_tmscores.csv')

OUT_PNG = Path('figure3A_query_similarity_swapped_vs_single_50pct_text.png')
OUT_PDF = Path('figure3A_query_similarity_swapped_vs_single_50pct_text.pdf')
OUT_SVG = Path('figure3A_query_similarity_swapped_vs_single_50pct_text.svg')


def read_float_column(path, col):
    vals = []
    with open(path, newline='') as f:
        for row in csv.DictReader(f):
            v = row.get(col, '')
            if v not in ('', None):
                vals.append(float(v))
    return np.asarray(vals, dtype=float)

# 1. AF2 swapped: query-side scores from AF2_TMs.csv
# AF2_A is the TM-score relative to the query protein.
af2_swapped = read_float_column(AF2_SWAPPED, 'AF2_A')

# 2. AF2 single sequence: AF2_A
# This file contains a leading unnamed index column plus AF2_A.
af2_single = read_float_column(AF2_SINGLE, 'AF2_A')

# 3. AF3 swapped: query-side scores from junk/swapped-MSA jobs.
af3_swapped_vals = []
with open(AF3_SWAPPED, newline='') as f:
    for row in csv.DictReader(f):
        if str(row['Job_Type']).lower() != 'junk':
            continue
        query_pdb = row['Job_Label'].split('_')[0]
        if str(row['Compared_Protein']) == query_pdb:
            af3_swapped_vals.append(float(row['TM_Score']))
af3_swapped = np.asarray(af3_swapped_vals, dtype=float)

# 4. AF3 single sequence: last column, TM_Score.
af3_single = read_float_column(AF3_SINGLE, 'TM_Score')

for name, vals in [
    ('AF2 swapped', af2_swapped),
    ('AF2 single sequence', af2_single),
    ('AF3 swapped', af3_swapped),
    ('AF3 single sequence', af3_single),
]:
    print(f'{name}: n={len(vals):,}, median={np.median(vals):.3f}, mean={np.mean(vals):.3f}, >=0.6={np.mean(vals>=0.6)*100:.1f}%')

fig, ax = plt.subplots(figsize=(9.2, 6.8))

# Match the query-similarity styling from Figure 2.
ax.axhspan(0.6, 1.02, facecolor='#ae0639', alpha=0.10, zorder=0)
ax.text(
    0.02, 0.965, 'Similar to Query',
    transform=ax.transAxes,
    ha='left', va='top',
    fontsize=24,
    fontweight='bold',
    color='#ae0639',
    alpha=1.0,
)

parts = ax.violinplot(
    [af2_swapped, af2_single, af3_swapped, af3_single],
    positions=[1, 2, 3, 4],
    widths=0.72,
    showmeans=False,
    showmedians=True,
    showextrema=True,
    bw_method='scott',
)

for body in parts['bodies']:
    body.set_facecolor('0.75')
    body.set_edgecolor('black')
    body.set_alpha(1.0)
    body.set_zorder(2)

for key in ['cmedians', 'cbars', 'cmins', 'cmaxes']:
    parts[key].set_color('black')
    parts[key].set_linewidth(1.5)
    parts[key].set_zorder(3)

ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(
    [
        f'AlphaFold2\nSwapped MSA\n(n={len(af2_swapped):,})',
        f'AlphaFold2\nSingle sequence\n(n={len(af2_single):,})',
        f'AlphaFold3\nSwapped MSA\n(n={len(af3_swapped):,})',
        f'AlphaFold3\nSingle sequence\n(n={len(af3_single):,})',
    ],
    fontsize=17,
)
#ax.tick_params(axis="x", labelrotation=45)

ax.set_ylabel('Query similarity (TM-score)', fontsize=21)
ax.set_ylim(0, 1.02)
ax.tick_params(axis='y', labelsize=18)
ax.grid(axis='y', linestyle='--', alpha=0.22, zorder=1)

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight')
fig.savefig(OUT_PDF, bbox_inches='tight')
fig.savefig(OUT_SVG, bbox_inches='tight')
print(OUT_PNG)
print(OUT_PDF)
print(OUT_SVG)
