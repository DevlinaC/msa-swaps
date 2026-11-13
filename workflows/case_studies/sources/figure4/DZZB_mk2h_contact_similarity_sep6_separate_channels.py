#!/usr/bin/env python3
"""DZZB/mk2h structural and MSA-derived contact-map similarity.

Pair universe: upper-triangle residue pairs with |i-j| >= 6.
Structural contacts: Cbeta-Cbeta <= 8 A (Calpha for Gly).
For each homodimer, intramolecular contacts are the OR of chains A and B;
intermolecular contacts are AB OR AB.T. Structural similarity is evaluated
with intra- and intermolecular contacts treated as distinct channels.

MSA-derived similarity uses the single 43x43 residue-pair probability map.
Pearson R is computed over eligible pairs. Jaccard uses the same rank-tolerant
top-L definition as the preceding analyses, with L = 43 and a wider set of
floor(1.5 L) = 64; score ties at cutoffs are retained.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import pearsonr

BASE = Path('.')
PDB_DZZB = BASE / 'DZZB_sm2_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb'
PDB_MK2H = BASE / 'mk2h_sm2_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_002.pdb'
PAIR_DZZB = BASE / 'DZZB_predicted_cb_contacts.txt'
PAIR_MK2H = BASE / 'mk2h_predicted_cb_contacts.txt'
OUT_CSV = BASE / 'DZZB_mk2h_contact_similarity_sep6_separate_channels_metrics.csv'
OUT_PNG = BASE / 'DZZB_mk2h_contact_similarity_sep6_separate_channels_table.png'
OUT_SVG = BASE / 'DZZB_mk2h_contact_similarity_sep6_separate_channels_table.svg'

MIN_SEQUENCE_SEPARATION = 6
STRUCTURAL_CONTACT_CUTOFF_A = 8.0
TOP_L_MULTIPLIER = 1.5


def parse_chain_contact_coords(path):
    residues = {}
    order = {}
    with open(path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if line[16] not in (' ', 'A'):
                continue
            ch = line[21].strip() or '_'
            resnum = int(line[22:26])
            icode = line[26]
            resn = line[17:20].strip()
            atom = line[12:16].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            key = (resnum, icode)
            residues.setdefault(ch, {})
            order.setdefault(ch, [])
            if key not in residues[ch]:
                residues[ch][key] = {'resn': resn, 'atoms': {}}
                order[ch].append(key)
            residues[ch][key]['atoms'][atom] = xyz
    coords = {}
    for ch in order:
        arr = []
        for key in order[ch]:
            entry = residues[ch][key]
            atom = 'CA' if entry['resn'] == 'GLY' else 'CB'
            if atom not in entry['atoms']:
                atom = 'CA'
            arr.append(entry['atoms'][atom])
        coords[ch] = np.vstack(arr)
    return coords


def eligible_mask(n):
    i, j = np.indices((n, n))
    return np.triu(np.abs(i-j) >= MIN_SEQUENCE_SEPARATION, k=1)


def structural_maps(path):
    c = parse_chain_contact_coords(path)
    if set(c) < {'A', 'B'}:
        raise ValueError(f'Expected chains A and B in {path}')
    A, B = c['A'], c['B']
    if len(A) != len(B):
        raise ValueError('A/B chain length mismatch')

    def contact(X, Y):
        d = np.linalg.norm(X[:, None, :] - Y[None, :, :], axis=2)
        return d <= STRUCTURAL_CONTACT_CUTOFF_A

    AA = contact(A, A)
    BB = contact(B, B)
    AB = contact(A, B)
    np.fill_diagonal(AA, False)
    np.fill_diagonal(BB, False)

    intra = AA | BB
    inter = AB | AB.T
    np.fill_diagonal(intra, False)
    np.fill_diagonal(inter, False)
    return intra, inter


def binary_metrics(A, B, mask):
    a = A[mask].astype(float)
    b = B[mask].astype(float)
    r = float(pearsonr(a, b).statistic) if np.std(a) and np.std(b) else float('nan')
    sa = a.astype(bool)
    sb = b.astype(bool)
    union = sa | sb
    j = float((sa & sb).sum() / union.sum()) if union.sum() else float('nan')
    return {
        'pearson_r': r,
        'jaccard': j,
        'DZZB_contacts': int(sa.sum()),
        'mk2h_contacts': int(sb.sum()),
        'shared': int((sa & sb).sum()),
        'union': int(union.sum()),
    }


def distinct_channel_metrics(mapsA, mapsB, mask):
    a = np.concatenate([mapsA[0][mask].astype(float), mapsA[1][mask].astype(float)])
    b = np.concatenate([mapsB[0][mask].astype(float), mapsB[1][mask].astype(float)])
    r = float(pearsonr(a, b).statistic) if np.std(a) and np.std(b) else float('nan')
    sa = a.astype(bool)
    sb = b.astype(bool)
    union = sa | sb
    j = float((sa & sb).sum() / union.sum()) if union.sum() else float('nan')
    return {
        'pearson_r': r,
        'jaccard': j,
        'DZZB_contacts': int(sa.sum()),
        'mk2h_contacts': int(sb.sum()),
        'shared': int((sa & sb).sum()),
        'union': int(union.sum()),
    }


def ranked_set_with_ties(scores, mask, k):
    vals = scores[mask]
    threshold = np.sort(vals)[::-1][k - 1]
    selected = mask & (scores >= threshold)
    return {tuple(x) for x in np.argwhere(selected)}, float(threshold)


def rank_tolerant_jaccard(a, b, mask, L):
    wide_k = int(np.floor(TOP_L_MULTIPLIER * L))
    S, s_thr = ranked_set_with_ties(a, mask, L)
    T, t_thr = ranked_set_with_ties(b, mask, L)
    Sw, sw_thr = ranked_set_with_ties(a, mask, wide_k)
    Tw, tw_thr = ranked_set_with_ties(b, mask, wide_k)
    credited = (S & Tw) | (T & Sw)
    denom = S | T
    strict = len(S & T) / len(denom)
    tolerant = len(credited) / len(denom)
    return {
        'tolerant_jaccard': float(tolerant),
        'strict_topL_jaccard': float(strict),
        'credited_count': len(credited),
        'topL_union_count': len(denom),
        'topL_DZZB_count_with_ties': len(S),
        'topL_mk2h_count_with_ties': len(T),
        'wide_DZZB_count_with_ties': len(Sw),
        'wide_mk2h_count_with_ties': len(Tw),
        'topL_DZZB_threshold': s_thr,
        'topL_mk2h_threshold': t_thr,
        'wide_DZZB_threshold': sw_thr,
        'wide_mk2h_threshold': tw_thr,
    }


def make_table(rows):
    fig, ax = plt.subplots(figsize=(10.8, 4.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    left = 0.03
    bottom = 0.17
    width = 0.94
    height = 0.72
    col_edges = [left, left + width * 0.50, left + width * 0.75, left + width]
    row_edges = [bottom, bottom + height / 3, bottom + 2 * height / 3, bottom + height]

    #left, bottom, width, height = 0.03, 0.17, 0.94, 0.72
    #col_edges = [left, left + 0.50 * width, left + 0.75 * width, left + width]
    nrows = len(rows) + 1
    #row_edges = np.linspace(bottom, bottom + height, nrows + 1)

    for rr in range(nrows):
        y0 = row_edges[nrows - 1 - rr]
        y1 = row_edges[nrows - rr]
        for c in range(3):
            if rr == 0 or c == 0:
                gray = 1.0
            else:
                val = rows[rr - 1][c]
                gray = max(0.0, 1.0 - float(val))
            ax.add_patch(Rectangle(
                (col_edges[c], y0),
                col_edges[c + 1] - col_edges[c],
                y1 - y0,
                facecolor=str(gray),
                edgecolor='black',
                linewidth=2.1,
            ))

    headers = ['Representation', 'Pearson R', 'Jaccard']
    yh = (row_edges[-2] + row_edges[-1]) / 2
    for c, text in enumerate(headers):
        ax.text(
            (col_edges[c] + col_edges[c + 1]) / 2,
            yh,
            text,
            ha='center',
            va='center',
            fontsize=24,
            fontweight='bold',
        )

    for ridx, row in enumerate(rows):
        y0 = row_edges[nrows - 2 - ridx]
        y1 = row_edges[nrows - 1 - ridx]
        y = (y0 + y1) / 2
        ax.text(
            (col_edges[0] + col_edges[1]) / 2,
            y,
            row[0],
            ha='center',
            va='center',
            fontsize=24,
            fontweight='bold',
        )
        for c, val in ((1, row[1]), (2, row[2])):
            gray = max(0.0, 1.0 - float(val))
            text_color = 'white' if gray < 0.48 else 'black'
            ax.text(
                (col_edges[c] + col_edges[c + 1]) / 2,
                y,
                f'{val:.2f}',
                ha='center',
                va='center',
                fontsize=24,
                fontweight='bold',
                color=text_color,
            )

    #ax.text(
    #    0.5,
    #    0.075,
    #    r'Eligible residue pairs: $|i-j|\geq6$; structural contact cutoff = 8 Å.',
    #    ha='center',
    #    fontsize=11.4,
    #)
    #ax.text(
    #    0.5,
    #    0.025,
    #    r'Structural intra/inter contacts treated as distinct channels; MSA-derived Jaccard: rank-tolerant top-$L$ ($L=43$).',
    #    ha='center',
    #    fontsize=11.1,
    #)

    fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(OUT_SVG, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    P_D = np.loadtxt(PAIR_DZZB)
    P_M = np.loadtxt(PAIR_MK2H)
    if P_D.shape != P_M.shape or P_D.ndim != 2 or P_D.shape[0] != P_D.shape[1]:
        raise ValueError('Contact matrices must be square and the same shape')

    n = P_D.shape[0]
    mapsD = structural_maps(PDB_DZZB)
    mapsM = structural_maps(PDB_MK2H)
    if mapsD[0].shape != (n, n) or mapsM[0].shape != (n, n):
        raise ValueError('PDB/matrix length mismatch')

    mask = eligible_mask(n)
    intra = binary_metrics(mapsD[0], mapsM[0], mask)
    inter = binary_metrics(mapsD[1], mapsM[1], mask)
    distinct = distinct_channel_metrics(mapsD, mapsM, mask)

    msa_r = float(pearsonr(P_D[mask], P_M[mask]).statistic)
    msa_j = rank_tolerant_jaccard(P_D, P_M, mask, n)

    metrics = {
        'length': n,
        'minimum_sequence_separation_inclusive': MIN_SEQUENCE_SEPARATION,
        'eligible_upper_triangle_pairs_per_channel': int(mask.sum()),
        'structural_contact_cutoff_A': STRUCTURAL_CONTACT_CUTOFF_A,
        'intra_only_Pearson_R': intra['pearson_r'],
        'intra_only_Jaccard': intra['jaccard'],
        'intra_only_DZZB_contacts': intra['DZZB_contacts'],
        'intra_only_mk2h_contacts': intra['mk2h_contacts'],
        'intra_only_shared': intra['shared'],
        'intra_only_union': intra['union'],
        'inter_only_Pearson_R': inter['pearson_r'],
        'inter_only_Jaccard': inter['jaccard'],
        'inter_only_DZZB_contacts': inter['DZZB_contacts'],
        'inter_only_mk2h_contacts': inter['mk2h_contacts'],
        'inter_only_shared': inter['shared'],
        'inter_only_union': inter['union'],
        'distinct_channels_Pearson_R': distinct['pearson_r'],
        'distinct_channels_Jaccard': distinct['jaccard'],
        'distinct_channels_DZZB_contacts': distinct['DZZB_contacts'],
        'distinct_channels_mk2h_contacts': distinct['mk2h_contacts'],
        'distinct_channels_shared': distinct['shared'],
        'distinct_channels_union': distinct['union'],
        'MSA_derived_Pearson_R': msa_r,
        'MSA_derived_top_L': n,
        'MSA_derived_rank_tolerant_Jaccard': msa_j['tolerant_jaccard'],
        **{f'MSA_derived_{k}': v for k, v in msa_j.items() if k != 'tolerant_jaccard'},
    }

    with open(OUT_CSV, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])
        for key, value in metrics.items():
            writer.writerow([key, value])

    rows = [
        ('Structural contacts*', distinct['pearson_r'], distinct['jaccard']),
        ('MSA contacts', msa_r, msa_j['tolerant_jaccard']),
    ]
    make_table(rows)

    print('DZZB vs mk2h, |i-j| >= 6')
    print('Eligible residue pairs per channel:', mask.sum())
    print('Intra only:', intra)
    print('Inter only:', inter)
    print('Distinct channels:', distinct)
    print('MSA Pearson R:', msa_r)
    print('MSA tolerant Jaccard:', msa_j['tolerant_jaccard'])
    print('MSA strict top-L Jaccard:', msa_j['strict_topL_jaccard'])
    print('Tolerant credited/union:', msa_j['credited_count'], msa_j['topL_union_count'])
    print('Outputs:', OUT_CSV, OUT_PNG, OUT_SVG)


if __name__ == '__main__':
    main()
