#!/usr/bin/env python3
"""BCCIP alpha/beta structural and MSA-derived contact similarity after MSA registration.

Registration
------------
The two supplied gapped MSA query strings have lengths 322 (BCCIP alpha)
and 314 (BCCIP beta), matching the corresponding PDB/contact-matrix coordinate
systems. Their ungapped sequences are identical (166 residues). Residues are
therefore registered by ordinal position in this shared ungapped MSA sequence:
the kth non-gap residue in the alpha string is paired with the kth non-gap
residue in the beta string.

Comparison
----------
Eligible pairs are upper-triangle pairs in the 166-residue registered sequence
with |i-j| >= 6.

Structural contacts:
    Cbeta-Cbeta <= 8.0 A; Calpha is used for Gly.
    Pearson R compares binary contact vectors over all eligible registered pairs.
    Jaccard is the exact set overlap of structural contacts.

MSA-derived contacts:
    Pairformer contact scores are sampled at the registered full-length indices.
    Pearson R compares continuous scores over all eligible registered pairs.
    Jaccard uses the manuscript rank-tolerant top-L definition with L=166 and
    top-floor(1.5L)=249; ties at rank cutoffs are retained.
"""

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import pearsonr

BASE = Path('./')
PDB_ALPHA = BASE / 'exp_deep_balanced_sweep_1c95b88606b4__alpha_full_model_4_seed_2.pdb'
PDB_BETA = BASE / 'exp_deep_balanced_sweep_1c95b88606b4__beta_full_model_3_seed_3.pdb'
PAIR_ALPHA = BASE / 'bccip_alpha_predicted_cb_contacts.txt'
PAIR_BETA = BASE / 'bccip_beta_predicted_cb_contacts.txt'

OUT_CSV = BASE / 'BCCIP_alpha_beta_contact_similarity_MSA_registered_sep6_metrics.csv'
OUT_PNG = BASE / 'BCCIP_alpha_beta_contact_similarity_MSA_registered_sep6_table.png'
OUT_SVG = BASE / 'BCCIP_alpha_beta_contact_similarity_MSA_registered_sep6_table.svg'

ALPHA_GAPPED = (
    '---------------------------------------------------------'
    'TLIFRREDLTEEEYKGLKEILEEIREEAPQLEE'
    '-------------------'
    'GE'
    '-----------'
    'EVVGFVLKIDLSTEKGKAILKEVEALIERLKEEGVEEEIVELFKEFLKDETKPVGKALKENTKVVPPEVAKPIVEKLKEELANP'
    '-------'
    'KLTYLYLLLSAEH'
    '--------------------------------------------------------------'
    'GSIKEIILEMGVIEGSDIFDGKGEYFETFMIDLS'
)

BETA_GAPPED = (
    '---------------------------------------------------------'
    'TLIFRREDLTEEEYKGLKEILEEIREEAPQLEE'
    '-------------------'
    'GE'
    '-----------'
    'EVVGFVLKIDLSTEKGKAILKEVEALIERLKEEGVEEEIVELFKEFLKDETKPVGKALKENTKVVPPEVAKPIVEKLKEELANP'
    '-------'
    'KLTYLYLLLSAEH'
    '------------------------------------------------------'
    'GSIKEIILEMGVIEGSDIFDGKGEYFETFMIDLS'
)

MIN_SEQUENCE_SEPARATION = 6
STRUCTURAL_CONTACT_CUTOFF_A = 8.0
TOP_L_MULTIPLIER = 1.5


def parse_pdb_contact_coordinates(path: Path):
    """Return ordered Cbeta coordinates (Calpha for Gly) for the first chain."""
    residues = {}
    order = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith('ATOM'):
                continue
            if line[16] not in (' ', 'A'):
                continue
            chain = line[21].strip() or '_'
            # These predictions contain one chain; use the first encountered chain.
            if order and chain != order[0][0]:
                continue
            resnum = int(line[22:26])
            icode = line[26]
            resname = line[17:20].strip()
            atom = line[12:16].strip()
            key = (chain, resnum, icode)
            if key not in residues:
                residues[key] = {'resname': resname, 'atoms': {}}
                order.append(key)
            residues[key]['atoms'][atom] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )

    coords = []
    for key in order:
        entry = residues[key]
        atom_name = 'CA' if entry['resname'] == 'GLY' else 'CB'
        if atom_name not in entry['atoms']:
            atom_name = 'CA'
        if atom_name not in entry['atoms']:
            raise ValueError(f'Missing contact coordinate for {key} in {path}')
        coords.append(entry['atoms'][atom_name])
    return np.vstack(coords)


def build_registration(alpha_gapped: str, beta_gapped: str):
    """Map shared ungapped-sequence positions to zero-based full-length indices."""
    alpha_ungapped = alpha_gapped.replace('-', '')
    beta_ungapped = beta_gapped.replace('-', '')
    if alpha_ungapped != beta_ungapped:
        raise ValueError('The two gapped MSA strings do not encode the same ungapped sequence.')

    alpha_idx = np.array([i for i, aa in enumerate(alpha_gapped) if aa != '-'], dtype=int)
    beta_idx = np.array([i for i, aa in enumerate(beta_gapped) if aa != '-'], dtype=int)
    if len(alpha_idx) != len(beta_idx):
        raise ValueError('Different numbers of non-gap residues in registration strings.')

    return alpha_ungapped, alpha_idx, beta_idx


def registered_eligible_mask(n_registered: int):
    """Upper triangle in registered sequence coordinates, |i-j| >= 6."""
    i, j = np.indices((n_registered, n_registered))
    return np.triu(np.abs(i - j) >= MIN_SEQUENCE_SEPARATION, k=1)


def structural_contact_map(coords_registered: np.ndarray):
    delta = coords_registered[:, None, :] - coords_registered[None, :, :]
    distances = np.sqrt(np.sum(delta * delta, axis=2))
    contacts = distances <= STRUCTURAL_CONTACT_CUTOFF_A
    np.fill_diagonal(contacts, False)
    return contacts


def pair_set(binary_map: np.ndarray, mask: np.ndarray):
    return {tuple(x) for x in np.argwhere(binary_map & mask)}


def exact_jaccard(a: set, b: set):
    union = a | b
    return len(a & b) / len(union) if union else float('nan')


def ranked_set_with_ties(scores: np.ndarray, mask: np.ndarray, k: int):
    vals = scores[mask]
    if not (1 <= k <= vals.size):
        raise ValueError(f'Invalid top-k={k}; eligible pair count={vals.size}')
    threshold = np.sort(vals)[::-1][k - 1]
    selected = mask & (scores >= threshold)
    return {tuple(x) for x in np.argwhere(selected)}, float(threshold)


def rank_tolerant_jaccard(a: np.ndarray, b: np.ndarray, mask: np.ndarray, L: int):
    wide_k = int(np.floor(TOP_L_MULTIPLIER * L))
    S, s_thr = ranked_set_with_ties(a, mask, L)
    T, t_thr = ranked_set_with_ties(b, mask, L)
    S_w, sw_thr = ranked_set_with_ties(a, mask, wide_k)
    T_w, tw_thr = ranked_set_with_ties(b, mask, wide_k)

    credited = (S & T_w) | (T & S_w)
    denom = S | T
    strict = len(S & T) / len(denom) if denom else float('nan')
    tolerant = len(credited) / len(denom) if denom else float('nan')
    return {
        'tolerant_jaccard': float(tolerant),
        'strict_topL_jaccard': float(strict),
        'credited_count': len(credited),
        'topL_union_count': len(denom),
        'topL_alpha_count_with_ties': len(S),
        'topL_beta_count_with_ties': len(T),
        'wide_alpha_count_with_ties': len(S_w),
        'wide_beta_count_with_ties': len(T_w),
        'topL_alpha_threshold': s_thr,
        'topL_beta_threshold': t_thr,
        'wide_alpha_threshold': sw_thr,
        'wide_beta_threshold': tw_thr,
    }


def make_table_figure(struct_r, struct_j, pair_r, pair_j, L):
    fig, ax = plt.subplots(figsize=(10.8, 4.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')

    left, bottom, width, height = 0.03, 0.17, 0.94, 0.72
    col_edges = [left, left + 0.50 * width, left + 0.75 * width, left + width]
    row_edges = [bottom, bottom + height / 3, bottom + 2 * height / 3, bottom + height]

    fills = {
        (0, 0): 1.0, (0, 1): 1.0, (0, 2): 1.0,
        (1, 0): 1.0, (1, 1): max(0.0, 1.0 - struct_r), (1, 2): max(0.0, 1.0 - struct_j),
        (2, 0): 1.0, (2, 1): max(0.0, 1.0 - pair_r), (2, 2): max(0.0, 1.0 - pair_j),
    }

    for r in range(3):
        y0, y1 = row_edges[2-r], row_edges[3-r]
        for c in range(3):
            x0, x1 = col_edges[c], col_edges[c+1]
            gray = fills[(r, c)]
            ax.add_patch(Rectangle(
                (x0, y0), x1-x0, y1-y0,
                facecolor=str(gray), edgecolor='black', linewidth=2.3,
            ))

    headers = ['Representation', 'Pearson R', 'Jaccard']
    for c, txt in enumerate(headers):
        ax.text(
            (col_edges[c]+col_edges[c+1])/2,
            (row_edges[2]+row_edges[3])/2,
            txt, ha='center', va='center', fontsize=24, fontweight='bold',
        )

    rows = [
        ('Structural contacts', struct_r, struct_j),
        ('MSA contacts', pair_r, pair_j),
    ]
    y_centers = [(row_edges[1]+row_edges[2])/2, (row_edges[0]+row_edges[1])/2]
    for ridx, (label, rval, jval) in enumerate(rows):
        y = y_centers[ridx]
        ax.text(
            (col_edges[0]+col_edges[1])/2, y, label,
            ha='center', va='center', fontsize=24, fontweight='bold',
        )
        for c, value in ((1, rval), (2, jval)):
            bg_gray = fills[(ridx+1, c)]
            color = 'white' if bg_gray < 0.48 else 'black'
            ax.text(
                (col_edges[c]+col_edges[c+1])/2, y, f'{value:.2f}',
                ha='center', va='center', fontsize=27, fontweight='bold', color=color,
            )

    ax.text(
        0.5, 0.07,
        rf'Registered through shared {L}-residue MSA sequence; eligible pairs: $|i-j|\geq6$; structural cutoff = 8 Å.',
        ha='center', fontsize=11.2,
    )
    ax.text(
        0.5, 0.025,
        rf'Structural Jaccard: exact overlap; MSA-derived Jaccard: rank-tolerant top-$L$ ($L={L}$).',
        ha='center', fontsize=11.2,
    )

    fig.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(OUT_SVG, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    pair_alpha = np.loadtxt(PAIR_ALPHA)
    pair_beta = np.loadtxt(PAIR_BETA)
    coords_alpha = parse_pdb_contact_coordinates(PDB_ALPHA)
    coords_beta = parse_pdb_contact_coordinates(PDB_BETA)

    if pair_alpha.shape != (len(ALPHA_GAPPED), len(ALPHA_GAPPED)):
        raise ValueError(f'Alpha contact matrix shape {pair_alpha.shape} does not match gapped registration length {len(ALPHA_GAPPED)}')
    if pair_beta.shape != (len(BETA_GAPPED), len(BETA_GAPPED)):
        raise ValueError(f'Beta contact matrix shape {pair_beta.shape} does not match gapped registration length {len(BETA_GAPPED)}')
    if len(coords_alpha) != len(ALPHA_GAPPED):
        raise ValueError(f'Alpha PDB length {len(coords_alpha)} does not match gapped registration length {len(ALPHA_GAPPED)}')
    if len(coords_beta) != len(BETA_GAPPED):
        raise ValueError(f'Beta PDB length {len(coords_beta)} does not match gapped registration length {len(BETA_GAPPED)}')

    shared_seq, alpha_idx, beta_idx = build_registration(ALPHA_GAPPED, BETA_GAPPED)
    L = len(shared_seq)
    mask = registered_eligible_mask(L)

    # Reindex structures and Pairformer matrices into the common 166-residue frame.
    alpha_coords_reg = coords_alpha[alpha_idx]
    beta_coords_reg = coords_beta[beta_idx]
    alpha_pair_reg = pair_alpha[np.ix_(alpha_idx, alpha_idx)]
    beta_pair_reg = pair_beta[np.ix_(beta_idx, beta_idx)]

    alpha_struct = structural_contact_map(alpha_coords_reg)
    beta_struct = structural_contact_map(beta_coords_reg)

    struct_r = float(pearsonr(alpha_struct[mask].astype(float), beta_struct[mask].astype(float)).statistic)
    alpha_set = pair_set(alpha_struct, mask)
    beta_set = pair_set(beta_struct, mask)
    struct_j = float(exact_jaccard(alpha_set, beta_set))

    pair_r = float(pearsonr(alpha_pair_reg[mask], beta_pair_reg[mask]).statistic)
    pair_j_details = rank_tolerant_jaccard(alpha_pair_reg, beta_pair_reg, mask, L)
    pair_j = pair_j_details['tolerant_jaccard']

    metrics = {
        'registration': 'kth non-gap residue in alpha MSA string paired to kth non-gap residue in beta MSA string',
        'alpha_full_length': len(ALPHA_GAPPED),
        'beta_full_length': len(BETA_GAPPED),
        'registered_shared_length_L': L,
        'registered_shared_sequence': shared_seq,
        'minimum_registered_sequence_separation_inclusive': MIN_SEQUENCE_SEPARATION,
        'eligible_upper_triangle_registered_pairs': int(mask.sum()),
        'structural_contact_cutoff_A': STRUCTURAL_CONTACT_CUTOFF_A,
        'alpha_structural_contact_count': len(alpha_set),
        'beta_structural_contact_count': len(beta_set),
        'shared_structural_contact_count': len(alpha_set & beta_set),
        'union_structural_contact_count': len(alpha_set | beta_set),
        'structural_Pearson_R': struct_r,
        'structural_Jaccard': struct_j,
        'MSA_derived_Pearson_R': pair_r,
        'MSA_derived_top_L': L,
        'MSA_derived_rank_tolerant_Jaccard': pair_j,
        **{f'MSA_derived_{k}': v for k, v in pair_j_details.items() if k != 'tolerant_jaccard'},
    }

    with OUT_CSV.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['metric', 'value'])
        for k, v in metrics.items():
            writer.writerow([k, v])

    make_table_figure(struct_r, struct_j, pair_r, pair_j, L)

    print('BCCIP alpha vs beta, registered through supplied gapped MSA sequence')
    print('Full lengths alpha/beta:', len(ALPHA_GAPPED), len(BETA_GAPPED))
    print('Shared registered sequence length L:', L)
    print('Eligible registered pairs:', int(mask.sum()))
    print(f'Structural contacts alpha={len(alpha_set)}, beta={len(beta_set)}, shared={len(alpha_set & beta_set)}, union={len(alpha_set | beta_set)}')
    print(f'Structural Pearson R = {struct_r:.9f}')
    print(f'Structural exact Jaccard = {struct_j:.9f}')
    print(f'MSA-derived Pearson R = {pair_r:.9f}')
    print(f'MSA-derived rank-tolerant Jaccard = {pair_j:.9f}')
    print(f'MSA-derived strict top-L Jaccard = {pair_j_details["strict_topL_jaccard"]:.9f}')
    print(f'Credited/union for tolerant Jaccard = {pair_j_details["credited_count"]}/{pair_j_details["topL_union_count"]}')
    print('Outputs:', OUT_CSV, OUT_PNG, OUT_SVG)


if __name__ == '__main__':
    main()
