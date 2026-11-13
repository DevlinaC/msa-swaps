#!/usr/bin/env python3
"""Compare RfaH/NusG 512-sequence subsampled MSA Pairformer contact maps
and representative AlphaFold structures using the Figure 4 definitions.

Figure 4 definitions
--------------------
* Sequence length: 162 residues.
* CTD: residues 110-162, inclusive (1-based); L = 53.
* Eligible pairs: upper triangle, |i-j| >= 6, with at least one residue in CTD.
* Structural contacts: C-beta <= 8 A; C-alpha is used for Gly.
* Structural similarity: Pearson R on binary contact vectors and exact Jaccard.
* Pairformer similarity: Pearson R on continuous Pairformer scores and
  rank-tolerant top-L Jaccard.
* Rank-tolerant Jaccard: each top-L contact from one map is credited when it
  appears within the top floor(1.5L) contacts of the other map; denominator is
  the union of the two top-L sets. Ties at rank cutoffs are retained.

The script reports all 2x2 RfaH-vs-NusG comparisons, the two natural paired
comparisons (contact map 1 with 1, 2 with 2; structure 1 with 1, 2 with 2),
and within-target replicate similarity.
"""

from pathlib import Path
import csv
import numpy as np
from scipy.stats import pearsonr
from Bio.PDB import MMCIFParser

BASE = Path('./')

PAIRFORMER = {
    'RfaH_1': BASE / 'RfaH_contacts_1.txt',
    'RfaH_2': BASE / 'RfaH_contacts_2.txt',
    'NusG_1': BASE / 'NusG_contacts_1.txt',
    'NusG_2': BASE / 'NusG_contacts_2.txt',
}

# The two representative structures supplied for each target are treated as
# replicate 1 and replicate 2 in the order provided.
STRUCTURES = {
    'RfaH_1': BASE / 'RfaH_512_sample_1_model.cif',
    'RfaH_2': BASE / 'RfaH_512_sample_2_model.cif',
    'NusG_1': BASE / 'NusG_512_sample_1_model.cif',
    'NusG_2': BASE / 'NusG_512_sample_2_model.cif',
}

OUT = BASE / 'RfaH_NusG_512_overlap_results.csv'

CTD_START = 110
CTD_END = 162
CTD_LENGTH = CTD_END - CTD_START + 1  # 53
MIN_SEP = 6
DIST_CUTOFF_A = 8.0
TOP_L = CTD_LENGTH
WIDE_MULTIPLIER = 1.5


def eligible_mask(n: int) -> np.ndarray:
    """Upper-triangle Figure 4 pair universe."""
    i, j = np.indices((n, n))
    ctd0 = CTD_START - 1
    ctd1 = CTD_END
    ctd_involving = ((i >= ctd0) & (i < ctd1)) | ((j >= ctd0) & (j < ctd1))
    separated = np.abs(i - j) >= MIN_SEP
    return np.triu(ctd_involving & separated, k=1)


def parse_cif_contact_coordinates(path: Path):
    """Return one contact coordinate/residue: C-beta, or C-alpha for Gly."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    model = next(structure.get_models())
    chain = next(model.get_chains())

    coords = []
    residue_numbers = []
    residue_names = []
    for residue in chain.get_residues():
        if residue.id[0] != ' ':
            continue
        atom_name = 'CA' if residue.resname == 'GLY' else 'CB'
        if atom_name not in residue:
            raise ValueError(f'Missing {atom_name} for {residue.resname} {residue.id[1]} in {path}')
        coords.append(np.asarray(residue[atom_name].coord, dtype=float))
        residue_numbers.append(residue.id[1])
        residue_names.append(residue.resname)

    return np.vstack(coords), residue_numbers, residue_names


def structural_contact_map(coords: np.ndarray) -> np.ndarray:
    delta = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt(np.sum(delta * delta, axis=2))
    contacts = dist <= DIST_CUTOFF_A
    np.fill_diagonal(contacts, False)
    return contacts


def pair_set(binary_map: np.ndarray, mask: np.ndarray):
    return {tuple(x) for x in np.argwhere(binary_map & mask)}


def exact_jaccard(a: set, b: set) -> float:
    union = a | b
    return len(a & b) / len(union) if union else float('nan')


def ranked_set_with_ties(scores: np.ndarray, mask: np.ndarray, k: int):
    vals = scores[mask]
    if not (1 <= k <= vals.size):
        raise ValueError(f'Invalid k={k}; only {vals.size} eligible pairs')
    threshold = np.sort(vals)[::-1][k - 1]
    selected = mask & (scores >= threshold)
    return {tuple(x) for x in np.argwhere(selected)}, float(threshold)


def rank_tolerant_jaccard(a: np.ndarray, b: np.ndarray, mask: np.ndarray, L: int):
    wide_k = int(np.floor(WIDE_MULTIPLIER * L))
    A, a_thr = ranked_set_with_ties(a, mask, L)
    B, b_thr = ranked_set_with_ties(b, mask, L)
    A_wide, aw_thr = ranked_set_with_ties(a, mask, wide_k)
    B_wide, bw_thr = ranked_set_with_ties(b, mask, wide_k)

    credited = (A & B_wide) | (B & A_wide)
    denom = A | B
    strict_shared = A & B

    return {
        'rank_tolerant_jaccard': len(credited) / len(denom) if denom else float('nan'),
        'strict_topL_jaccard': len(strict_shared) / len(denom) if denom else float('nan'),
        'credited_count': len(credited),
        'topL_union_count': len(denom),
        'topL_A_count_with_ties': len(A),
        'topL_B_count_with_ties': len(B),
        'wide_A_count_with_ties': len(A_wide),
        'wide_B_count_with_ties': len(B_wide),
        'topL_A_threshold': a_thr,
        'topL_B_threshold': b_thr,
        'wide_A_threshold': aw_thr,
        'wide_B_threshold': bw_thr,
    }


def compare_pairformer(name_a, a, name_b, b, mask):
    r = float(pearsonr(a[mask], b[mask]).statistic)
    rt = rank_tolerant_jaccard(a, b, mask, TOP_L)
    return {
        'representation': 'MSA Pairformer',
        'A': name_a,
        'B': name_b,
        'pearson_R': r,
        'jaccard': rt['rank_tolerant_jaccard'],
        'jaccard_type': 'rank-tolerant top-L',
        'strict_topL_jaccard': rt['strict_topL_jaccard'],
        'A_contact_count': rt['topL_A_count_with_ties'],
        'B_contact_count': rt['topL_B_count_with_ties'],
        'shared_or_credited_count': rt['credited_count'],
        'union_count': rt['topL_union_count'],
    }


def compare_structures(name_a, a, name_b, b, mask):
    va = a[mask].astype(float)
    vb = b[mask].astype(float)
    r = float(pearsonr(va, vb).statistic)
    A = pair_set(a, mask)
    B = pair_set(b, mask)
    return {
        'representation': 'AlphaFold structure',
        'A': name_a,
        'B': name_b,
        'pearson_R': r,
        'jaccard': exact_jaccard(A, B),
        'jaccard_type': 'exact structural contact-set',
        'strict_topL_jaccard': '',
        'A_contact_count': len(A),
        'B_contact_count': len(B),
        'shared_or_credited_count': len(A & B),
        'union_count': len(A | B),
    }


def main():
    pair = {k: np.loadtxt(v) for k, v in PAIRFORMER.items()}
    shapes = {k: x.shape for k, x in pair.items()}
    if len(set(shapes.values())) != 1:
        raise ValueError(f'Pairformer matrices differ in shape: {shapes}')
    n = next(iter(pair.values())).shape[0]
    if any(x.shape != (n, n) for x in pair.values()):
        raise ValueError('All Pairformer matrices must be square')

    mask = eligible_mask(n)

    struct = {}
    for k, path in STRUCTURES.items():
        coords, nums, names = parse_cif_contact_coordinates(path)
        if len(coords) != n:
            raise ValueError(f'{k}: structure length {len(coords)} != matrix length {n}')
        if nums != list(range(1, n + 1)):
            raise ValueError(f'{k}: expected sequential residue numbering 1..{n}')
        struct[k] = structural_contact_map(coords)

    rows = []

    # Cross-target 2x2 comparisons.
    for r in ('RfaH_1', 'RfaH_2'):
        for ng in ('NusG_1', 'NusG_2'):
            rows.append(compare_pairformer(r, pair[r], ng, pair[ng], mask))
            rows.append(compare_structures(r, struct[r], ng, struct[ng], mask))

    # Within-target replicate comparisons: useful robustness checks.
    rows.append(compare_pairformer('RfaH_1', pair['RfaH_1'], 'RfaH_2', pair['RfaH_2'], mask))
    rows.append(compare_pairformer('NusG_1', pair['NusG_1'], 'NusG_2', pair['NusG_2'], mask))
    rows.append(compare_structures('RfaH_1', struct['RfaH_1'], 'RfaH_2', struct['RfaH_2'], mask))
    rows.append(compare_structures('NusG_1', struct['NusG_1'], 'NusG_2', struct['NusG_2'], mask))

    fieldnames = [
        'representation', 'A', 'B', 'pearson_R', 'jaccard', 'jaccard_type',
        'strict_topL_jaccard', 'A_contact_count', 'B_contact_count',
        'shared_or_credited_count', 'union_count'
    ]
    with OUT.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f'n={n}; eligible pairs={int(mask.sum())}; CTD={CTD_START}-{CTD_END}; L={TOP_L}')
    print('\nCROSS-TARGET COMPARISONS')
    for row in rows[:8]:
        strict = '' if row['strict_topL_jaccard'] == '' else f", strict top-L J={row['strict_topL_jaccard']:.6f}"
        print(f"{row['representation']:18s} {row['A']:7s} vs {row['B']:7s}: "
              f"R={row['pearson_R']:.6f}, J={row['jaccard']:.6f}{strict}; "
              f"counts {row['A_contact_count']}/{row['B_contact_count']}, "
              f"shared/credited={row['shared_or_credited_count']}, union={row['union_count']}")

    print('\nWITHIN-TARGET REPLICATE COMPARISONS')
    for row in rows[8:]:
        strict = '' if row['strict_topL_jaccard'] == '' else f", strict top-L J={row['strict_topL_jaccard']:.6f}"
        print(f"{row['representation']:18s} {row['A']:7s} vs {row['B']:7s}: "
              f"R={row['pearson_R']:.6f}, J={row['jaccard']:.6f}{strict}")

    print(f'\nWrote {OUT}')


if __name__ == '__main__':
    main()
