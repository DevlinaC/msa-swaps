#!/usr/bin/env python3
"""GB98 T25V/T25E contact-map similarity with |i-j| >= 6.

This script compares:
  1) AF2 structural contact maps for GB98-T25V and GB98-T25E.
  2) MSA-derived contact-probability maps for the same two query-MSA inputs.

Definitions
-----------
Eligible residue pairs:
    Upper-triangle pairs with |i-j| >= 6.

Structural contacts:
    Cbeta-Cbeta <= 8.0 A; Calpha is used for Gly.
    Pearson R is calculated between the two binary structural-contact vectors
    over all eligible residue pairs.
    Jaccard is the ordinary set Jaccard of structural contacts.

MSA-derived contacts:
    Pearson R is calculated between the two contact-probability vectors over
    all eligible residue pairs.
    Jaccard uses the manuscript's rank-tolerant definition: let S and T be the
    top-L contact sets and S_w and T_w be the top-floor(1.5L) sets. A top-L
    pair receives agreement credit if the identical pair is present in the
    other map's top-floor(1.5L) set. Thus the credited set is
        (S intersect T_w) union (T intersect S_w)
    and tolerant Jaccard = |credited set| / |S union T|.
    Score ties at the rank cutoffs are retained.

Outputs
-------
  GB98_T25V_T25E_contact_similarity_sep6_metrics.csv
  GB98_T25V_T25E_contact_similarity_sep6_table.png
  GB98_T25V_T25E_contact_similarity_sep6_table.svg
"""

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import pearsonr

BASE = Path(__file__).resolve().parent

PDB_T25V = BASE / "T25V_unrelaxed_rank_001_alphafold2_model_2_seed_003.pdb"
PDB_T25E = BASE / "T25E_unrelaxed_rank_001_alphafold2_model_3_seed_003.pdb"
PAIR_T25V = BASE / "T52V.txt"
PAIR_T25E = BASE / "predicted_T25E_contacts.txt"

OUT_CSV = BASE / "GB98_T25V_T25E_contact_similarity_sep6_metrics.csv"
OUT_PNG = BASE / "GB98_T25V_T25E_contact_similarity_sep6_table.png"
OUT_SVG = BASE / "GB98_T25V_T25E_contact_similarity_sep6_table.svg"

MIN_SEQUENCE_SEPARATION = 6
STRUCTURAL_CONTACT_CUTOFF_A = 8.0
TOP_L_MULTIPLIER = 1.5

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_pdb_contact_coordinates(path: Path):
    """Return sequence, residue keys, and Cbeta coords (Calpha for Gly)."""
    residues = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            chain = line[21]
            resnum = int(line[22:26])
            icode = line[26]
            resname = line[17:20].strip()
            atom = line[12:16].strip()
            key = (chain, resnum, icode)
            residues.setdefault(key, {"resname": resname, "atoms": {}})
            residues[key]["atoms"][atom] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )

    sequence, keys, coords = [], [], []
    for key, entry in residues.items():
        resname = entry["resname"]
        atom_name = "CA" if resname == "GLY" else "CB"
        if atom_name not in entry["atoms"]:
            raise ValueError(f"Missing {atom_name} for residue {key} ({resname}) in {path}")
        sequence.append(THREE_TO_ONE.get(resname, "X"))
        keys.append(key)
        coords.append(entry["atoms"][atom_name])

    return "".join(sequence), keys, np.vstack(coords)


def eligible_pair_mask(n: int, min_sep: int = MIN_SEQUENCE_SEPARATION):
    """Upper-triangle mask for |i-j| >= min_sep."""
    i, j = np.indices((n, n))
    return np.triu(np.ones((n, n), dtype=bool), k=1) & (np.abs(i - j) >= min_sep)


def structural_contact_map(coords: np.ndarray, cutoff: float = STRUCTURAL_CONTACT_CUTOFF_A):
    delta = coords[:, None, :] - coords[None, :, :]
    distances = np.sqrt(np.sum(delta * delta, axis=2))
    return distances <= cutoff, distances


def pair_set(binary_map: np.ndarray, mask: np.ndarray):
    return {tuple(x) for x in np.argwhere(binary_map & mask)}


def jaccard(a: set, b: set):
    union = a | b
    return len(a & b) / len(union) if union else float("nan")


def ranked_set_with_ties(scores: np.ndarray, mask: np.ndarray, k: int):
    """Top-k eligible pairs, retaining all score ties at the kth rank."""
    values = scores[mask]
    if k < 1 or k > values.size:
        raise ValueError(f"Invalid k={k}; eligible pair count={values.size}")
    threshold = np.sort(values)[::-1][k - 1]
    selected = mask & (scores >= threshold)
    return {tuple(x) for x in np.argwhere(selected)}, float(threshold)


def rank_tolerant_jaccard(a: np.ndarray, b: np.ndarray, mask: np.ndarray, L: int):
    """Rank-tolerant top-L Jaccard exactly as specified in the Methods."""
    wide_k = int(np.floor(TOP_L_MULTIPLIER * L))
    S, s_thr = ranked_set_with_ties(a, mask, L)
    T, t_thr = ranked_set_with_ties(b, mask, L)
    S_w, sw_thr = ranked_set_with_ties(a, mask, wide_k)
    T_w, tw_thr = ranked_set_with_ties(b, mask, wide_k)

    credited = (S & T_w) | (T & S_w)
    denom = S | T
    score = len(credited) / len(denom) if denom else float("nan")

    # Strict top-L Jaccard is retained as a useful diagnostic.
    strict = len(S & T) / len(denom) if denom else float("nan")
    return {
        "tolerant_jaccard": score,
        "strict_topL_jaccard": strict,
        "credited_count": len(credited),
        "topL_union_count": len(denom),
        "topL_A_count_with_ties": len(S),
        "topL_B_count_with_ties": len(T),
        "wide_A_count_with_ties": len(S_w),
        "wide_B_count_with_ties": len(T_w),
        "topL_A_threshold": s_thr,
        "topL_B_threshold": t_thr,
        "wide_A_threshold": sw_thr,
        "wide_B_threshold": tw_thr,
    }


def make_table_figure(struct_r, struct_j, predicted_r, predicted_j):
    """Compact grayscale table matching the style used in the manuscript."""
    fig, ax = plt.subplots(figsize=(10.8, 4.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Table geometry.
    left = 0.03
    bottom = 0.17
    width = 0.94
    height = 0.72
    col_edges = [left, left + width * 0.50, left + width * 0.75, left + width]
    row_edges = [bottom, bottom + height / 3, bottom + 2 * height / 3, bottom + height]

    # Cell fills. Similarity cells use grayscale intensity to preserve print compatibility.
    fills = {
        (0, 0): 1.0, (0, 1): 1.0, (0, 2): 1.0,
        (1, 0): 1.0, (1, 1): max(0.0, 1.0 - struct_r), (1, 2): max(0.0, 1.0 - struct_j),
        (2, 0): 1.0, (2, 1): max(0.0, 1.0 - predicted_r), (2, 2): max(0.0, 1.0 - predicted_j),
    }

    # row index here: 0 header, 1 structural, 2 MSA-derived; map to y from top.
    for r in range(3):
        y0, y1 = row_edges[2-r], row_edges[3-r]
        for c in range(3):
            x0, x1 = col_edges[c], col_edges[c+1]
            gray = fills[(r, c)]
            ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0,
                                   facecolor=str(gray), edgecolor="black", linewidth=2.3))

    headers = ["Representation", "Pearson R", "Jaccard"]
    for c, txt in enumerate(headers):
        ax.text((col_edges[c]+col_edges[c+1])/2, (row_edges[2]+row_edges[3])/2,
                txt, ha="center", va="center", fontsize=24, fontweight="bold")

    rows = [
        ("Structural contacts", struct_r, struct_j),
        ("MSA contacts", predicted_r, predicted_j),
    ]
    y_centers = [(row_edges[1]+row_edges[2])/2, (row_edges[0]+row_edges[1])/2]
    for ridx, (label, rvalue, jvalue) in enumerate(rows):
        y = y_centers[ridx]
        ax.text((col_edges[0]+col_edges[1])/2, y, label,
                ha="center", va="center", fontsize=24, fontweight="bold")
        for c, value in [(1, rvalue), (2, jvalue)]:
            bg_gray = fills[(ridx+1, c)]
            text_color = "white" if bg_gray < 0.48 else "black"
            ax.text((col_edges[c]+col_edges[c+1])/2, y, f"{value:.2f}",
                    ha="center", va="center", fontsize=27, fontweight="bold", color=text_color)

    #ax.text(0.5, 0.045,
    #        r"Structural Jaccard: exact contact-set overlap; MSA-derived Jaccard: rank-tolerant top-$L$ ($L=56$).  $|i-j|\geq6$.",
    #        ha="center", va="center", fontsize=11.5)

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_SVG, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    predicted_v = np.loadtxt(PAIR_T25V)
    predicted_e = np.loadtxt(PAIR_T25E)
    if predicted_v.shape != predicted_e.shape or predicted_v.ndim != 2 or predicted_v.shape[0] != predicted_v.shape[1]:
        raise ValueError("The two predicted contact matrices must be square and have identical dimensions.")

    n = predicted_v.shape[0]
    seq_v, keys_v, coords_v = parse_pdb_contact_coordinates(PDB_T25V)
    seq_e, keys_e, coords_e = parse_pdb_contact_coordinates(PDB_T25E)
    if len(coords_v) != n or len(coords_e) != n:
        raise ValueError(f"Length mismatch: matrices are {n}x{n}; PDB lengths are {len(coords_v)} and {len(coords_e)}")
    if keys_v != keys_e:
        raise ValueError("PDB residue indexing differs between T25V and T25E.")

    mask = eligible_pair_mask(n)
    struct_v, _ = structural_contact_map(coords_v)
    struct_e, _ = structural_contact_map(coords_e)

    # Structural-map similarity.
    struct_r = float(pearsonr(struct_v[mask].astype(float), struct_e[mask].astype(float)).statistic)
    struct_set_v = pair_set(struct_v, mask)
    struct_set_e = pair_set(struct_e, mask)
    struct_j = float(jaccard(struct_set_v, struct_set_e))

    # MSA-derived map similarity.
    predicted_r = float(pearsonr(predicted_v[mask].astype(float), predicted_e[mask].astype(float)).statistic)
    tolerant = rank_tolerant_jaccard(predicted_v, predicted_e, mask, n)
    predicted_j = float(tolerant["tolerant_jaccard"])

    sequence_differences = [
        f"{i+1}:{a}>{b}" for i, (a, b) in enumerate(zip(seq_v, seq_e)) if a != b
    ]

    metrics = {
        "sequence_T25V": seq_v,
        "sequence_T25E": seq_e,
        "sequence_differences": ";".join(sequence_differences),
        "length_L": n,
        "minimum_sequence_separation": MIN_SEQUENCE_SEPARATION,
        "eligible_upper_triangle_pairs": int(mask.sum()),
        "structural_contact_cutoff_A": STRUCTURAL_CONTACT_CUTOFF_A,
        "T25V_structural_contact_count": len(struct_set_v),
        "T25E_structural_contact_count": len(struct_set_e),
        "shared_structural_contact_count": len(struct_set_v & struct_set_e),
        "structural_Pearson_R": struct_r,
        "structural_Jaccard": struct_j,
        "MSA_derived_Pearson_R": predicted_r,
        "MSA_derived_rank_tolerant_Jaccard": predicted_j,
        **{f"MSA_derived_{k}": v for k, v in tolerant.items() if k != "tolerant_jaccard"},
    }

    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        for k, v in metrics.items():
            writer.writerow([k, v])

    make_table_figure(struct_r, struct_j, predicted_r, predicted_j)

    print("GB98 T25V vs T25E, |i-j| >= 6")
    print(f"Eligible upper-triangle pairs: {mask.sum()}")
    print(f"Structural contacts: T25V={len(struct_set_v)}, T25E={len(struct_set_e)}, shared={len(struct_set_v & struct_set_e)}")
    print(f"Structural Pearson R = {struct_r:.9f}")
    print(f"Structural Jaccard = {struct_j:.9f}")
    print(f"MSA-derived Pearson R = {predicted_r:.9f}")
    print(f"MSA-derived rank-tolerant Jaccard = {predicted_j:.9f}")
    print(f"MSA-derived strict top-L Jaccard (diagnostic) = {tolerant['strict_topL_jaccard']:.9f}")
    print(f"Credited/union pairs for tolerant Jaccard = {tolerant['credited_count']}/{tolerant['topL_union_count']}")
    print(f"Wrote {OUT_CSV}")
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
