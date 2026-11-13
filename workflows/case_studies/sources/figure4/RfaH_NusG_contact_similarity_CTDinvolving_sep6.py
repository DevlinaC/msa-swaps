#!/usr/bin/env python3
"""RfaH/NusG CTD-involving contact-map similarity with |i-j| >= 6.

Compares:
  1) AF2 structural contact maps for RfaH and NusG.
  2) MSA-derived contact-probability maps for the same query-MSA inputs.

Pair universe
-------------
Upper-triangle residue pairs satisfying BOTH:
  * |i-j| >= 6
  * at least one residue is in the CTD, residues 110-162 (1-based inclusive)

Structural contacts
-------------------
C-beta--C-beta <= 8.0 A; C-alpha is used for Gly.
Pearson R is calculated between the two binary structural-contact vectors over
all eligible residue pairs. Jaccard is ordinary exact set Jaccard.

MSA-derived contacts
--------------------
Pearson R is calculated between the two continuous contact-probability vectors
over all eligible residue pairs.

Jaccard uses the manuscript's rank-tolerant definition. Let S and T be the
top-L contact sets and S_w and T_w the top-floor(1.5L) sets. A top-L pair
receives agreement credit if the identical pair is present in the other map's
wider set. Thus:

    credited = (S intersect T_w) union (T intersect S_w)
    tolerant Jaccard = |credited| / |S union T|

Score ties at rank cutoffs are retained. For this CTD-focused analysis,
L = CTD length = 53 residues.

Outputs
-------
  RfaH_NusG_contact_similarity_CTDinvolving_sep6_metrics.csv
  RfaH_NusG_contact_similarity_CTDinvolving_sep6_table.png
  RfaH_NusG_contact_similarity_CTDinvolving_sep6_table.svg
"""

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import pearsonr

BASE = Path(__file__).resolve().parent

PDB_RFAH = BASE / "RfaH_MSA_94a41_unrelaxed_rank_001_alphafold2_model_2_seed_002.pdb"
PDB_NUSG = BASE / "NusG_nogaps_unrelaxed_rank_001_alphafold2_model_3_seed_004.pdb"
PAIR_RFAH = BASE / "RfaH_unbiased_predicted_cb_contacts.txt"
PAIR_NUSG = BASE / "NusG_predicted_cb_contacts.txt"

OUT_CSV = BASE / "RfaH_NusG_contact_similarity_CTDinvolving_sep6_metrics.csv"
OUT_PNG = BASE / "RfaH_NusG_contact_similarity_CTDinvolving_sep6_table.png"
OUT_SVG = BASE / "RfaH_NusG_contact_similarity_CTDinvolving_sep6_table.svg"

CTD_START = 110
CTD_END = 162
CTD_LENGTH = CTD_END - CTD_START + 1  # 53
MIN_SEQUENCE_SEPARATION = 6
STRUCTURAL_CONTACT_CUTOFF_A = 8.0
TOP_L = CTD_LENGTH
TOP_L_MULTIPLIER = 1.5

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_pdb_contact_coordinates(path: Path):
    """Return sequence, residue keys, and C-beta coords (C-alpha for Gly)."""
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


def eligible_pair_mask(n: int):
    """Upper triangle, |i-j|>=6, and >=1 residue in CTD 110-162."""
    i, j = np.indices((n, n))
    ctd0 = CTD_START - 1
    ctd1 = CTD_END  # Python-exclusive upper bound equals 162
    ctd_involving = ((i >= ctd0) & (i < ctd1)) | ((j >= ctd0) & (j < ctd1))
    separated = np.abs(i - j) >= MIN_SEQUENCE_SEPARATION
    return np.triu(ctd_involving & separated, k=1)


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
    """Rank-tolerant top-L Jaccard used in the manuscript Methods."""
    wide_k = int(np.floor(TOP_L_MULTIPLIER * L))
    S, s_thr = ranked_set_with_ties(a, mask, L)
    T, t_thr = ranked_set_with_ties(b, mask, L)
    S_w, sw_thr = ranked_set_with_ties(a, mask, wide_k)
    T_w, tw_thr = ranked_set_with_ties(b, mask, wide_k)

    credited = (S & T_w) | (T & S_w)
    denom = S | T
    tolerant = len(credited) / len(denom) if denom else float("nan")
    strict = len(S & T) / len(denom) if denom else float("nan")

    return {
        "tolerant_jaccard": tolerant,
        "strict_topL_jaccard": strict,
        "credited_count": len(credited),
        "topL_union_count": len(denom),
        "topL_RfaH_count_with_ties": len(S),
        "topL_NusG_count_with_ties": len(T),
        "wide_RfaH_count_with_ties": len(S_w),
        "wide_NusG_count_with_ties": len(T_w),
        "topL_RfaH_threshold": s_thr,
        "topL_NusG_threshold": t_thr,
        "wide_RfaH_threshold": sw_thr,
        "wide_NusG_threshold": tw_thr,
    }


def make_table_figure(struct_r, struct_j, predicted_r, predicted_j):
    """Compact grayscale manuscript-style similarity table."""
    fig, ax = plt.subplots(figsize=(10.8, 4.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left = 0.03
    bottom = 0.17
    width = 0.94
    height = 0.72
    col_edges = [left, left + width * 0.50, left + width * 0.75, left + width]
    row_edges = [bottom, bottom + height / 3, bottom + 2 * height / 3, bottom + height]

    fills = {
        (0, 0): 1.0, (0, 1): 1.0, (0, 2): 1.0,
        (1, 0): 1.0, (1, 1): max(0.0, 1.0 - struct_r), (1, 2): max(0.0, 1.0 - struct_j),
        (2, 0): 1.0, (2, 1): max(0.0, 1.0 - predicted_r), (2, 2): max(0.0, 1.0 - predicted_j),
    }

    for r in range(3):
        y0, y1 = row_edges[2-r], row_edges[3-r]
        for c in range(3):
            x0, x1 = col_edges[c], col_edges[c+1]
            gray = fills[(r, c)]
            ax.add_patch(Rectangle(
                (x0, y0), x1-x0, y1-y0,
                facecolor=str(gray), edgecolor="black", linewidth=2.3,
            ))

    headers = ["Representation", "Pearson R", "Jaccard"]
    for c, txt in enumerate(headers):
        ax.text(
            (col_edges[c]+col_edges[c+1])/2,
            (row_edges[2]+row_edges[3])/2,
            txt, ha="center", va="center", fontsize=24, fontweight="bold",
        )

    rows = [
        ("Structural contacts", struct_r, struct_j),
        ("MSA contacts", predicted_r, predicted_j),
    ]
    y_centers = [(row_edges[1]+row_edges[2])/2, (row_edges[0]+row_edges[1])/2]
    for ridx, (label, rvalue, jvalue) in enumerate(rows):
        y = y_centers[ridx]
        ax.text(
            (col_edges[0]+col_edges[1])/2, y, label,
            ha="center", va="center", fontsize=24, fontweight="bold",
        )
        for c, value in [(1, rvalue), (2, jvalue)]:
            bg_gray = fills[(ridx+1, c)]
            text_color = "white" if bg_gray < 0.48 else "black"
            ax.text(
                (col_edges[c]+col_edges[c+1])/2, y, f"{value:.2f}",
                ha="center", va="center", fontsize=27, fontweight="bold",
                color=text_color,
            )

    #ax.text(
    #    0.5, 0.073,
    #    r"CTD-involving pairs (residues 110-162), $|i-j|\geq6$.",
    #    ha="center", va="center", fontsize=11.5,
    #)
    #ax.text(
    #    0.5, 0.027,
    #    r"Structural Jaccard: exact contact-set overlap; MSA-derived Jaccard: rank-tolerant top-$L$ ($L=53$).",
    #    ha="center", va="center", fontsize=11.2,
    #)

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_SVG, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    predicted_rfah = np.loadtxt(PAIR_RFAH)
    predicted_nusg = np.loadtxt(PAIR_NUSG)
    if (
        predicted_rfah.shape != predicted_nusg.shape
        or predicted_rfah.ndim != 2
        or predicted_rfah.shape[0] != predicted_rfah.shape[1]
    ):
        raise ValueError("The two predicted contact matrices must be square and have identical dimensions.")

    n = predicted_rfah.shape[0]
    seq_rfah, keys_rfah, coords_rfah = parse_pdb_contact_coordinates(PDB_RFAH)
    seq_nusg, keys_nusg, coords_nusg = parse_pdb_contact_coordinates(PDB_NUSG)

    if len(coords_rfah) != n or len(coords_nusg) != n:
        raise ValueError(
            f"Length mismatch: matrices are {n}x{n}; PDB lengths are "
            f"{len(coords_rfah)} and {len(coords_nusg)}"
        )
    if [k[1:] for k in keys_rfah] != [k[1:] for k in keys_nusg]:
        raise ValueError("PDB residue numbering/insertion-code indexing differs between RfaH and NusG.")
    if n < CTD_END:
        raise ValueError(f"Sequence length {n} is shorter than CTD_END={CTD_END}")

    mask = eligible_pair_mask(n)
    struct_rfah, _ = structural_contact_map(coords_rfah)
    struct_nusg, _ = structural_contact_map(coords_nusg)

    # Structural-map similarity over the CTD-involving pair universe.
    struct_r = float(
        pearsonr(
            struct_rfah[mask].astype(float),
            struct_nusg[mask].astype(float),
        ).statistic
    )
    struct_set_rfah = pair_set(struct_rfah, mask)
    struct_set_nusg = pair_set(struct_nusg, mask)
    struct_j = float(jaccard(struct_set_rfah, struct_set_nusg))

    # MSA-derived map similarity over the exact same pair universe.
    predicted_r = float(
        pearsonr(
            predicted_rfah[mask].astype(float),
            predicted_nusg[mask].astype(float),
        ).statistic
    )
    tolerant = rank_tolerant_jaccard(
        predicted_rfah, predicted_nusg, mask, TOP_L
    )
    predicted_j = float(tolerant["tolerant_jaccard"])

    sequence_identity = sum(a == b for a, b in zip(seq_rfah, seq_nusg)) / n

    metrics = {
        "sequence_RfaH": seq_rfah,
        "sequence_NusG": seq_nusg,
        "length": n,
        "sequence_identity_fraction": sequence_identity,
        "CTD_start_1based": CTD_START,
        "CTD_end_1based": CTD_END,
        "CTD_length": CTD_LENGTH,
        "minimum_sequence_separation_inclusive": MIN_SEQUENCE_SEPARATION,
        "eligible_upper_triangle_pairs": int(mask.sum()),
        "structural_contact_cutoff_A": STRUCTURAL_CONTACT_CUTOFF_A,
        "RfaH_structural_contact_count": len(struct_set_rfah),
        "NusG_structural_contact_count": len(struct_set_nusg),
        "shared_structural_contact_count": len(struct_set_rfah & struct_set_nusg),
        "structural_contact_union_count": len(struct_set_rfah | struct_set_nusg),
        "structural_Pearson_R": struct_r,
        "structural_Jaccard": struct_j,
        "MSA_derived_Pearson_R": predicted_r,
        "MSA_derived_top_L": TOP_L,
        "MSA_derived_rank_tolerant_Jaccard": predicted_j,
        **{f"MSA_derived_{k}": v for k, v in tolerant.items() if k != "tolerant_jaccard"},
    }

    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        for key, value in metrics.items():
            writer.writerow([key, value])

    make_table_figure(struct_r, struct_j, predicted_r, predicted_j)

    print("RfaH vs NusG: CTD-involving pairs, |i-j| >= 6")
    print(f"CTD: residues {CTD_START}-{CTD_END} (L={CTD_LENGTH})")
    print(f"Eligible upper-triangle pairs: {mask.sum()}")
    print(
        "Structural contacts: "
        f"RfaH={len(struct_set_rfah)}, NusG={len(struct_set_nusg)}, "
        f"shared={len(struct_set_rfah & struct_set_nusg)}, "
        f"union={len(struct_set_rfah | struct_set_nusg)}"
    )
    print(f"Structural Pearson R = {struct_r:.9f}")
    print(f"Structural Jaccard = {struct_j:.9f}")
    print(f"MSA-derived Pearson R = {predicted_r:.9f}")
    print(f"MSA-derived rank-tolerant Jaccard = {predicted_j:.9f}")
    print(f"MSA-derived strict top-L Jaccard (diagnostic) = {tolerant['strict_topL_jaccard']:.9f}")
    print(
        "Credited/union pairs for tolerant Jaccard = "
        f"{tolerant['credited_count']}/{tolerant['topL_union_count']}"
    )
    print(f"Top-L counts with ties: RfaH={tolerant['topL_RfaH_count_with_ties']}, NusG={tolerant['topL_NusG_count_with_ties']}")
    print(f"Wide-set counts with ties: RfaH={tolerant['wide_RfaH_count_with_ties']}, NusG={tolerant['wide_NusG_count_with_ties']}")
    print(f"Wrote {OUT_CSV}")
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_SVG}")


if __name__ == "__main__":
    main()
