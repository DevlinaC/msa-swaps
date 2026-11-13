#!/usr/bin/env python3
"""
compute_msa_transformer_jaccard.py

Computes a Jaccard index between MSA Transformer coupling
matrices for swapped (junk) MSAs vs. proper MSAs, using the corrected/
deeper af_sampling data. For each target/partner pair (folder named
"{target}_{partner}_aligned", where target/Protein A is the sequence
source and partner/Protein B is the swapped-in MSA source), computes:

  Jaccard_vs_ProteinA: swapped matrix vs. target's OWN proper matrix
    Swapped: original_junk_msas_af_sampling/{target}/{target}_{partner}_aligned/
                 full_{target}_{partner}_aligned.csv
             (OR, for the 8 corrected 2mm4 pairs, 2mm4_fixes/{pair}/full_{pair}.csv)
    Proper:  original_proper_msas_af_sampling/{target}/full_{target}.csv

  Jaccard_vs_ProteinB: swapped matrix vs. partner's proper matrix
    Swapped: (same file as above)
    Proper:  original_proper_msas_af_sampling/{partner}/full_{partner}.csv

Both comparisons use the exact same tolerant-Jaccard method (center-crop
to shared size, top-L vs top-1.5L agreement, upper triangle, |i-j| >= 6) --
no offset search, matching the original single-comparison script.

The 8 target/2mm4 pairs whose original swapped .csv was corrupted never
got a full_*.csv under original_junk_msas_af_sampling/; their corrected
version lives in 2mm4_fixes/{pair}/full_{pair}.csv instead. This script
redirects to that location automatically for those 8 pairs only -- every
other pair uses the normal path, unchanged.

Writes one row per pair to a CSV with columns matching the GREMLIN Jaccard
output style: Aligned_File, Protein_A, Protein_B, Jaccard_vs_ProteinA,
Jaccard_vs_ProteinB.

Usage:
    python3 compute_msa_transformer_jaccard.py \\
        --base /data/AF_interpretation/msa_swaps \\
        --out-csv /data/schwartzlc/jaccard_results_msatransformer.csv
"""

import argparse
import csv
import glob
import os
import sys

import numpy as np

FIXED_2MM4_TARGETS = {"1bor", "1bw6", "1cqu", "1gb1", "1vex", "2fs1", "2nc2", "2ovo"}
FIXES_DIR = "/data/AF_interpretation/msa_swaps/2mm4_fixes"


def crop_center(m: np.ndarray, n: int) -> np.ndarray:
    """Crop symmetric matrix m to n x n, removing extras from both termini.
    Odd differences take one extra residue off the front (N-terminus)."""
    d = m.shape[0] - n
    if d == 0:
        return m
    front = (d + 1) // 2  # ceil -> 1 extra off the front when d is odd
    back = d // 2         # floor
    return m[front:m.shape[0] - back, front:m.shape[0] - back]


def top(m: np.ndarray, mask: np.ndarray, k: int) -> np.ndarray:
    thr = np.sort(m[mask])[-k]
    return (m >= thr) & mask


def tolerant_jaccard(matrix_1: np.ndarray, matrix_2: np.ndarray) -> float:
    """Crop both matrices to shared size and compute the tolerant Jaccard
    index (top-L vs top-1.5L agreement, upper triangle, |i-j| >= 6)."""
    n = min(matrix_1.shape[0], matrix_2.shape[0])
    matrix_1 = crop_center(matrix_1, n)
    matrix_2 = crop_center(matrix_2, n)
    assert matrix_1.shape == matrix_2.shape == (n, n)

    mask = np.triu(np.ones((n, n), bool), k=6)
    L = n

    A, B = top(matrix_1, mask, L), top(matrix_2, mask, L)
    Awide = top(matrix_1, mask, int(1.5 * L))
    Bwide = top(matrix_2, mask, int(1.5 * L))

    # agreed = top-L in one and top-1.5L in the other
    agreed = (A & Bwide) | (B & Awide)
    return agreed.sum() / (A | B).sum()


def get_swapped_csv_path(junk_dir: str, target: str, swap_path: str, swap: str, partner: str) -> str:
    """Return the path to the swapped matrix CSV, redirecting to
    2mm4_fixes/ for the 8 corrected target/2mm4 pairs."""
    if partner == "2mm4" and target in FIXED_2MM4_TARGETS:
        return os.path.join(FIXES_DIR, swap, f"full_{swap}.csv")
    return os.path.join(swap_path, f"full_{swap}.csv")


def main():
    ap = argparse.ArgumentParser(
        description="Compute tolerant Jaccard index (vs Protein A and vs Protein B) "
                    "between swapped and proper MSA Transformer coupling matrices "
                    "(af_sampling data, with 2mm4-fix redirect)."
    )
    ap.add_argument("--base", required=True,
                     help="Path to msa_swaps directory (contains "
                          "original_junk_msas_af_sampling / "
                          "original_proper_msas_af_sampling)")
    ap.add_argument("--out-csv", required=True,
                     help="Path to write per-pair results CSV")
    args = ap.parse_args()

    junk_dir = os.path.join(args.base, "original_junk_msas_af_sampling")
    proper_dir = os.path.join(args.base, "original_proper_msas_af_sampling")

    if not os.path.isdir(junk_dir):
        sys.exit(f"[ERROR] junk dir not found: {junk_dir}")
    if not os.path.isdir(proper_dir):
        sys.exit(f"[ERROR] proper dir not found: {proper_dir}")

    rows = []
    n_skipped = 0
    count = 0

    # --- normal pairs, discovered from junk_dir as before ---
    pdbs = glob.glob(os.path.join(junk_dir, "*/"))
    seen_pairs = set()

    for pdb_path in sorted(pdbs):
        target = pdb_path.rstrip("/").split("/")[-1]  # Protein A

        swaps = glob.glob(os.path.join(junk_dir, target, "*/"))
        for swap_path in sorted(swaps):
            swap = swap_path.rstrip("/").split("/")[-1]  # e.g. "1bor_1bw6_aligned"
            parts = swap.split("_")
            if len(parts) < 2:
                print(f"[WARN] unexpected swap folder name: {swap}, skipping", file=sys.stderr)
                n_skipped += 1
                continue
            partner = parts[1]  # Protein B (source of the swapped-in MSA)
            seen_pairs.add((target, partner))

            swapped_csv = get_swapped_csv_path(junk_dir, target, swap_path, swap, partner)
            proper_a_csv = os.path.join(proper_dir, target, f"full_{target}.csv")
            proper_b_csv = os.path.join(proper_dir, partner, f"full_{partner}.csv")

            missing = False
            for path in (swapped_csv, proper_a_csv, proper_b_csv):
                if not os.path.isfile(path):
                    print(f"[WARN] missing {path}, skipping {swap}", file=sys.stderr)
                    missing = True
            if missing:
                n_skipped += 1
                continue

            try:
                swapped_msa = np.loadtxt(swapped_csv, delimiter=",")
                proper_a_msa = np.loadtxt(proper_a_csv, delimiter=",")
                proper_b_msa = np.loadtxt(proper_b_csv, delimiter=",")

                jaccard_vs_a = tolerant_jaccard(swapped_msa, proper_a_msa)
                jaccard_vs_b = tolerant_jaccard(swapped_msa, proper_b_msa)
            except Exception as e:
                print(f"[WARN] failed on {swap} (A={target}, B={partner}): {e}", file=sys.stderr)
                n_skipped += 1
                continue

            rows.append([swap, target, partner, jaccard_vs_a, jaccard_vs_b])
            count += 1
            print(f"{swap}: vs_A={jaccard_vs_a:.5f}, vs_B={jaccard_vs_b:.5f}")

    # --- the 8 corrected 2mm4 pairs, not discoverable via junk_dir glob
    # (no folder exists there), so handle them explicitly here. ---
    for target in sorted(FIXED_2MM4_TARGETS):
        partner = "2mm4"
        if (target, partner) in seen_pairs:
            continue  # Already handled above (shouldn't happen but just in case).

        swap = f"{target}_{partner}_aligned"
        swapped_csv = os.path.join(FIXES_DIR, swap, f"full_{swap}.csv")
        proper_a_csv = os.path.join(proper_dir, target, f"full_{target}.csv")
        proper_b_csv = os.path.join(proper_dir, partner, f"full_{partner}.csv")

        missing = False
        for path in (swapped_csv, proper_a_csv, proper_b_csv):
            if not os.path.isfile(path):
                print(f"[WARN] missing {path}, skipping {swap}", file=sys.stderr)
                missing = True
        if missing:
            n_skipped += 1
            continue

        try:
            swapped_msa = np.loadtxt(swapped_csv, delimiter=",")
            proper_a_msa = np.loadtxt(proper_a_csv, delimiter=",")
            proper_b_msa = np.loadtxt(proper_b_csv, delimiter=",")

            jaccard_vs_a = tolerant_jaccard(swapped_msa, proper_a_msa)
            jaccard_vs_b = tolerant_jaccard(swapped_msa, proper_b_msa)
        except Exception as e:
            print(f"[WARN] failed on {swap} (A={target}, B={partner}): {e}", file=sys.stderr)
            n_skipped += 1
            continue

        rows.append([swap, target, partner, jaccard_vs_a, jaccard_vs_b])
        count += 1
        print(f"{swap}: vs_A={jaccard_vs_a:.5f}, vs_B={jaccard_vs_b:.5f}  [from 2mm4_fixes]")

    if count == 0:
        sys.exit("[ERROR] No pairs successfully processed; check paths/data.")

    with open(args.out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Aligned_File", "Protein_A", "Protein_B",
                          "Jaccard_vs_ProteinA", "Jaccard_vs_ProteinB"])
        writer.writerows(rows)

    a_scores = [r[3] for r in rows]
    b_scores = [r[4] for r in rows]

    print(f"\nWrote {count} pair result(s) to {args.out_csv}")
    if n_skipped:
        print(f"Skipped {n_skipped} pair(s) due to missing files or errors.", file=sys.stderr)
    print(f"vs Protein A -- Min: {min(a_scores):.5f}, Max: {max(a_scores):.5f}, "
          f"Avg: {sum(a_scores)/count:.5f}")
    print(f"vs Protein B -- Min: {min(b_scores):.5f}, Max: {max(b_scores):.5f}, "
          f"Avg: {sum(b_scores)/count:.5f}")


if __name__ == "__main__":
    main()
