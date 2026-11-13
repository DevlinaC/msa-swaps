# Import libraries.
import glob # Finds files matching a filename pattern.
import os # Works w/ file paths & checks if files exist.

import numpy as np # Used for fast array calc.
import pandas as pd # Used for reading & manipulating tables.


# --- --- ---
# GREMLIN pipeline. Unchanged from dual_jaccard_v2.py, EXCEPT:
#   - shared_alignment_aligned() now redirects to 2mm4_fixes/ for the 8
#     target/2mm4 pairs whose original .a3m/.fas were corrupted & never
#     existed under original_junk_msas/. All other pairs are unaffected.
#
# This code compares predicted residue contacts BT:
#   1. An aligned (swapped) protein (A_B_aligned).
#   2. Protein A's native sequence.
#   3. Protein B's native MSA.
# Simularity is measured using the Jaccard index.
#
# gap-column filtering (25% cutoff) removes columns
# w/ too many gaps, z-score floor (min_score=1.0),
# fixed top_n=100 (100 highest-scoring contacts),
# seq_sep=5 (ignores residues closer than 5 positions),
# offset search -5 to 5 when comparing to protein B.
# --- --- ---

# Targets whose original {target}_2mm4_aligned.a3m/.fas were corrupted &
# corrected; corrected files live in 2mm4_fixes/{target}_2mm4_aligned/
# instead of original_junk_msas/{target}/{target}_2mm4_aligned/.
FIXED_2MM4_TARGETS = {"1bor", "1bw6", "1cqu", "1gb1", "1vex", "2fs1", "2nc2", "2ovo"}
FIXES_DIR = "/data/AF_interpretation/msa_swaps/2mm4_fixes/"


def parse_gremlin(file_path):
    """
    Reads a GREMLIN output CSV. Keeps only the columns needed
    for the comparison.

    Returns:
        DataFrame w/ i (1st residue position), j (2nd residue position),
        score (GREMLIN z-score).
    """
    df = pd.read_csv(file_path)
    df.columns = [c.strip() for c in df.columns] # Remove accidental spaces from column names.
    out = df[["i", "j", "zscore"]].copy()
    out = out.rename(columns={"zscore": "score"})
    return out


def read_fasta_sequences(fasta_path):
    """
    Read every sequence from a FASTA alignment file.

    Returns:
        A list of protein sequences as strings.
    """
    sequences = []
    current = []
    with open(fasta_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line: # Skip empty lines.
                continue
            if line.startswith(">"): # Line starting w/ ">" marks a new sequence.
                if current:
                    sequences.append("".join(current))
                    current = []
            else:
                current.append(line.upper()) # Add sequence letters.
        if current:
            sequences.append("".join(current)) # Save the final sequence.
    return sequences


def compute_gap_mask(fasta_path, gap_cutoff=0.25):
    """
    Determine which alignment columns have acceptable #s of gaps. Columns where >
    gap_cutoff (25%) of sequences have '-' are removed.

    Returns:
        Boolean array.
        True = keep this column.
        False = ignore this column.
    """
    try:
        sequences = read_fasta_sequences(fasta_path)
    except FileNotFoundError:
        print(f"Warning: alignment file not found at {fasta_path}, skipping gap filtering.", flush=True)
        return None
    if not sequences:
        return None
    lengths = {len(s) for s in sequences}
    if len(lengths) != 1:
        print(f"Warning: sequences in {fasta_path} aren't aligned (uneven lengths {sorted(lengths)}), "
              f"skipping gap filtering.", flush=True)
        return None
    arr = np.array([list(s) for s in sequences])
    gap_fraction = np.mean(arr == "-", axis=0)
    return gap_fraction < gap_cutoff


def apply_gap_mask(df, gap_mask, file_path=""):
    """
    Remove predicted contacts involving residues located in high-gap
    alignment columns.
    """
    if gap_mask is None or df.empty:
        return df
    max_index = int(max(df["i"].max(), df["j"].max()))
    if max_index >= len(gap_mask):
        print(f"Warning: index {max_index} exceeds gap_mask length {len(gap_mask)} "
              f"for {file_path}; skipping gap filter for this file.", flush=True)
        return df
    keep = gap_mask[df["i"].to_numpy()] & gap_mask[df["j"].to_numpy()]
    return df[keep]


def get_gremlin_contact_set(file_path, seq_sep=5, top_n=100, offset=0, min_score=1.0, gap_mask=None):
    """
    Convert a GREMLIN prediction file into a set of residue pairs.
    """
    if not os.path.isfile(file_path):
        print(f"Warning: File not found at {file_path}", flush=True)
        return set()
    df = parse_gremlin(file_path) # Read predictions.
    if gap_mask is not None:
        df = apply_gap_mask(df, gap_mask, file_path=file_path) # Remove high-gap postions.
    if min_score is not None:
        df = df[df["score"] >= min_score] # Keep only scores above the threshold.
    df_filtered = df[abs(df["i"] - df["j"]) >= seq_sep].copy() # Remove contacts BT nearby residues.
    df_sorted = df_filtered.sort_values(by="score", ascending=False)
    df_top_n = df_sorted.head(top_n) # Keep the top N predictions.
    contact_pairs = set()
    for _, row in df_top_n.iterrows(): # Apply the optional residue offset.
        i_offset = int(row["i"]) + offset
        j_offset = int(row["j"]) + offset
        contact_pairs.add(tuple(sorted([i_offset, j_offset])))
    return contact_pairs # Return a set of residue pairs.


def calculate_jaccard(set_a, set_b):
    """
    Compute the Jaccard simularity BT 2 sets.

    Jaccard = (# shared contacts) / (# unique contacts)

    Values range from:
        0 = no overlap.
        1 = identical sets.
    """
    if not set_a or not set_b:
        return 0.0
    intersection = len(set_a.intersection(set_b))
    union = len(set_a.union(set_b))
    return intersection / union


def parse_aligned_filename(filename):
    stem = filename.replace(".csv", "")
    parts = stem.split("_")
    return parts[1], parts[2]


def process_one_gremlin_file(file_aligned, native_path_fn_A, native_path_fn_B,
                              alignment_path_fn_aligned, alignment_path_fn_A, alignment_path_fn_B,
                              seq_sep, top_n, min_score, gap_cutoff, offset_range):
    """
    Process 1 aligned GREMLIN prediction file.

    Returns:
        A dict w/ the Jaccard scores.
    """
    filename = os.path.basename(file_aligned)
    protein_A, protein_B = parse_aligned_filename(filename)
    file_native_A = native_path_fn_A(protein_A, protein_B)
    file_native_B = native_path_fn_B(protein_A, protein_B)

    gap_mask_aligned = gap_mask_A = gap_mask_B = None
    if gap_cutoff is not None:
        gap_mask_aligned = compute_gap_mask(alignment_path_fn_aligned(protein_A, protein_B), gap_cutoff)
        gap_mask_A = compute_gap_mask(alignment_path_fn_A(protein_A, protein_B), gap_cutoff)
        gap_mask_B = compute_gap_mask(alignment_path_fn_B(protein_A, protein_B), gap_cutoff)

    set_aligned_vs_A = get_gremlin_contact_set(file_aligned, seq_sep, top_n, 0, min_score, gap_mask_aligned)
    set_native_A = get_gremlin_contact_set(file_native_A, seq_sep, top_n, 0, min_score, gap_mask_A)
    jaccard_vs_A = calculate_jaccard(set_aligned_vs_A, set_native_A)

    set_native_B = get_gremlin_contact_set(file_native_B, seq_sep, top_n, 0, min_score, gap_mask_B)
    best_score, best_offset = 0.0, 0
    for offset in offset_range:
        set_aligned = get_gremlin_contact_set(file_aligned, seq_sep, top_n, offset, min_score, gap_mask_aligned)
        score = calculate_jaccard(set_aligned, set_native_B)
        if score > best_score:
            best_score, best_offset = score, offset

    return {
        "Aligned_File": filename, "Protein_A": protein_A, "Protein_B": protein_B,
        "Jaccard_vs_ProteinA": jaccard_vs_A, "Jaccard_vs_ProteinB": best_score,
        "Best_Offset_ProteinB": best_offset,
    }


def run_gremlin_pipeline(aligned_files, native_path_fn_A, native_path_fn_B,
                          alignment_path_fn_aligned, alignment_path_fn_A, alignment_path_fn_B,
                          output_csv, seq_sep=5, top_n=100, min_score=1.0,
                          gap_cutoff=0.25, offset_range=range(-5, 6), label="GREMLIN"):
    """
    Runs the complete GREMLIN analysis on every aligned file.

    Results are written to a CSV file, & a summary of the average Jaccard scores is
    returned.
    """
    results, failures = [], []
    for f in aligned_files:
        try:
            results.append(process_one_gremlin_file(
                f, native_path_fn_A, native_path_fn_B,
                alignment_path_fn_aligned, alignment_path_fn_A, alignment_path_fn_B,
                seq_sep, top_n, min_score, gap_cutoff, offset_range))
        except Exception as exc:
            failures.append((os.path.basename(f), str(exc)))

    output_df = pd.DataFrame(results)
    output_df.to_csv(output_csv, index=False)
    mean_A = output_df["Jaccard_vs_ProteinA"].mean() if not output_df.empty else float("nan")
    mean_B = output_df["Jaccard_vs_ProteinB"].mean() if not output_df.empty else float("nan")
    median_B = output_df["Jaccard_vs_ProteinB"].median() if not output_df.empty else float("nan")

    summary = (f"[{label}] {len(results)}/{len(aligned_files)} rows ({len(failures)} failed) | "
               f"seq_sep={seq_sep}, top_n={top_n}, min_score={min_score}, gap_cutoff={gap_cutoff} "
               f"(strict Jaccard, unchanged method) | "
               f"A_mean={mean_A:.4f} | B_mean={mean_B:.4f}, B_median={median_B:.4f} | "
               f"separation(B-A)={mean_B - mean_A:.4f}")
    if failures:
        summary += f"\n[{label}] First few failures: {failures[:5]}"
    return summary

# --- --- ---
# Dataset-specific configuration
# --- --- ---

GREMLIN_SEQ_SEP = 5 # Ignore contacts BT residues closer than this distance.
GREMLIN_TOP_N = 100 # Compare only the 100 strongest GREMLIN predictions.
GREMLIN_MIN_SCORE = 1.0 # Min GREMLIN confidence score (z-score).
GREMLIN_GAP_CUTOFF = 0.25 # Remove alignement coolumns w/ > 25% gaps.
GREMLIN_OFFSET_RANGE = range(-5, 6) # Try residue shifts from -5 to +5 when matching Protein B (MSA protein).

GREMLIN_JUNK_DIR = "/data/AF_interpretation/msa_swaps/gremlin/results_gremlin_junk/"
GREMLIN_PROPER_DIR = "/data/AF_interpretation/msa_swaps/gremlin/results_gremlin_proper/"
ALIGNMENT_JUNK_DIR = "/data/AF_interpretation/msa_swaps/original_junk_msas/"
ALIGNMENT_PROPER_DIR = "/data/AF_interpretation/msa_swaps/original_proper_msas/"


def gremlin_path_A(protein_A, protein_B):
    return os.path.join(GREMLIN_PROPER_DIR, f"full_{protein_A}.csv")


def gremlin_path_B(protein_A, protein_B):
    return os.path.join(GREMLIN_PROPER_DIR, f"full_{protein_B}.csv")


def shared_alignment_aligned(protein_A, protein_B):
    # The 8 target/2mm4 pairs whose original .a3m/.fas were corrupted never
    # got a folder under original_junk_msas/; their corrected .fas lives in
    # 2mm4_fixes/ instead. Every other pair uses the normal location.
    if protein_B == "2mm4" and protein_A in FIXED_2MM4_TARGETS:
        pair_name = f"{protein_A}_{protein_B}_aligned"
        return os.path.join(FIXES_DIR, pair_name, f"{pair_name}.fas")
    return os.path.join(
        ALIGNMENT_JUNK_DIR,
        protein_A,
        f"{protein_A}_{protein_B}_aligned",
        f"{protein_A}_{protein_B}_aligned.fas",
    )


def shared_alignment_A(protein_A, protein_B):
    return os.path.join(ALIGNMENT_PROPER_DIR, protein_A, f"{protein_A}.fas")


def shared_alignment_B(protein_A, protein_B):
    return os.path.join(ALIGNMENT_PROPER_DIR, protein_B, f"{protein_B}.fas")


def run_gremlin_job():
    """
    Find every aligned GREMLIN result & run the full pipeline.
    """
    aligned_files = sorted(
        glob.glob(os.path.join(GREMLIN_JUNK_DIR, "full_*_aligned.csv"))
    )
    return run_gremlin_pipeline(
        aligned_files,
        gremlin_path_A,
        gremlin_path_B,
        shared_alignment_aligned,
        shared_alignment_A,
        shared_alignment_B,
        output_csv="jaccard_results_gremlin.csv",
        seq_sep=GREMLIN_SEQ_SEP,
        top_n=GREMLIN_TOP_N,
        min_score=GREMLIN_MIN_SCORE,
        gap_cutoff=GREMLIN_GAP_CUTOFF,
        offset_range=GREMLIN_OFFSET_RANGE,
    )


if __name__ == "__main__":
    print(run_gremlin_job(), flush=True)
