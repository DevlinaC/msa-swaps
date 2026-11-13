"""Analysis functions preserved from gremlin_jaccard_v2.py.

Dataset discovery and fixed cluster paths are replaced by compute_contacts.py.
See docs/file_sources.csv for original source checksums.
"""
import os
import numpy as np
import pandas as pd


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
    return contact_pairs


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
