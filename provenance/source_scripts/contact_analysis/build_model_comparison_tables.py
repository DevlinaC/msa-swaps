#!/usr/bin/env python3
"""
build_model_comparison_tables.py

Builds AF2-vs-AF3 comparison tables for the single-sequence (no-MSA) experiment.

For each experiment, writes 3 CSVs:
  - {experiment}_all.csv           every prediction, seed-matched
  - {experiment}_above_0.5.csv     rows where AF3's TM-score > 0.5
  - {experiment}_around_0.5.csv    rows where AF3's TM-score is in [0.4, 0.6]

PAIRING RULE ("seed-to-seed"): AF2 labels each prediction as
model_M_seed_S (5 models x 5 seeds = 25). AF3 labels each prediction as
seed-S_sample-N (5 seeds x 5 samples = 25). Both schemes produce 5
structural variants per seed value, AF2's "model" & AF3's "sample"
are the corresponding sub-index w/in a seed. Predictions are paired by:

    AF2 model_M_seed_S  <-->  AF3 seed-S_sample-(M-1)

Ex: AF2 "model_4_seed_001" <-> AF3 "seed-2_sample-3"
  (AF2 seed_001 -> AF3 seed-2; AF2 model_4 -> AF3 sample-3)

Every row includes:
  - Target protein
  - AF2_Label (real model/seed identifier, Ex: "model_4_seed_001")
  - AF2_TM_Score
  - AF3_Label (real seed/sample identifier, Ex: "seed-2_sample-3")
  - AF3_TM_Score
  - AF3_minus_AF2 (signed difference)
  - AF3_Higher (True/False)
  - PDB_Deposition_Date
  - Training_Access (which model(s) had this protein's structure available
    during training)

Usage:
    python3 build_model_comparison_tables.py \\
        --msa-swaps-root /data/AF_interpretation/msa_swaps \\
        --af3-singleseq-csv /data/schwartzlc/singleseq_af3_tmscores.csv \\
        --out-dir /data/schwartzlc
"""

# Import libraries.
import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.request
from datetime import datetime
from typing import Optional

# Structures deposited before this date are considered available to AF2
# during training.
AF2_CUTOFF = datetime(2018, 4, 30)
# Structures deposited before this date are considered available to AF3
# during training.
AF3_CUTOFF = datetime(2021, 9, 30)

# AF2 result directories containing the single-sequence results.
AF2_ROOTS = ["targets_60", "targets_100", "targets_150"]
# Directory names that shouldn't be interpreted as target proteins.
NON_TARGET_NAMES = {"temp_pdbs"}

# Targets included in the single-sequence/no-MSA experiment.
SINGLESEQ_TARGETS = [
    "1bor", "1bw6", "1cqu", "1dtk", "1fex", "1gb1", "1h7v", "1vex", "2fs1",
    "2jvv", "2k9h", "2mm4", "2nc2", "2ovo",
    "1c7m", "3l1x", "4g69", "4uz2", "5aim", "6ls1", "7qdv", "7y8h", "8ss1", "8vrc",
    "1de3", "1stb", "1z6u", "3h3i", "4ga0", "4zj9", "7cqf", "7h4z", "7yuj", "8jx4",
]

# Matches AF2 structure labels such as:
#   alphafold2_ptm_model_1_seed_004
#   alphafold2_model_1_seed_004
ENTRY_HEADER = re.compile(r"alphafold2_(?:ptm_)?model_(\d)_seed_(\d+)")
# Matches standalone TM-score lines in the AF2 TM-align output.
TMSCORE_STANDALONE = re.compile(r"^\s*(\d+\.\d+)\s*$", re.MULTILINE)


def is_probable_target_dir(name: str) -> bool:
    """
    Determine whether a directory name looks like an AF2 target directory.
    """
    # Ignore known temporary/non-target directories.
    if name in NON_TARGET_NAMES:
        return False
    # Ignore files & archives if they appear in the directory listing.
    if name.endswith((".fa", ".fasta", ".tar.gz", ".zip", ".png", ".csv")):
        return False
    # Ignore directories w/ the proper-MSA data.
    if "proper_msa" in name:
        return False
    return True


def fetch_deposition_date(pdb_id: str, max_retries: int = 5) -> str:
    """
    Query RCSB PDB REST API for a structure's deposition date, retrying
    on timeout/transient errors w/ exponential backoff before giving up.
    """
    # Build the RCSB API URL for this PDB entry.
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
    last_error = None
    # Try the API request up to max_retries times.
    for attempt in range(1, max_retries + 1):
        try:
            # Request the PDB metadata & decode the JSON response.
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read().decode())
            # Extract only the YYYY-MM-DD portion of the deposition date.
            return data.get("rcsb_accession_info", {}).get("deposit_date", "UNKNOWN")[:10]
        except Exception as e:
            last_error = e
            # Wait progressively longer BT failed attempts.
            if attempt < max_retries:
                wait = 2 ** attempt  # 2, 4, 8, 16 seconds.
                print(f"[RETRY {attempt}/{max_retries}] {pdb_id}: {e} -- "
                      f"retrying in {wait}s...", file=sys.stderr)
                time.sleep(wait)
    # If every attempt failed, record the problem & return UNKNOWN.
    print(f"[WARN] Could not fetch deposition date for {pdb_id} after "
          f"{max_retries} attempts: {last_error}", file=sys.stderr)
    return "UNKNOWN"


def training_access(date_str: str) -> str:
    """
    Classify whether the target structure was available during AF2 &/or
    AF3 training based on its PDB deposition date.

    Returns:
        "AF2 and AF3"  -> deposited before the AF2 cutoff
        "AF3 only"     -> deposited after AF2 but before AF3 cutoff
        "neither"      -> deposited after the AF3 cutoff
        "UNKNOWN"      -> missing/invalid deposition date
    """
    if date_str == "UNKNOWN":
        return "UNKNOWN"
    try:
        dt = datetime.strptime(date_str, "%Y-%m-%d")
    except ValueError:
        return "UNKNOWN"
    if dt < AF2_CUTOFF:
        return "AF2 and AF3"
    elif dt < AF3_CUTOFF:
        return "AF3 only"
    else:
        return "neither"


def af2_to_af3_key(model_num: str, seed_num: str) -> tuple[int, int]:
    """
    Map AF2's (model_M, seed_S) to the (af3_seed, af3_sample) key that
    should match it: af3_seed = AF2_seed + 1 (0-idx -> 1-idx),
    af3_sample = AF2_model - 1 (1-idx -> 0-idx).
    """
    af3_seed = int(seed_num) + 1
    af3_sample = int(model_num) - 1
    return (af3_seed, af3_sample)


def parse_af2_singleseq_labeled(path: str) -> dict:
    """
    Parse tmscores_ss_{target}_3.txt, returning
    {(af3_seed, af3_sample): (af2_label, score)}, keyed by the AF3
    seed/sample tuple that should match each AF2 entry.

    Handles 2 different header styles seen across target size classes:
      targets_100/150: ['1de3', 'rcsb_pdb_1DE3_..._model_1_seed_000']
      targets_60:       1bor,rcsb_pdb_1BOR_..._model_1_seed_002
    """
    # Missing files simply produce no AF2 entries for this target.
    if not os.path.isfile(path):
        return {}
    # Read the complete TM-align output file.
    with open(path, "r", errors="replace") as fh:
        content = fh.read()

    # Split the file into individual AF2 prediction entries.
    entries = re.split(r"(?=\['[\w]+', 'rcsb_pdb|^[\w]+,rcsb_pdb)", content, flags=re.MULTILINE)
    results = {}
    for entry in entries:
        # Find the AF2 model & seed numbers in this entry.
        header_match = ENTRY_HEADER.search(entry)
        if not header_match:
            continue
        model_num, seed_num = header_match.groups()
        # Create a clean AF2 label for the output CSV.
        af2_label = f"model_{model_num}_seed_{seed_num.zfill(3)}"
        # Extract TM-score values from the entry.
        score_lines = TMSCORE_STANDALONE.findall(entry)
        if not score_lines:
            continue
        # Use the final TM-score found in this entry.
        score = float(score_lines[-1])
        # Convert AF2's model/seed numbering to the corresponding
        # AF3 seed/sample numbering.
        key = af2_to_af3_key(model_num, seed_num)
        results[key] = (af2_label, score)
    return results


def load_af2_singleseq_by_target(msa_swaps_root: str) -> dict:
    """
    Load all AF2 single-sequence results.

    Return:
        {target: {(af3_seed, af3_sample): (af2_label, score)}}.
    """
    result = {}
    # Search all AF2 target-sized directories.
    for root_name in AF2_ROOTS:
        root_path = os.path.join(msa_swaps_root, root_name)
        if not os.path.isdir(root_path):
            continue
        for target in sorted(os.listdir(root_path)):
            target_path = os.path.join(root_path, target)
            if not os.path.isdir(target_path) or not is_probable_target_dir(target):
                continue
            # Locate the single-sequence TM-align output for this target.
            ss_file = os.path.join(target_path, "tmscores_ss", f"tmscores_ss_{target}_3.txt")
            entries = parse_af2_singleseq_labeled(ss_file)
            if entries:
                result[target] = entries
    return result


def load_af3_singleseq_by_target(csv_path: str) -> dict:
    """
    Load AF3 single-sequence TM-align results from CSV.

    Return:
        {target: {(seed, sample): (af3_label, score)}}.
    """
    result = {}
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            target = row["Target"]
            # Skip rows where the TM-score isn't numeric.
            try:
                score = float(row["TM_Score"])
            except (TypeError, ValueError):
                continue
            # Read AF3's seed/sample identifiers.
            seed_num = int(row["Seed_Number"])
            sample_num = int(row["Sample_Number"])
            label = f"seed-{seed_num}_sample-{sample_num}"
            # Store this score under its seed/sample key.
            result.setdefault(target, {})[(seed_num, sample_num)] = (label, score)
    return result


def build_comparison_rows(af2_by_target: dict, af3_by_target: dict,
                           dep_dates: dict) -> list[dict]:
    """
    Match AF2 & AF3 predictions target-by-target & create 1
    comparison row for every successfully seed-matched prediction.

    Only targets present in both AF2 & AF3 datasets are included.
    """
    rows = []
    # Only compare targets that have results from both models.
    for target in sorted(set(af2_by_target) & set(af3_by_target)):
        af2_dict = af2_by_target[target]  # {(seed,sample): (label, score)}
        af3_dict = af3_by_target[target]  # {(seed,sample): (label, score)}

        # Find predictions that exist in both datasets.
        matched_keys = sorted(set(af2_dict) & set(af3_dict))
        unmatched_af2 = set(af2_dict) - set(af3_dict) # Identify AF2 predictions with no AF3 counterpart.
        unmatched_af3 = set(af3_dict) - set(af2_dict) # Identify AF3 predictions with no AF2 counterpart.
        
        # Warn if the expected 25-to-25 matching is incomplete.
        if unmatched_af2 or unmatched_af3:
            print(f"[WARN] {target}: {len(unmatched_af2)} AF2 entries and "
                  f"{len(unmatched_af3)} AF3 entries had no seed-matched "
                  f"counterpart (matched {len(matched_keys)}/25).", file=sys.stderr)
        
        # Get the PDB deposition date and infer training access.
        dep_date = dep_dates.get(target, "UNKNOWN")
        access = training_access(dep_date)

        # Create 1 row for every AF2/AF3 prediction pair.
        for key in matched_keys:
            af2_label, af2_score = af2_dict[key]
            af3_label, af3_score = af3_dict[key]
            rows.append({
                "Target": target,
                "AF2_Label": af2_label,
                "AF2_TM_Score": af2_score,
                "AF3_Label": af3_label,
                "AF3_TM_Score": af3_score,
                "AF3_minus_AF2": round(af3_score - af2_score, 5),
                "AF3_Higher": af3_score > af2_score,
                "PDB_Deposition_Date": dep_date,
                "Training_Access": access,
            })
    return rows


def write_tiered_csvs(rows: list[dict], out_dir: str, experiment_name: str):
    """
    Write 3 CSVs for an experiment:

        {experiment}_all.csv
        {experiment}_above_0.5.csv
        {experiment}_around_0.5.csv

    The filtering is based on AF3's TM-score.
    """
    fieldnames = ["Target", "AF2_Label", "AF2_TM_Score",
                  "AF3_Label", "AF3_TM_Score", "AF3_minus_AF2", "AF3_Higher",
                  "PDB_Deposition_Date", "Training_Access"]

    def write_csv(subset_rows, suffix):
        """
        Write 1 filtered subset to its corresponding CSV file.
        """
        path = os.path.join(out_dir, f"{experiment_name}_{suffix}.csv")
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(subset_rows)
        print(f"  Wrote {len(subset_rows)} row(s): {path}")

    write_csv(rows, "all")
    write_csv([r for r in rows if r["AF3_TM_Score"] > 0.5], "above_0.5")
    write_csv([r for r in rows if 0.4 <= r["AF3_TM_Score"] <= 0.6], "around_0.5")

    if rows:
        n_af3_higher = sum(1 for r in rows if r["AF3_Higher"])
        print(f"  AF3 higher than AF2 (seed-matched): {n_af3_higher}/{len(rows)} "
              f"({100*n_af3_higher/len(rows):.1f}%)")
    else:
        print("  No rows.")


def main():
    # Define command-line arguments.
    ap = argparse.ArgumentParser()
    ap.add_argument("--msa-swaps-root", required=True)
    ap.add_argument("--af3-singleseq-csv", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True) # Create the output directory if it doesn't already exist.

    # Fetch each single-sequence target's PDB deposition date once.
    all_targets = sorted(set(SINGLESEQ_TARGETS))
    print(f"Fetching PDB deposition dates for {len(all_targets)} unique targets...")
    dep_dates = {}
    # Query RCSB for each target's deposition date.
    for target in all_targets:
        dep_dates[target] = fetch_deposition_date(target)
        time.sleep(0.2)

    # Safety net: re-attempt any target that still came back UNKNOWN after
    # its own retries, in case of a longer-lasting network issue.
    still_unknown = [t for t, d in dep_dates.items() if d == "UNKNOWN"]
    if still_unknown:
        print(f"\n{len(still_unknown)} target(s) still UNKNOWN after first pass; "
              f"retrying: {still_unknown}")
        time.sleep(3)
        for target in still_unknown:
            dep_dates[target] = fetch_deposition_date(target, max_retries=5)
            time.sleep(0.5)

    # Report any targets whose deposition dates remain unavailable.
    final_unknown = [t for t, d in dep_dates.items() if d == "UNKNOWN"]
    if final_unknown:
        print(f"\n[WARN] Still could not fetch deposition dates for: {final_unknown} "
              f"(will show as UNKNOWN/UNKNOWN in output CSVs)")
    print("Done.\n")

    # --- Single-sequence experiment ---
    print("=== Single-sequence (no-MSA) experiment ===")
    af2_ss = load_af2_singleseq_by_target(args.msa_swaps_root)
    af3_ss = load_af3_singleseq_by_target(args.af3_singleseq_csv)
    print(f"AF2 targets with parsed labels: {len(af2_ss)}")
    print(f"AF3 targets: {len(af3_ss)}")

    ss_rows = build_comparison_rows(af2_ss, af3_ss, dep_dates)
    write_tiered_csvs(ss_rows, args.out_dir, "singleseq")

    print("\nDone. All CSVs include real AF2/AF3 labels, PDB deposition date, "
          "& which model(s) had training access to that structure.")


if __name__ == "__main__":
    main()
