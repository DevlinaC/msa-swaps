#!/usr/bin/env python3
"""
plot_4_violins.py

Builds 1 figure w/ 4 violins, comparing swapped-MSA
(junk) predictions against Protein B (the MSA source) across 4
different metrics:

  1. AF3 TM-score vs Protein B
     from tm_align_full_results.csv: Job_Type == "junk",
     TM_Score where Compared_Protein == partner (Protein B).

  2. AF2 TM-score vs Protein B
     from targets_60/100/150 tmscores_junk_{target}_{partner}/
         tmscores_junk_{partner}_3.txt files (bare-float lines,
         25 scores per pair; partner's file = vs Protein B).

  3. GREMLIN Jaccard_vs_ProteinB
     from jaccard_results_gremlin.csv.

  4. MSA Transformer Jaccard_vs_ProteinB
     from jaccard_results_msatransformer.csv.

All 4 metrics are bounded in [0, 1], so they share 1 y-axis.

Usage:
    python3 plot_4_violins.py \
        --msa-swaps-root /data/AF_interpretation/msa_swaps \
        --tm-csv /data/schwartzlc/tm_align_full_results.csv \
        --gremlin-csv /data/schwartzlc/jaccard_results_gremlin.csv \
        --msatransformer-csv /data/schwartzlc/jaccard_results_msatransformer.csv \
        --out-dir /data/schwartzlc/plots
"""

import argparse # For command-line arguments.
import csv # For reading CSV files.
import os # For working w/ files & folders.
import re # For using regular expressions (pattern matching).
import sys # For printing error messages to stderr.

import matplotlib # For making graphs.
matplotlib.use("Agg")  # Headless-safe backend for Biowulf.
import matplotlib.pyplot as plt # Get plotting funcs from matplotlib.
import numpy as np # For numerical operations & random # generation.

AF2_ROOTS = ["targets_60", "targets_100", "targets_150"] # 3 folders w/ the AF2 prediction results.
TMSCORE_LINE = re.compile(r"^\s*(\d+\.\d+)\s*$", re.MULTILINE)
NON_TARGET_NAMES = {"temp_pdbs", "tmscores_ss"} # Ignore bc not protein tm score folders.


def is_probable_target_dir(name: str) -> bool:
    """
    Decide whether a folder is probably a real target directory.
    Returns True if it looks like one. Else, returns False.
    """
    # Ignore known folders that aren't protein targets.
    if name in NON_TARGET_NAMES:
        return False
    # Ignore common files instead of folders.
    if name.endswith((".fa", ".fasta", ".tar.gz", ".zip", ".png", ".csv")):
        return False
    # Ignore folders related to proper MSA processing.
    if "proper_msa" in name:
        return False
    # Otherwise it looks like a valid target folder.
    return True


def parse_tmscores_from_file(path: str) -> list[float]:
    """
    Reads a text file & extracts every TM-score found inside.
    Returns a list of floating-pt #s.
    """
    # If the file doesn't exist, return an empty list.
    if not os.path.isfile(path):
        return []
    # Open the file for reading.
    # errors'"replace" prevents crashes if strange characters appear.
    with open(path, "r", errors="replace") as fh:
        content = fh.read()
    # Find every decimal # matching my regex, convert each string into a float,
    # & return the list.
    return [float(m) for m in TMSCORE_LINE.findall(content)]


def load_af2_vs_b(msa_swaps_root: str) -> list[float]:
    """
    Walk targets_60/100/150 & collect AF2 TM-scores vs Protein B only
    (the tmscores_junk_{partner}_3.txt file in each pair folder).
    """
    af2_b = [] # Store all TM-scores here.
    n_pairs_found = 0 # Count successful protein pairs.
    n_pairs_missing_files = 0 # Count missing files.

    # Loop over: targets_60, targets_100, targets_150.
    for root_name in AF2_ROOTS:
        # Build the complete folder path.
        root_path = os.path.join(msa_swaps_root, root_name)
        # Skip if the folder doesn't exist.
        if not os.path.isdir(root_path):
            print(f"[WARN] AF2 root not found: {root_path}", file=sys.stderr)
            continue

        # Look at every folder inside.
        for target in sorted(os.listdir(root_path)):
            target_path = os.path.join(root_path, target)
            # Skip anything that isn't a valid target directory.
            if not os.path.isdir(target_path) or not is_probable_target_dir(target):
                continue
            # Folder names start w/ tmscores_junk<target>_.
            tm_prefix = f"tmscores_junk_{target}_"
            # Search every folder inside this target.
            for entry in sorted(os.listdir(target_path)):
                entry_path = os.path.join(target_path, entry)
                # Skip folders that don't match the naming pattern.
                if not os.path.isdir(entry_path) or not entry.startswith(tm_prefix):
                    continue

                # Everything after the prefix is the partner protein.
                partner = entry[len(tm_prefix):]
                # Skip if somehow empty.
                if not partner:
                    continue
                
                # Build the expected filename.
                file_b = os.path.join(entry_path, f"tmscores_junk_{partner}_3.txt")
                scores_b = parse_tmscores_from_file(file_b)

                # Skip if nothing was found.
                if not scores_b:
                    n_pairs_missing_files += 1
                    continue

                af2_b.extend(scores_b) # Add every score into the master list.
                n_pairs_found += 1 # Count successful pair.

    # Print summary.
    print(f"AF2: found {n_pairs_found} target/partner pairs "
          f"({n_pairs_missing_files} skipped, missing files).")
    print(f"AF2: {len(af2_b)} protein-B scores loaded.")
    return af2_b


def load_af3_vs_b(tm_csv_path: str) -> list[float]:
    """
    Extract AF3 TM-scores for junk jobs, Compared_Protein == partner (Protein B).
    """
    scores = []
    n_unmatched = 0
    
    # Open CSV file.
    with open(tm_csv_path, newline="") as fh:
        # Read each row as a dict.
        reader = csv.DictReader(fh)
        # Process every row.
        for row in reader:
            # Only use junk jobs.
            if row.get("Job_Type") != "junk":
                continue
            
            # Ignore unexpected labels.
            job_label = row["Job_Label"]
            suffix = "_aligned"
            if not job_label.endswith(suffix):
                n_unmatched += 1
                continue
            # Remove "_aligned".
            core = job_label[: -len(suffix)]
            # Split: target_partner.
            target, _, partner = core.partition("_")
            # Skip malformed labels.
            if not target or not partner:
                n_unmatched += 1
                continue
            
            # Keep only comparisons against Protein B.
            if row["Compared_Protein"] != partner:
                continue  # This row is the vs-A comparison; skip it.
            
            # Convert score into float.
            try:
                scores.append(float(row["TM_Score"]))
            except (TypeError, ValueError):
                continue

    print(f"AF3: {len(scores)} scores loaded vs Protein B "
          f"({n_unmatched} rows unmatched/skipped).")
    return scores


def load_jaccard_vs_b(csv_path: str, label: str) -> list[float]:
    """
    Read a CSV file w/ Jaccard scores.
    Return:
        Every valid Jaccard score.
    """
    scores = []
    n_skipped = 0

    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            raw = row.get("Jaccard_vs_ProteinB", "")
            try:
                scores.append(float(raw))
            except (TypeError, ValueError):
                n_skipped += 1

    print(f"{label}: {len(scores)} scores loaded vs Protein B "
          f"({n_skipped} rows skipped/unparseable).")
    return scores


def make_violin_figure(af3_scores, af2_scores, gremlin_scores, msat_scores, out_path):
    """
    Make 1 violin plot comparing the 4 different scoring methods.
    """
    plt.rcParams.update({
        "font.size": 12,
        "font.family": "sans-serif",
        "axes.edgecolor": "#333333",
        "axes.linewidth": 1.0,
    })

    data = [af3_scores, af2_scores, gremlin_scores, msat_scores]
    # Labels shown under each violin.
    labels = [
        f"AF3 TM-score\n(n={len(af3_scores)})",
        f"AF2 TM-score\n(n={len(af2_scores)})",
        f"GREMLIN\nJaccard (n={len(gremlin_scores)})",
        f"MSA Transformer\nJaccard (n={len(msat_scores)})",
    ]
    # Color for each violin.
    colors = ["#DD8452", "#4C72B0", "#55A868", "#8172B2"]

    # Create figure.
    fig, ax = plt.subplots(figsize=(9, 6), dpi=300)

    positions = [1, 2, 3, 4]
    # Draw violin plots.
    parts = ax.violinplot(
        data, positions=positions, widths=0.7,
        showmeans=False, showmedians=True, showextrema=True,
    )

    # Color each violin.
    for pc, color in zip(parts["bodies"], colors):
        pc.set_facecolor(color)
        pc.set_alpha(0.55)
        pc.set_edgecolor("#333333")
        pc.set_linewidth(1.0)

    # Style median & whiskers.
    for key in ("cbars", "cmins", "cmaxes", "cmedians"):
        if key in parts:
            parts[key].set_edgecolor("#333333")
            parts[key].set_linewidth(1.2)

    # Random generator for jitter.
    rng = np.random.default_rng(seed=42)
    # Plot every individual point.
    for i, scores in zip(positions, data):
        if not scores:
            continue
        # Small horizontal random offsets.
        jitter = rng.uniform(-0.08, 0.08, size=len(scores))
        ax.scatter(
            np.full(len(scores), i) + jitter, scores,
            s=6, alpha=0.15, color="#333333", edgecolors="none",
            zorder=2,
        )

    # Label x-axis.
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    # Label -axis.
    ax.set_ylabel("Simularity to Native MSA Structure")
    # Figure title.
    ax.set_title("Swapped-MSA (junk) predictions vs. Protein B\nacross four metrics")
    # Force y-axis BT 0 & 1.
    ax.set_ylim(-0.02, 1.02)
    # Add horizontal grid lines.
    ax.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
    # Remove unnecessary borders.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Prevent labels from overlapping.
    fig.tight_layout()
    # Save figure.
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    # Close figure to free memory.
    plt.close(fig)
    print(f"\nSaved: {out_path}")


def main():
    # Create command-line parser.
    ap = argparse.ArgumentParser(
        description="Build a 4-violin figure comparing AF3 TM-score, AF2 TM-score, "
                    "GREMLIN Jaccard, and MSA Transformer Jaccard, all vs Protein B."
    )
    # Required command-line arguments.
    ap.add_argument("--msa-swaps-root", required=True,
                     help="Path containing targets_60/, targets_100/, targets_150/ "
                          "(AF2 results)")
    ap.add_argument("--tm-csv", required=True,
                     help="Path to tm_align_full_results.csv (AF3 results)")
    ap.add_argument("--gremlin-csv", required=True,
                     help="Path to jaccard_results_gremlin.csv")
    ap.add_argument("--msatransformer-csv", required=True,
                     help="Path to jaccard_results_msatransformer.csv")
    ap.add_argument("--out-dir", required=True,
                     help="Directory to save the output PNG")
    args = ap.parse_args()

    # Make output folder if needed.
    os.makedirs(args.out_dir, exist_ok=True)

    # Load every dataset.
    af3_scores = load_af3_vs_b(args.tm_csv)
    af2_scores = load_af2_vs_b(args.msa_swaps_root)
    gremlin_scores = load_jaccard_vs_b(args.gremlin_csv, "GREMLIN")
    msat_scores = load_jaccard_vs_b(args.msatransformer_csv, "MSA Transformer")

    # Warn if any dataset is empty.
    if not (af3_scores and af2_scores and gremlin_scores and msat_scores):
        print("[ERROR] One or more metrics returned no data; check CSV/path inputs.",
              file=sys.stderr)

    # Output filename.
    out_path = os.path.join(args.out_dir, "4_violins_proteinB.png")
    # Make the plot.
    make_violin_figure(af3_scores, af2_scores, gremlin_scores, msat_scores, out_path)


if __name__ == "__main__":
    main()