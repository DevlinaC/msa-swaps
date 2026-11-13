#!/usr/bin/env python3
"""Item 2: GREMLIN contact overlap at sequence separations 6, 12, and 24.

Reanalyze saved scores, retaining historical top-100/z-score/gap/offset rules.
Separation 5 is a mandatory regression check against all 362 saved pairs.
Run locally with one thread or within an allocated Biowulf compute job.
"""
import os
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import jaccard as scipy_jaccard

FIXED = {"1bor", "1bw6", "1cqu", "1gb1", "1vex", "2fs1", "2nc2", "2ovo"}
SEPARATIONS = (5, 6, 12, 24)
ROOT = Path(__file__).resolve().parents[2]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gap_mask(path):
    """Stream aligned FASTA records, retaining the historical '-' gap rule."""
    counts = None
    n = 0
    sequence = []
    header_seen = False

    def consume(parts):
        nonlocal counts, n
        text = "".join(parts).upper()
        require(bool(text), f"Empty FASTA record: {path}")
        require(text.isascii(), f"Non-ASCII alignment: {path}")
        values = np.frombuffer(text.encode("ascii"), dtype=np.uint8)
        if counts is None:
            counts = np.zeros(len(values), dtype=np.int64)
        require(len(values) == len(counts), f"Unequal alignment lengths: {path}")
        counts += values == ord("-")
        n += 1

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header_seen:
                    consume(sequence)
                sequence = []
                header_seen = True
            else:
                require(header_seen, f"Sequence before FASTA header: {path}")
                sequence.append(line)
    require(header_seen, f"No FASTA records: {path}")
    consume(sequence)
    return counts / n < 0.25, n


def read_scores(path, mask):
    frame = pd.read_csv(path)
    require({"i", "j", "zscore"}.issubset(frame.columns), f"Bad schema: {path}")
    require(len(frame) > 0, f"Empty score table: {path}")
    values = frame[["i", "j", "zscore"]].to_numpy(dtype=float)
    indices = values[:, :2]
    require(np.isfinite(values).all(), f"Nonfinite scores/indices: {path}")
    require(np.equal(indices, np.floor(indices)).all(), f"Noninteger indices: {path}")
    require((indices >= 0).all() and (indices < len(mask)).all(), f"Indices outside alignment: {path}")
    require((indices[:, 0] < indices[:, 1]).all(), f"Expected upper triangle: {path}")
    require(not frame.duplicated(["i", "j"]).any(), f"Duplicate residue pairs: {path}")
    return frame


def select_contacts(frame, mask, separation):
    """Match pandas' historical descending sort, including cutoff tie order."""
    idx = frame[["i", "j"]].to_numpy(dtype=int)
    keep = mask[idx[:, 0]] & mask[idx[:, 1]]
    filtered = frame.loc[keep & (frame.zscore >= 1.0)]
    eligible = filtered.loc[abs(filtered.i - filtered.j) >= separation]
    ordered = eligible.sort_values("zscore", ascending=False, kind="quicksort")
    selected = ordered.head(100)
    contacts = set(map(tuple, selected[["i", "j"]].to_numpy(dtype=int)))
    cutoff_ties = 0
    if len(ordered) > 100:
        boundary = ordered.iloc[99].zscore
        if ordered.iloc[100].zscore == boundary:
            cutoff_ties = int((ordered.zscore == boundary).sum())
    return contacts, len(eligible), cutoff_ties


def historical_jaccard(a, b):
    return len(a & b) / len(a | b) if a and b else 0.0


def compare(a, b, optimize):
    selected, offset, best = a, 0, 0.0
    if optimize:
        for candidate in range(-5, 6):
            shifted = {(i + candidate, j + candidate) for i, j in a}
            value = historical_jaccard(shifted, b)
            if value > best:
                selected, offset, best = shifted, candidate, value
    union = selected | b
    historical = historical_jaccard(selected, b)
    library_error = np.nan
    if union:
        keys = sorted(union)
        library = 1 - scipy_jaccard([k in selected for k in keys], [k in b for k in keys])
        library_error = abs(historical - library)
        require(library_error < 1e-12, "Manual and SciPy Jaccard disagree")
    status = ("both_empty" if not union else "swapped_empty" if not selected
              else "native_empty" if not b else "ok")
    return dict(jaccard=historical if union else np.nan,
                historical_jaccard=historical, zero_offset_jaccard=historical_jaccard(a, b) if a | b else np.nan,
                offset=offset, contacts_swapped=len(selected), contacts_native=len(b),
                intersection=len(selected & b), union=len(union), status=status,
                scipy_abs_error=library_error)


def summarize(frame):
    rows = []
    for (separation, reference), group in frame.groupby(["separation", "reference"]):
        for length_class in ("all", 60, 100, 150):
            sub = group if length_class == "all" else group[group.length_class == length_class]
            valid = sub.jaccard.dropna()
            rows.append(dict(separation=separation, reference=reference, length_class=length_class,
                n_pairs=len(sub), n_defined=len(valid), n_both_empty=int((sub.status == "both_empty").sum()),
                n_swapped_empty=int((sub.contacts_swapped == 0).sum()),
                n_native_empty=int((sub.contacts_native == 0).sum()),
                mean=valid.mean(), median=valid.median(), q25=valid.quantile(.25), q75=valid.quantile(.75),
                mean_selected_swapped=sub.contacts_swapped.mean(), mean_selected_native=sub.contacts_native.mean(),
                n_nonzero_offset=int((sub.offset != 0).sum())))
    return pd.DataFrame(rows)


def make_plot(frame, out):
    plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.7), sharey=True)
    for ax, reference, title, color in zip(axes, ("A", "B"),
            ("Query A", "MSA donor B"), ("#31688e", "#b45a2a")):
        values = [frame.loc[(frame.separation == s) & (frame.reference == reference), "jaccard"].dropna().to_numpy()
                  for s in (6, 12, 24)]
        boxes = ax.boxplot(values, tick_labels=["≥6", "≥12", "≥24"], patch_artist=True,
                           widths=.45, showfliers=True,
                           flierprops=dict(marker=".", markersize=3, alpha=.35),
                           medianprops=dict(color="black", linewidth=1.6))
        for box in boxes["boxes"]:
            box.set_facecolor(color)
            box.set_alpha(.5)
        for x, v in enumerate(values, 1):
            ax.text(x, 1.035, f"n={len(v)}", ha="center", fontsize=9)
        ax.set_title(title, pad=22)
        ax.set_xlabel("Minimum sequence separation |i − j|")
        ax.set_ylim(-.04, 1.09)
        ax.set_yticks(np.arange(0, 1.01, .2))
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=.15)
    axes[0].set_ylabel("GREMLIN contact-set Jaccard")
    fig.suptitle("GREMLIN overlap across sequence-separation cutoffs", y=.99, fontsize=14)
    fig.text(.5, .025, "Up to 100 contacts; z-score ≥1; column gap fraction <25%.\n"
             "A: offset 0. B: best offset −5…+5. Boxes: IQR; line: median; whiskers: 1.5×IQR.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .11, 1, .93))
    fig.savefig(out / "gremlin_separation.png", dpi=300)
    fig.savefig(out / "gremlin_separation.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT)
    parser.add_argument("--native-alignments", type=Path, help="Override candidate native gap-mask alignments")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "for_github_repo/analysis/gremlin_separation")
    args = parser.parse_args()
    root, out = args.data_root.resolve(), args.out_dir.resolve()
    require(not out.exists(), f"Output directory already exists; choose a new --out-dir: {out}")
    native_root = args.native_alignments or root / "original_proper_msas_af_sampling"
    sources = {}

    def track(path):
        path = Path(path).resolve()
        require(path.is_file() and path.stat().st_size > 0, f"Missing or empty input: {path}")
        key = str(path.relative_to(root)) if path.is_relative_to(root) else str(path)
        if key not in sources:
            sources[key] = dict(sha256=sha256(path), bytes=path.stat().st_size)
        return path

    saved = pd.read_csv(track(root / "for_github_repo/TM_Jaccard_codes/jaccard_results_gremlin.csv"))
    targets = pd.read_csv(track(root / "for_github_repo/msa-swaps-github/data/targets.csv"))
    targets["pdb_id"] = targets.pdb_id.str.lower()
    require(len(targets) == 34 and targets.pdb_id.nunique() == 34, "Expected 34 unique targets")
    classes = dict(zip(targets.pdb_id, targets.length_class))
    require({int(k): int(v) for k, v in targets.length_class.value_counts().items()} == {60: 14, 100: 10, 150: 10},
            "Unexpected size-class design")
    expected = {(a, b) for a in classes for b in classes if a != b and classes[a] == classes[b]}
    actual = set(zip(saved.Protein_A, saved.Protein_B))
    require(len(saved) == len(actual) == 362 and actual == expected, "Pair manifest is incomplete or duplicated")
    require(np.isfinite(saved[["Jaccard_vs_ProteinA", "Jaccard_vs_ProteinB", "Best_Offset_ProteinB"]].to_numpy()).all(),
            "Invalid historical result values")
    audit, baseline, results = [], [], []

    def load(score_path, alignment_path, identifier, condition):
        mask, n_sequences = gap_mask(track(alignment_path))
        frame = read_scores(track(score_path), mask)
        selections = {}
        for sep in SEPARATIONS:
            contacts, n_eligible, ties = select_contacts(frame, mask, sep)
            selections[sep] = contacts
            audit.append(dict(identifier=identifier, condition=condition, separation=sep,
                score_file=str(score_path.relative_to(root)), alignment_file=str(alignment_path),
                n_sequences=n_sequences, alignment_columns=len(mask), retained_columns=int(mask.sum()),
                score_rows=len(frame), n_eligible=n_eligible, n_selected=len(contacts), boundary_tie_rows=ties))
        return selections

    native = {a: load(root / f"gremlin/results_gremlin_proper/full_{a}.csv",
                     native_root / a / f"{a}.fas", a, "native") for a in sorted(classes)}
    print("Validated and selected contacts for 34 native inputs.", flush=True)
    for number, row in enumerate(saved.sort_values(["Protein_A", "Protein_B"]).itertuples(index=False), 1):
        a, b = row.Protein_A, row.Protein_B
        pair = f"{a}_{b}_aligned"
        fixed = b == "2mm4" and a in FIXED
        alignment = root / (f"2mm4_fixes/{pair}/{pair}.fas" if fixed else f"original_junk_msas/{a}/{pair}/{pair}.fas")
        swapped = load(root / f"gremlin/results_gremlin_junk/full_{pair}.csv", alignment, pair, "swapped")
        for sep in SEPARATIONS:
            for ref, protein in (("A", a), ("B", b)):
                comparison = compare(swapped[sep], native[protein][sep], ref == "B")
                results.append(dict(pair=pair, protein_a=a, protein_b=b, length_class=classes[a],
                                    corrected=fixed, separation=sep, reference=ref, **comparison))
                if sep == 5:
                    previous = getattr(row, f"Jaccard_vs_Protein{ref}")
                    error = abs(comparison["historical_jaccard"] - previous)
                    old_offset = int(row.Best_Offset_ProteinB) if ref == "B" else 0
                    baseline.append(dict(pair=pair, reference=ref, saved=previous,
                        reproduced=comparison["historical_jaccard"], absolute_error=error,
                        saved_offset=old_offset, reproduced_offset=comparison["offset"]))
                    require(error <= 1e-10 and comparison["offset"] == old_offset,
                            f"Historical regression failed: {pair} {ref}")
        if number % 60 == 0 or number == 362:
            print(f"Processed {number}/362 swapped inputs.", flush=True)

    frame = pd.DataFrame(results)
    require(len(frame) == 362 * 4 * 2, "Incomplete output comparisons")
    require(not frame.duplicated(["pair", "separation", "reference"]).any(), "Duplicate output comparisons")
    summary = summarize(frame)
    directions = []
    for sep in SEPARATIONS:
        wide = frame[frame.separation == sep].pivot(index="pair", columns="reference", values="jaccard")
        valid = wide.dropna()
        directions.append(dict(separation=sep, n_both_defined=len(valid),
            n_B_gt_A=int((valid.B > valid.A).sum()), n_equal=int((valid.B == valid.A).sum()),
            n_A_gt_B=int((valid.A > valid.B).sum()), n_undefined=362-len(valid)))
    track(Path(__file__))
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(), data_root=str(root),
        native_alignment_root=str(native_root), python=platform.python_version(), platform=platform.platform(),
        versions={"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        parameters=dict(separations=list(SEPARATIONS), baseline=5, requested=[6,12,24], top_n=100,
            min_zscore=1.0, gap_fraction_strictly_below=.25, offsets_A=[0], offsets_B=list(range(-5,6)),
            selection_order="gap mask, z-score, separation, descending pandas quicksort, first 100",
            offset_ties="start at zero; update only for strictly greater score while scanning -5 through +5",
            empty_union="undefined Jaccard; historical_jaccard retains legacy zero for regression",
            one_empty="Jaccard zero", threads=1),
        limitations=["Native alignment candidates match prior audit but complete historical row identity is unconfirmed.",
            "Donor B offsets are optimized; A is not. No new biological residue alignment is inferred.",
            "Offsets retain historical behavior, including shifted indices outside a native alignment.",
            "Pairs share proteins; summaries and boxplots are descriptive, without independent-replicate inference."],
        source_files=sources, validation=dict(baseline_comparisons=len(baseline),
            maximum_baseline_error=max(r["absolute_error"] for r in baseline), offset_mismatches=0,
            maximum_scipy_error=float(frame.scipy_abs_error.max())))
    out.mkdir(parents=True)
    for name, table in (("comparisons.csv", frame), ("summary.csv", summary),
                        ("contact_selection_audit.csv", pd.DataFrame(audit)),
                        ("baseline_validation.csv", pd.DataFrame(baseline)),
                        ("donor_query_direction.csv", pd.DataFrame(directions))):
        table.to_csv(out / name, index=False)
    make_plot(frame, out)
    provenance["output_sha256"] = {p.name: sha256(p) for p in sorted(out.iterdir()) if p.is_file()}
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    # Written last: an interrupted/failed plot or output write cannot mark the run complete.
    (out / "COMPLETE.json").write_text(json.dumps({"comparisons": len(frame), "pairs": 362,
        "provenance_sha256": sha256(out / "provenance.json")}, indent=2) + "\n")
    print(summary[summary.length_class == "all"].to_string(index=False))
    print(pd.DataFrame(directions).to_string(index=False))
    print(f"Results: {out}")


if __name__ == "__main__":
    main()
