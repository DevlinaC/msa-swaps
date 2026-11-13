#!/usr/bin/env python3
"""Item 3: paired-score Pearson R over the union of top-L contact sets.

Both GREMLIN and MSA Transformer; separations 6, 12, 24. GREMLIN
reuses item 2 offsets (never maximizes R), with zero-offset diagnostics.
"""
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"

import argparse
import gzip
import csv
import json
import platform
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy
from scipy.stats import pearsonr, NearConstantInputWarning

import analyze_gremlin_separation as g

SEPARATIONS = (6, 12, 24)


def correlation(x, y):
    """Return undefined, never zero, for undefined Pearson coefficients."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    g.require(x.shape == y.shape and x.ndim == 1, "Unpaired Pearson inputs")
    g.require(np.isfinite(x).all() and np.isfinite(y).all(), "Nonfinite Pearson inputs")
    if len(x) < 2:
        return np.nan, "too_few_pairs", np.nan
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan, "constant_vector", np.nan
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", NearConstantInputWarning)
        r = float(pearsonr(x, y).statistic)
    if any(issubclass(w.category, NearConstantInputWarning) for w in caught):
        return np.nan, "near_constant_vector", np.nan
    independent = float(np.corrcoef(x, y)[0, 1])
    error = abs(r - independent)
    g.require(np.isfinite(r) and error <= 1e-10, "Pearson validation failed")
    return r, "ok", error


def top_union(x, y, length):
    """Rank within a shared finite-pair universe; retain all cutoff ties."""
    g.require(length > 0 and len(x) == len(y), "Invalid top-L inputs")
    all_r, all_status, all_error = correlation(x, y)
    result = dict(L=length, n_eligible=len(x), all_eligible_r=all_r,
                  all_eligible_status=all_status, all_eligible_numpy_error=all_error)
    if len(x) < length:
        return dict(result, pearson_r=np.nan, status="insufficient_eligible_for_L",
                    n_top_swapped=0, n_top_native=0, n_union=0, n_intersection=0,
                    cutoff_swapped=np.nan, cutoff_native=np.nan, numpy_abs_error=np.nan), None
    cut_x, cut_y = np.sort(x)[-length], np.sort(y)[-length]
    top_x, top_y = x >= cut_x, y >= cut_y
    union = top_x | top_y
    r, status, error = correlation(x[union], y[union])
    return dict(result, pearson_r=r, status=status, n_top_swapped=int(top_x.sum()),
                n_top_native=int(top_y.sum()), n_union=int(union.sum()),
                n_intersection=int((top_x & top_y).sum()), cutoff_swapped=float(cut_x),
                cutoff_native=float(cut_y), numpy_abs_error=error), (union, top_x, top_y)


def gremlin_pairs(swapped, native, swap_len, native_len, offset, separation):
    """Match identical coordinates after a fixed shift; no missing-score fill."""
    left = swapped.copy()
    left[["i", "j"]] += offset
    joined = left.merge(native, on=["i", "j"], how="inner", suffixes=("_swapped", "_native"), validate="one_to_one")
    joined = joined[joined.j - joined.i >= separation].sort_values(["i", "j"])
    start, end = max(0, offset), min(native_len, swap_len + offset)
    length = end - start
    g.require(length > 0, "No GREMLIN shared coordinate span")
    if len(joined):
        g.require((joined.i >= start).all() and (joined.j < end).all(), "Joined pair outside shared span")
    return (joined[["i", "j"]].to_numpy(dtype=int),
            joined.zscore_swapped.to_numpy(), joined.zscore_native.to_numpy(), length)


def matrix_pairs(a, b, separation):
    length = min(len(a), len(b))
    start_a, start_b = (len(a) - length + 1) // 2, (len(b) - length + 1) // 2
    aa = a[start_a:start_a+length, start_a:start_a+length]
    bb = b[start_b:start_b+length, start_b:start_b+length]
    i, j = np.triu_indices(length, k=separation)
    return np.column_stack((i, j)), aa[i, j], bb[i, j], length, start_a, start_b


def summary_table(frame):
    records = []
    for (method, separation, reference), group in frame.groupby(["method", "separation", "reference"]):
        for length_class in ("all", 60, 100, 150):
            sub = group if length_class == "all" else group[group.length_class == length_class]
            values = sub.pearson_r.dropna()
            records.append(dict(method=method, separation=separation, reference=reference,
                length_class=length_class, n_pairs=len(sub), n_defined=len(values),
                n_undefined=len(sub)-len(values), mean_r=values.mean(), median_r=values.median(),
                q25=values.quantile(.25), q75=values.quantile(.75), min_r=values.min(), max_r=values.max(),
                n_negative=int((values < 0).sum()), mean_union=sub.n_union.mean(),
                mean_all_eligible_r=sub.all_eligible_r.mean(),
                mean_zero_offset_r=sub.zero_offset_r.mean()))
    return pd.DataFrame(records)


def figure(frame, item2, out):
    """Integrated figure: item 2 Jaccard above item 3 Pearson distributions."""
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(3, 2, figsize=(9, 10.5))
    panels = [(item2, "jaccard", "GREMLIN Jaccard · up to 100 contacts"),
              (frame[frame.method == "GREMLIN"], "pearson_r", "GREMLIN Pearson R · union of top-L sets"),
              (frame[frame.method == "MSA Transformer"], "pearson_r", "MSA Transformer Pearson R · union of top-L sets")]
    for row, (data, metric, label) in enumerate(panels):
        for col, ref in enumerate(("A", "B")):
            ax = axes[row, col]
            values = [data.loc[(data.separation == sep) & (data.reference == ref), metric].dropna().to_numpy()
                      for sep in SEPARATIONS]
            box = ax.boxplot(values, tick_labels=["≥6", "≥12", "≥24"], widths=.45, patch_artist=True,
                flierprops=dict(marker=".", markersize=2.5, alpha=.35),
                medianprops=dict(color="black", linewidth=1.5))
            for patch in box["boxes"]:
                patch.set_facecolor(("#31688e", "#b45a2a")[col]); patch.set_alpha(.5)
            for x, vals in enumerate(values, 1):
                ax.text(x, 1.055 if row else 1.035, f"n={len(vals)}", ha="center", fontsize=8)
            ax.set_ylim((-.05,1.12) if row == 0 else (-1.05,1.18))
            ax.set_yticks([0,.25,.5,.75,1] if row == 0 else [-1,-.5,0,.5,1])
            if row:
                ax.axhline(0, color="#888888", linewidth=.7, linestyle="--")
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="y", alpha=.15)
            ax.set_title(f"{'ABCDEF'[row*2+col]}  {'Query A' if ref=='A' else 'MSA donor B'}", loc="left", fontsize=12)
            if col == 0:
                ax.set_ylabel("Jaccard" if row == 0 else "Pearson R")
            if row == 2:
                ax.set_xlabel("Minimum sequence separation |i − j|")
        axes[row,0].text(0,1.22,label,transform=axes[row,0].transAxes,fontsize=12,fontweight="bold")
    fig.subplots_adjust(left=.10, right=.97, bottom=.12, top=.89, hspace=.67, wspace=.22)
    fig.suptitle("Contact agreement across sequence-separation cutoffs", fontsize=15, y=.99)
    fig.text(.5,.035,"GREMLIN uses item 2 offsets; MSA Transformer uses center-cropping.\n"
        "Pearson uses continuous scores on shared residue pairs, with no z-score floor.\n"
        "Boxes: IQR; line: median; whiskers: 1.5×IQR. Repeated proteins; descriptive distributions.",
        ha="center",fontsize=9)
    fig.savefig(out / "contact_agreement_combined.png",dpi=300)
    fig.savefig(out / "contact_agreement_combined.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=g.ROOT)
    parser.add_argument("--item2-dir", type=Path)
    parser.add_argument("--out-dir", type=Path, default=g.ROOT / "for_github_repo/analysis/pearson_topl")
    args = parser.parse_args()
    root, out = args.data_root.resolve(), args.out_dir.resolve()
    prior = (args.item2_dir or root / "for_github_repo/analysis/gremlin_separation").resolve()
    g.require(not out.exists(), f"Choose a new output directory: {out}")
    hashes = {}

    def track(path):
        path = Path(path).resolve()
        key = str(path.relative_to(root)) if path.is_relative_to(root) else str(path)
        if key not in hashes:
            g.require(path.is_file() and path.stat().st_size > 0, f"Missing/empty input: {path}")
            hashes[key] = dict(sha256=g.sha256(path), bytes=path.stat().st_size)
        return path

    prior_prov = json.loads(track(prior / "provenance.json").read_text())
    complete = json.loads(track(prior / "COMPLETE.json").read_text())
    g.require(complete["provenance_sha256"] == g.sha256(prior / "provenance.json"), "Item 2 provenance changed")
    g.require(prior_prov["output_sha256"]["comparisons.csv"] == g.sha256(track(prior / "comparisons.csv")), "Item 2 results changed")
    item2 = pd.read_csv(prior / "comparisons.csv")
    item2 = item2[item2.separation.isin(SEPARATIONS)]
    g.require(len(item2) == 2172 and not item2.duplicated(["pair", "separation", "reference"]).any(), "Incomplete item 2 offsets")
    offsets = item2.set_index(["pair", "separation", "reference"]).offset.to_dict()
    targets = pd.read_csv(track(root / "for_github_repo/msa-swaps-github/data/targets.csv"))
    classes = dict(zip(targets.pdb_id.str.lower(), targets.length_class))
    expected = {(a,b) for a in classes for b in classes if a != b and classes[a] == classes[b]}
    pairs = item2[["protein_a", "protein_b"]].drop_duplicates().sort_values(["protein_a", "protein_b"])
    g.require(len(classes)==34 and len(pairs)==362 and set(map(tuple,pairs.to_numpy()))==expected, "Incomplete pair manifest")
    audit = []

    def gremlin(score_path, alignment_path, identifier):
        track(score_path); track(alignment_path)
        # Reused inputs must be byte-identical to item 2, whose offsets are used.
        for path in (score_path, alignment_path):
            key = str(path.resolve().relative_to(root))
            g.require(hashes[key]["sha256"] == prior_prov["source_files"][key]["sha256"], f"Item 2 source changed: {path}")
        mask, n = g.gap_mask(alignment_path)
        frame = g.read_scores(score_path, mask)
        ij = frame[["i","j"]].to_numpy(dtype=int)
        retained = frame.loc[mask[ij[:,0]] & mask[ij[:,1]], ["i","j","zscore"]].copy()
        audit.append(dict(method="GREMLIN", identifier=identifier, columns=len(mask),
            retained_columns=int(mask.sum()), score_pairs=len(frame), usable_score_pairs=len(retained),
            n_sequences=n, max_asymmetry=np.nan))
        return retained, len(mask)

    def matrix(path, identifier):
        data = np.loadtxt(track(path),delimiter=",")
        g.require(data.ndim == 2 and data.shape[0] == data.shape[1] and np.isfinite(data).all(), f"Invalid matrix: {path}")
        asymmetry = float(np.max(abs(data-data.T)))
        g.require(asymmetry <= 1e-6, f"Asymmetric matrix: {path}")
        audit.append(dict(method="MSA Transformer", identifier=identifier, columns=len(data),
            retained_columns=len(data), score_pairs=len(data)*(len(data)-1)//2,
            usable_score_pairs=len(data)*(len(data)-1)//2, n_sequences=np.nan, max_asymmetry=asymmetry))
        return data

    native_g, native_m = {}, {}
    for protein in sorted(classes):
        native_g[protein] = gremlin(root/f"gremlin/results_gremlin_proper/full_{protein}.csv",
            root/f"original_proper_msas_af_sampling/{protein}/{protein}.fas",protein)
        native_m[protein] = matrix(root/f"original_proper_msas_af_sampling/{protein}/full_{protein}.csv",protein)
    # Stage all outputs in a new directory; COMPLETE is written only after validation.
    out.mkdir(parents=True)
    records = []
    member_count = 0
    membership_path = out / "selected_residue_pairs.csv.gz"
    print("Loaded 34 native GREMLIN and MSA Transformer inputs.",flush=True)
    with gzip.open(membership_path,"wt",newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method","pair","reference","separation","offset","i","j","swapped_score","native_score","top_swapped","top_native"])
        for count, (a,b) in enumerate(pairs.itertuples(index=False,name=None),1):
            pair=f"{a}_{b}_aligned"; fixed=b=="2mm4" and a in g.FIXED
            fas = root/(f"2mm4_fixes/{pair}/{pair}.fas" if fixed else f"original_junk_msas/{a}/{pair}/{pair}.fas")
            sg, sl = gremlin(root/f"gremlin/results_gremlin_junk/full_{pair}.csv",fas,pair)
            sm = matrix(root/(f"2mm4_fixes/{pair}/full_{pair}.csv" if fixed else f"original_junk_msas_af_sampling/{a}/{pair}/full_{pair}.csv"),pair)
            for separation in SEPARATIONS:
                for ref, protein in (("A",a),("B",b)):
                    ng,nl=native_g[protein]
                    off=int(offsets[(pair,separation,ref)])
                    for method in ("GREMLIN","MSA Transformer"):
                        if method == "GREMLIN":
                            ij,x,y,length=gremlin_pairs(sg,ng,sl,nl,off,separation)
                            offset=off; crop_a=crop_b=0
                        else:
                            ij,x,y,length,crop_a,crop_b=matrix_pairs(sm,native_m[protein],separation)
                            offset=0
                        result,selected=top_union(x,y,length)
                        if method=="GREMLIN" and off:
                            _,zx,zy,zl=gremlin_pairs(sg,ng,sl,nl,0,separation)
                            zero,_=top_union(zx,zy,zl)
                        else:
                            zero=result
                        records.append(dict(method=method,pair=pair,protein_a=a,protein_b=b,
                            reference=ref,length_class=classes[a],corrected=fixed,separation=separation,
                            offset=offset,crop_start_swapped=crop_a,crop_start_native=crop_b,
                            zero_offset_r=zero["pearson_r"],zero_offset_status=zero["status"],
                            zero_offset_numpy_error=zero["numpy_abs_error"],**result))
                        if selected is not None:
                            union,topx,topy=selected
                            for k in np.flatnonzero(union):
                                writer.writerow([method,pair,ref,separation,offset,int(ij[k,0]),int(ij[k,1]),
                                    float(x[k]),float(y[k]),int(topx[k]),int(topy[k])])
                                member_count+=1
            if count%60==0 or count==362:
                print(f"Processed {count}/362 pairs.",flush=True)
    frame=pd.DataFrame(records)
    g.require(len(frame)==4344 and not frame.duplicated(["method","pair","separation","reference"]).any(),"Incomplete Pearson output")
    g.require(int(frame.n_union.sum())==member_count,"Selected residue-pair export mismatch")
    summary=summary_table(frame)
    directions=[]
    for (method,sep),sub in frame.groupby(["method","separation"]):
        for score in ("pearson_r","zero_offset_r","all_eligible_r"):
            wide=sub.pivot(index="pair",columns="reference",values=score).dropna()
            directions.append(dict(method=method,separation=sep,score=score,n_both_defined=len(wide),
                n_B_gt_A=int((wide.B>wide.A).sum()),n_equal=int((wide.B==wide.A).sum()),
                n_A_gt_B=int((wide.A>wide.B).sum()),n_undefined=362-len(wide)))
    for name,table in (("comparisons.csv",frame),("summary.csv",summary),
            ("donor_query_direction.csv",pd.DataFrame(directions)),("input_audit.csv",pd.DataFrame(audit))):
        table.to_csv(out/name,index=False)
    figure(frame,item2,out)
    for path in (Path(__file__),Path(g.__file__),Path(__file__).with_name("test_topl_pearson.py")):
        track(path)
    errors=frame[["numpy_abs_error","all_eligible_numpy_error","zero_offset_numpy_error"]].to_numpy()
    provenance=dict(created_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),
        versions={"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__},
        data_root=str(root),item2_dir=str(prior),parameters=dict(separations=list(SEPARATIONS),
            top_selection="union of each input's top-L score-threshold set, retaining cutoff ties",
            score_universe="same finite residue pairs in both inputs, after separation and applicable gap filtering",
            gremlin_score="zscore, no minimum score floor",gremlin_gap_fraction="strictly <0.25 in both inputs",
            gremlin_L="length of shared original-coordinate interval after fixed item 2 offset",
            msatr_L="minimum matrix dimension after historical center-cropping",
            gremlin_offset="reuse item 2 Jaccard-selected offset for each cutoff; never optimize Pearson",
            missing_scores="exclude from shared universe; never fill with zero",
            insufficient_pairs="undefined primary R if fewer than L eligible pairs",
            constant_vectors="undefined R; never replace with zero",threads=1),
        validation=dict(expected_rows=4344,actual_rows=len(frame),selected_pair_rows=member_count,
            maximum_numpy_error=float(np.nanmax(errors)),status_counts=frame.status.value_counts().to_dict()),
        limitations=["Top-L union conditioning can induce negative correlations when different pairs score highly in each input.",
            "GREMLIN donor offsets were optimized using Jaccard, so correspondence is not independent of contact agreement.",
            "Shared native proteins and correlated residue pairs preclude treating all observations as independent replicates.",
            "Native gap-mask alignment historical identity remains unconfirmed; same files as item 2 are enforced."],
        source_files=hashes,output_sha256={p.name:g.sha256(p) for p in sorted(out.iterdir()) if p.is_file()})
    (out/"provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
    (out/"COMPLETE.json").write_text(json.dumps(dict(comparisons=len(frame),pairs=362,
        provenance_sha256=g.sha256(out/"provenance.json")),indent=2)+"\n")
    print(summary[summary.length_class=="all"].to_string(index=False))
    print(pd.DataFrame(directions).query("score=='pearson_r'").to_string(index=False))
    print("Statuses:",frame.status.value_counts().to_dict())
    print("Output:",out)


if __name__=="__main__":
    main()
