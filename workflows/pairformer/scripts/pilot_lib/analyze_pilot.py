#!/usr/bin/env python3
"""Compare Pairformer contact maps on a common, homolog-registered pair universe."""
import argparse
import csv
import io
import json
from pathlib import Path
import numpy as np
from common import atomic_bytes, require, sha256, write_json
from run_inference import validated_checkpoint


def pearson(x, y):
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def compare(a, b, length, separation):
    i, j = np.triu_indices(length, k=separation)
    require(len(i) >= length, "Fewer than L eligible pairs")
    x, y = a[i, j], b[i, j]
    sx = x >= np.partition(x, -length)[-length]
    sy = y >= np.partition(y, -length)[-length]
    union = sx | sy
    return dict(separation=separation, L=length, eligible_pairs=len(i),
                top_swap=int(sx.sum()), top_reference=int(sy.sum()),
                ordinary_jaccard=float((sx & sy).sum() / union.sum()),
                pearson_union_top_L=pearson(x[union], y[union]),
                pearson_all_eligible=pearson(x, y), ties="retained")


def analyze(root):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    manifest = json.loads((root / "inputs/manifest.json").read_text())
    out = root / "results/all_rows_cb"
    run = json.loads((out / "RUN.json").read_text())
    matrices, evidence = {}, {}
    for spec in manifest["inputs"]:
        folder = out / spec["id"]
        require(validated_checkpoint(folder, run["run_fingerprint"], spec["expected_length"]),
                f"Missing validated result: {spec['id']}")
        matrices[spec["role"]] = np.load(folder / "contacts.npy", allow_pickle=False)
        evidence[spec["id"]] = sha256(folder / "contacts.npy")
    mapping = manifest["donor_mapping"]
    idx = np.arange(mapping["start"], mapping["stop"])
    donor_idx = idx + mapping["offset"]
    maps = [matrices["query"][np.ix_(idx, idx)], matrices["swap"][np.ix_(idx, idx)],
            matrices["donor"][np.ix_(donor_idx, donor_idx)]]
    results = []
    for sep in (6, 12, 24):
        for label, reference in (("proper_1gb1", maps[0]), ("proper_2mm4", maps[2])):
            results.append(dict(reference=label, **compare(maps[1], reference, len(idx), sep)))
    dest = out / "analysis"
    dest.mkdir(exist_ok=True)
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(results[0]))
    writer.writeheader()
    writer.writerows(results)
    atomic_bytes(dest / "comparisons.csv", stream.getvalue().encode())
    figure, axes = plt.subplots(1, 3, figsize=(11, 3.6), constrained_layout=True)
    for axis, data, title in zip(axes, maps, ["Proper 1gb1", "1gb1 + 2mm4 MSA", "Proper 2mm4 (mapped)"]):
        display = data.copy()
        np.fill_diagonal(display, np.nan)
        im = axis.imshow(display, vmin=0, vmax=1, origin="lower", cmap="viridis")
        axis.set(title=title, xlabel="Mapped alignment column (0-based)", ylabel="Mapped alignment column (0-based)")
    figure.colorbar(im, ax=axes, label="Predicted Cβ contact probability", shrink=0.8)
    for suffix in ("png", "pdf"):
        buffer = io.BytesIO()
        figure.savefig(buffer, format=suffix, dpi=180)
        atomic_bytes(dest / ("contact_maps." + suffix), buffer.getvalue())
    plt.close(figure)
    write_json(dest / "COMPLETE.json", dict(run_fingerprint=run["run_fingerprint"],
               input_contact_hashes=evidence, donor_mapping=mapping, comparisons=results,
               analysis_code_sha256=sha256(Path(__file__)),
               interpretation="Technical pilot; not a structural-accuracy benchmark or an exact historical MSA Transformer replication."))
    print(stream.getvalue())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    analyze(parser.parse_args().root.resolve())
