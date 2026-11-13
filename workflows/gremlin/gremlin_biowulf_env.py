"""
Python 3 executable. Newer versions of Python may cause errors w/ the libraries.
"""
#!/usr/bin/env python3
"""
Required environment variables:
    GREMLIN_FASTA   path to one input .fas alignment
    GREMLIN_OUTDIR  directory where outputs for that input should be written

Common optional environment variables (all are set to defaults of the original GREMLIN script if not set):
    GREMLIN_SAMPLE_ID       output prefix; default is the FASTA filename stem
    GREMLIN_OPT_ITER        Adam iterations; default 100
    GREMLIN_OPT_RATE        Adam learning rate; default 1.0 
    GREMLIN_GAP_CUTOFF      remove columns with gap fraction >= this value; default 0.5 
    GREMLIN_EFF_CUTOFF      sequence identity cutoff for sequence weights; default 0.8 
    GREMLIN_BATCH_SIZE      mini-batch size; default 0 means use all sequences 
    GREMLIN_MIN_SEPARATION  minimum j-i separation for top contact table; default 5
    GREMLIN_SEED            random seed; default 0
    GREMLIN_DRY_RUN         1 means parse/preprocess only & skip TensorFlow; default 0
    GREMLIN_SAVE_MRF        1 means also save fitted V/W tensors as .npz; default 0
    GREMLIN_ALLOW_GPU       1 means allow TensorFlow to see GPUs; default 0
"""

# Import standard libraries
from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

# GREMLIN codes
ALPHABET = "ARNDCQEGHILKMFPSTWYV-"
STATES = len(ALPHABET)
AA_TO_NUM = {aa: index for index, aa in enumerate(ALPHABET)}


@dataclass
class Settings:
    """Container for all run settings read from environment variables."""

    fasta: Path
    outdir: Path
    sample_id: str
    limit: int
    gap_cutoff: float
    eff_cutoff: float
    opt_iter: int
    opt_rate: float
    batch_size: Optional[int]
    min_separation: int
    seed: int
    dry_run: bool
    save_mrf: bool
    allow_gpu: bool


def fail(message: str, exit_code: int = 2) -> None:
    """Print a clear error message to stderr & exit w/ a non-zero code."""
    print(f"ERROR: {message}", file=sys.stderr, flush=True)
    raise SystemExit(exit_code)


def env_text(name: str, default: Optional[str] = None, required: bool = False) -> str:
    """Read a string environment variable, optionally requiring a non-empty value."""
    value = os.environ.get(name, default)
    if value is None or str(value).strip() == "":
        if required:
            fail(f"required environment variable {name} is not set")
        return ""
    return str(value).strip()


def env_int(name: str, default: int, minimum: Optional[int] = None) -> int:
    """Read an integer environment variable & validate its minimum value."""
    text = env_text(name, str(default))
    try:
        value = int(text)
    except ValueError:
        fail(f"{name} must be an integer, but got {text!r}")
    if minimum is not None and value < minimum:
        fail(f"{name} must be >= {minimum}, but got {value}")
    return value


def env_float(name: str, default: float, minimum: Optional[float] = None, maximum: Optional[float] = None) -> float:
    """Read a float environment variable & validate an optional numeric range."""
    text = env_text(name, str(default))
    try:
        value = float(text)
    except ValueError:
        fail(f"{name} must be a number, but got {text!r}")
    if minimum is not None and value < minimum:
        fail(f"{name} must be >= {minimum}, but got {value}")
    if maximum is not None and value > maximum:
        fail(f"{name} must be <= {maximum}, but got {value}")
    return value


def env_bool(name: str, default: bool = False) -> bool:
    """Read a boolean environment variable from common true/false spellings."""
    default_text = "1" if default else "0"
    text = env_text(name, default_text).lower()
    if text in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "f", "no", "n", "off"}:
        return False
    fail(f"{name} must be boolean-like, for example 0/1 or true/false; got {text!r}")
    return False


def print_usage() -> None:
    """Print a minimal usage message using only environment variables."""
    print(
        "Usage:\n"
        "  GREMLIN_FASTA=/path/to/sample.fas \\\n"
        "  GREMLIN_OUTDIR=/path/to/output_dir \\\n"
        "  python -u gremlin_biowulf_env.py\n\n"
        "Required environment variables:\n"
        "  GREMLIN_FASTA                 input .fas alignment\n"
        "  GREMLIN_OUTDIR                output directory\n\n"
        "Optional environment variables:\n"
        "  GREMLIN_SAMPLE_ID             output prefix; default: FASTA stem\n"
        "  GREMLIN_OPT_ITER              default: 100\n"
        "  GREMLIN_OPT_RATE              default: 1.0\n"
        "  GREMLIN_GAP_CUTOFF            default: 0.5\n"
        "  GREMLIN_EFF_CUTOFF            default: 0.8\n"
        "  GREMLIN_BATCH_SIZE            default: 0, meaning all sequences\n"
        "  GREMLIN_MIN_SEPARATION        default: 5\n"
        "  GREMLIN_SEED                  default: 0\n"
        "  GREMLIN_DRY_RUN               default: 0\n"
        "  GREMLIN_SAVE_MRF              default: 0\n"
        "  GREMLIN_ALLOW_GPU             default: 0\n",
        flush=True,
    )


def read_settings() -> Settings:
    """Read all settings from environment variables & return a Settings object."""
    if env_bool("GREMLIN_HELP", False):
        print_usage()
        raise SystemExit(0)

    fasta = Path(env_text("GREMLIN_FASTA", required=True)).expanduser().resolve()
    outdir = Path(env_text("GREMLIN_OUTDIR", required=True)).expanduser().resolve()
    sample_id = env_text("GREMLIN_SAMPLE_ID", fasta.stem)
    limit = env_int("GREMLIN_LIMIT", -1, minimum=-1)
    gap_cutoff = env_float("GREMLIN_GAP_CUTOFF", 0.5, minimum=0.0, maximum=1.0)
    eff_cutoff = env_float("GREMLIN_EFF_CUTOFF", 0.8, minimum=0.0, maximum=1.0)
    opt_iter = env_int("GREMLIN_OPT_ITER", 100, minimum=1)
    opt_rate = env_float("GREMLIN_OPT_RATE", 1.0, minimum=0.0)
    raw_batch_size = env_int("GREMLIN_BATCH_SIZE", 0, minimum=0)
    batch_size = None if raw_batch_size == 0 else raw_batch_size
    min_separation = env_int("GREMLIN_MIN_SEPARATION", 5, minimum=0)
    seed = env_int("GREMLIN_SEED", 0)
    dry_run = env_bool("GREMLIN_DRY_RUN", False)
    save_mrf = env_bool("GREMLIN_SAVE_MRF", False)
    allow_gpu = env_bool("GREMLIN_ALLOW_GPU", False)

    if not fasta.is_file():
        fail(f"input FASTA does not exist or isn't a regular file: {fasta}")

    if not allow_gpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    return Settings(
        fasta=fasta,
        outdir=outdir,
        sample_id=sample_id,
        limit=limit,
        gap_cutoff=gap_cutoff,
        eff_cutoff=eff_cutoff,
        opt_iter=opt_iter,
        opt_rate=opt_rate,
        batch_size=batch_size,
        min_separation=min_separation,
        seed=seed,
        dry_run=dry_run,
        save_mrf=save_mrf,
        allow_gpu=allow_gpu,
    )


def load_tensorflow():
    """Import TensorFlow only after settings have been validated."""
    import tensorflow.compat.v1 as tf

    tf.disable_eager_execution()
    return tf


def aa2num(aa: str) -> int:
    """Convert one amino-acid character into the numeric state used by GREMLIN."""
    return AA_TO_NUM.get(aa.upper(), AA_TO_NUM["-"])


def parse_fasta(filename: Path, limit: int = -1):
    """Read a FASTA alignment & return arrays of headers & aligned sequences."""
    headers = []
    sequence_parts = []

    with filename.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if limit > 0 and len(headers) >= limit:
                    break
                headers.append(line[1:])
                sequence_parts.append([])
            else:
                if not sequence_parts:
                    raise ValueError(f"{filename}: sequence text appears before the 1st FASTA header at line {line_number}")
                sequence_parts[-1].append(line.upper())

    if not headers:
        raise ValueError(f"{filename}: no FASTA records were found")

    sequences = ["".join(parts) for parts in sequence_parts]
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError(f"{filename}: FASTA records aren't aligned; sequence lengths are {sorted(lengths)}")

    return np.array(headers), np.array(sequences)


def filt_gaps(msa: np.ndarray, gap_cutoff: float = 0.5):
    """Remove alignment columns whose gap fraction is >= to gap_cutoff."""
    gap_mask = (msa == STATES - 1).astype(float)
    gap_fraction = np.sum(gap_mask, axis=0) / float(msa.shape[0])
    non_gaps = np.where(gap_fraction < gap_cutoff)[0]
    if non_gaps.size < 2:
        raise ValueError("fewer than 2 alignment columns remain after gap filtering")
    return msa[:, non_gaps], non_gaps


def get_eff(msa: np.ndarray, eff_cutoff: float = 0.8) -> np.ndarray:
    """Compute sequence weights from pairwise identity among aligned sequences."""
    if msa.shape[0] == 1:
        return np.ones(1, dtype=np.float32)
    sequence_identity = 1.0 - squareform(pdist(msa, "hamming"))
    similar = (sequence_identity >= eff_cutoff).astype(float)
    weights = 1.0 / np.sum(similar, axis=-1)
    return weights.astype(np.float32)


def mk_msa(seqs: np.ndarray, gap_cutoff: float = 0.5, eff_cutoff: float = 0.8) -> dict:
    """Convert aligned sequence strings into the MSA dictionary expected by GREMLIN."""
    msa_ori = np.array([[aa2num(aa) for aa in sequence] for sequence in seqs], dtype=np.int32)
    msa, v_idx = filt_gaps(msa_ori, gap_cutoff=gap_cutoff)
    msa_weights = get_eff(msa, eff_cutoff=eff_cutoff)
    ncol = int(msa.shape[1])
    upper_pairs = np.stack(np.triu_indices(ncol, 1), axis=-1)
    w_idx = v_idx[upper_pairs]

    return {
        "msa_ori": msa_ori,
        "msa": msa,
        "weights": msa_weights,
        "neff": float(np.sum(msa_weights)),
        "v_idx": v_idx,
        "w_idx": w_idx,
        "nrow": int(msa.shape[0]),
        "ncol": ncol,
        "ncol_ori": int(msa_ori.shape[1]),
    }


def sym_w(tf, w):
    """Symmetrize a TensorFlow W tensor w/ shape (L, A, L, A)."""
    ncol = int(w.shape.as_list()[0])
    off_diagonal = np.reshape(1.0 - np.eye(ncol, dtype=np.float32), (ncol, 1, ncol, 1))
    off_diagonal = tf.constant(off_diagonal, dtype=tf.float32)
    w = w * off_diagonal
    w = w + tf.transpose(w, [2, 3, 0, 1])
    return w


def opt_adam(tf, loss, name: str, var_list=None, lr: float = 1.0, b1: float = 0.9, b2: float = 0.999, b_fix: bool = False):
    """Adam optimizer variant matching the original GREMLIN script structure."""
    if var_list is None:
        var_list = tf.trainable_variables()
    gradients = tf.gradients(loss, var_list)
    if b_fix:
        t = tf.Variable(0.0, name="t")

    opt_ops = []
    for number, (variable, gradient) in enumerate(zip(var_list, gradients)):
        if gradient is None:
            continue
        variable_shape = variable.shape.as_list()
        init = dict(initializer=tf.zeros_initializer(), trainable=False)
        mt = tf.get_variable(name + "_mt_" + str(number), shape=variable_shape, **init)
        vt = tf.get_variable(name + "_vt_" + str(number), shape=[], **init)
        mt_tmp = b1 * mt + (1.0 - b1) * gradient
        vt_tmp = b2 * vt + (1.0 - b2) * tf.reduce_sum(tf.square(gradient))
        lr_tmp = lr / (tf.sqrt(vt_tmp) + 1e-8)
        if b_fix:
            lr_tmp = lr_tmp * tf.sqrt(1.0 - tf.pow(b2, t)) / (1.0 - tf.pow(b1, t))
        opt_ops.append(variable.assign_add(-lr_tmp * mt_tmp))
        opt_ops.append(vt.assign(vt_tmp))
        opt_ops.append(mt.assign(mt_tmp))

    if b_fix:
        opt_ops.append(t.assign_add(1.0))
    return tf.group(opt_ops)


def GREMLIN(tf, msa: dict, opt_iter: int = 100, opt_rate: float = 1.0, batch_size: Optional[int] = None, seed: int = 0) -> dict:
    """Fit the GREMLIN Markov random field & return fitted V & W parameters."""
    tf.reset_default_graph()
    tf.set_random_seed(seed)
    np.random.seed(seed)

    ncol = msa["ncol"]
    MSA = tf.placeholder(tf.int32, shape=(None, ncol), name="msa")
    OH_MSA = tf.one_hot(MSA, STATES)
    MSA_weights = tf.placeholder(tf.float32, shape=(None,), name="msa_weights")

    V = tf.get_variable(name="V", shape=[ncol, STATES], initializer=tf.zeros_initializer())
    W = tf.get_variable(name="W", shape=[ncol, STATES, ncol, STATES], initializer=tf.zeros_initializer())
    W = sym_w(tf, W)

    def L2(x):
        return tf.reduce_sum(tf.square(x))

    VW = V + tf.tensordot(OH_MSA, W, 2)
    H = tf.reduce_sum(tf.multiply(OH_MSA, VW), axis=(1, 2))
    Z = tf.reduce_sum(tf.reduce_logsumexp(VW, axis=2), axis=1)
    PLL = H - Z

    L2_V = 0.01 * L2(V)
    L2_W = 0.01 * L2(W) * 0.5 * (ncol - 1) * (STATES - 1)
    loss = -tf.reduce_sum(PLL * MSA_weights) / tf.reduce_sum(MSA_weights)
    loss = loss + (L2_V + L2_W) / msa["neff"]
    opt = opt_adam(tf, loss, "adam", lr=opt_rate)

    def feed(feed_all: bool = False) -> dict:
        if batch_size is None or feed_all:
            return {MSA: msa["msa"], MSA_weights: msa["weights"]}
        idx = np.random.randint(0, msa["nrow"], size=batch_size)
        return {MSA: msa["msa"][idx], MSA_weights: msa["weights"][idx]}

    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("OMP_NUM_THREADS", "1")))
    session_config = tf.ConfigProto(intra_op_parallelism_threads=threads, inter_op_parallelism_threads=1)

    with tf.Session(config=session_config) as sess:
        sess.run(tf.global_variables_initializer())

        msa_cat = np.eye(STATES, dtype=np.float32)[msa["msa"]]
        pseudo_count = max(1e-6, 0.01 * np.log(max(msa["neff"], 1.000001)))
        weighted_counts = np.tensordot(msa["weights"], msa_cat, axes=(0, 0))
        V_ini = np.log(weighted_counts + pseudo_count).astype(np.float32)
        V_ini = V_ini - np.mean(V_ini, axis=-1, keepdims=True)
        sess.run(V.assign(V_ini))

        def get_loss() -> float:
            return round(float(sess.run(loss, feed(feed_all=True)) * msa["neff"]), 2)

        print("starting", get_loss(), flush=True)
        print_every = max(1, opt_iter // 10)
        for iteration in range(opt_iter):
            sess.run(opt, feed())
            if (iteration + 1) % print_every == 0 or (iteration + 1) == opt_iter:
                print("iter", iteration + 1, get_loss(), flush=True)

        V_out = sess.run(V)
        W_out = sess.run(W)

    tri = np.triu_indices(ncol, 1)
    W_out = W_out[tri[0], :, tri[1], :]

    return {
        "v": V_out,
        "w": W_out,
        "v_idx": msa["v_idx"],
        "w_idx": msa["w_idx"],
    }


def normalize(x: np.ndarray) -> np.ndarray:
    """Apply Box-Cox normalization & return z-scores w/ safeguards."""
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return x
    shifted = x - np.amin(x) + 1.0
    if np.allclose(shifted, shifted[0]):
        return np.zeros_like(shifted)
    y = stats.boxcox(shifted)[0]
    y_std = np.std(y)
    if y_std == 0.0:
        return np.zeros_like(y)
    return (y - np.mean(y)) / y_std


def get_mtx(mrf: dict) -> dict:
    """Convert fitted MRF weights into raw, APC-corrected, & z-scored pair scores."""
    raw = np.sqrt(np.sum(np.square(mrf["w"][:, :-1, :-1]), axis=(1, 2)))
    raw_sq = squareform(raw)
    denominator = np.sum(raw_sq)
    if denominator == 0.0:
        apc = np.zeros_like(raw)
    else:
        ap_sq = np.sum(raw_sq, axis=0, keepdims=True) * np.sum(raw_sq, axis=1, keepdims=True) / denominator
        apc = squareform(raw_sq - ap_sq, checks=False)
    return {
        "i": mrf["w_idx"][:, 0],
        "j": mrf["w_idx"][:, 1],
        "raw": raw,
        "apc": apc,
        "zscore": normalize(apc),
    }


def plot_mtx(mtx: dict, outfile: Path, key: str = "zscore", vmin: float = 1.0, vmax: float = 3.0, title: str = "") -> None:
    """Save a contact-score heatmap as a PNG file."""
    fig = plt.figure(figsize=(5, 5))
    plt.imshow(squareform(mtx[key]), cmap="Blues", interpolation="none", vmin=vmin, vmax=vmax)
    plt.grid(False)
    plt.xlabel("Sequence position")
    plt.ylabel("Sequence position")
    if title:
        plt.title(title)
    plt.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)


def contacts_to_dataframe(mtx: dict, msa: dict) -> pd.DataFrame:
    """Build the contact-score dataframe w/ zero-based & one-based labels."""
    table_data = dict(mtx)
    table_data["i_aa"] = np.array([ALPHABET[msa["msa_ori"][0][i]] + "_" + str(i + 1) for i in table_data["i"]])
    table_data["j_aa"] = np.array([ALPHABET[msa["msa_ori"][0][j]] + "_" + str(j + 1) for j in table_data["j"]])
    return pd.DataFrame(table_data, columns=["i", "j", "raw", "apc", "zscore", "i_aa", "j_aa"])


def write_summary(outfile: Path, settings: Settings, names: np.ndarray, msa: dict, full_csv: Path, top_csv: Path, heatmap_png: Path) -> None:
    """Write run settings & output paths to a tab-separated summary text file."""
    with outfile.open("w", encoding="utf-8") as handle:
        handle.write(f"sample_id\t{settings.sample_id}\n")
        handle.write(f"fasta\t{settings.fasta}\n")
        handle.write(f"records\t{len(names)}\n")
        handle.write(f"alignment_columns_original\t{msa['ncol_ori']}\n")
        handle.write(f"alignment_columns_after_gap_filter\t{msa['ncol']}\n")
        handle.write(f"effective_sequences_neff\t{msa['neff']:.6f}\n")
        handle.write(f"opt_iter\t{settings.opt_iter}\n")
        handle.write(f"opt_rate\t{settings.opt_rate}\n")
        handle.write(f"batch_size\t{0 if settings.batch_size is None else settings.batch_size}\n")
        handle.write(f"gap_cutoff\t{settings.gap_cutoff}\n")
        handle.write(f"eff_cutoff\t{settings.eff_cutoff}\n")
        handle.write(f"min_separation\t{settings.min_separation}\n")
        handle.write(f"seed\t{settings.seed}\n")
        handle.write(f"full_csv\t{full_csv}\n")
        handle.write(f"top_csv\t{top_csv}\n")
        handle.write(f"heatmap_png\t{heatmap_png}\n")


def run_one(settings: Settings) -> int:
    """Run the complete one-FASTA workflow & return a shell-style exit code."""
    settings.outdir.mkdir(parents=True, exist_ok=True)

    print(f"sample_id: {settings.sample_id}", flush=True)
    print(f"fasta: {settings.fasta}", flush=True)
    print(f"outdir: {settings.outdir}", flush=True)

    names, seqs = parse_fasta(settings.fasta, limit=settings.limit)
    msa = mk_msa(seqs, gap_cutoff=settings.gap_cutoff, eff_cutoff=settings.eff_cutoff)

    print(f"records: {len(names)}", flush=True)
    print(f"columns original/filtered: {msa['ncol_ori']}/{msa['ncol']}", flush=True)
    print(f"neff: {msa['neff']:.6f}", flush=True)

    if settings.dry_run:
        print("dry_run: parsed FASTA & built MSA successfully; TensorFlow wasn't imported.", flush=True)
        return 0

    tf = load_tensorflow()
    mrf = GREMLIN(tf, msa, opt_iter=settings.opt_iter, opt_rate=settings.opt_rate, batch_size=settings.batch_size, seed=settings.seed)
    mtx = get_mtx(mrf)
    contact_df = contacts_to_dataframe(mtx, msa)
    top_df = contact_df.loc[contact_df["j"] - contact_df["i"] > settings.min_separation].sort_values("zscore", ascending=False)

    full_csv = settings.outdir / f"full_{settings.sample_id}.csv"
    top_csv = settings.outdir / f"top_{settings.sample_id}.csv"
    heatmap_png = settings.outdir / f"{settings.sample_id}_heatmap.png"
    summary_txt = settings.outdir / f"{settings.sample_id}_summary.txt"

    contact_df.to_csv(full_csv, index=False)
    top_df.to_csv(top_csv, index=False)
    plot_mtx(mtx, heatmap_png, title=settings.sample_id)

    if settings.save_mrf:
        np.savez_compressed(
            settings.outdir / f"{settings.sample_id}_mrf.npz",
            v=mrf["v"],
            w=mrf["w"],
            v_idx=mrf["v_idx"],
            w_idx=mrf["w_idx"],
        )

    write_summary(summary_txt, settings, names, msa, full_csv, top_csv, heatmap_png)

    print(f"wrote: {full_csv}", flush=True)
    print(f"wrote: {top_csv}", flush=True)
    print(f"wrote: {heatmap_png}", flush=True)
    print(f"wrote: {summary_txt}", flush=True)
    return 0


def main() -> int:
    """Read settings from the environment & run 1 FASTA file."""
    settings = read_settings()
    return run_one(settings)


if __name__ == "__main__":
    raise SystemExit(main())
