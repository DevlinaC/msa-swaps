# Combines the 2-seq MSA (target + 1 extra) with random draws from the swapped MSA

import os
import glob
import sys
import numpy as np


def read_a3m(path):
    """Return list of (header, seq) preserving a3m case/insertions.
    Skips '#' comment lines (HHblits header) and blank lines."""
    records = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def combine_msa(two_seq_msa, swapped_msa, out_dir,
                total=512, seed=0):
    """
    Final MSA = target (row 0 of 2seq) + extra (row 1 of 2seq)
                + random sample from swapped_msa, up to `total` rows.
    `total` is INCLUSIVE of the target (AlphaFold convention).
    """
    two = read_a3m(two_seq_msa)
    if len(two) < 2:
        print(f"  skip: {two_seq_msa} has <2 sequences")
        return None

    target = two[0]          # (header, seq)
    extra = two[1]

    pool = read_a3m(swapped_msa)

    # Drop anything from the pool that duplicates target or extra
    # (swapped MSA usually leads with the query too).
    seen = {target[1], extra[1]}
    pool = [rec for rec in pool if rec[1] not in seen]

    n_sample = total - 2                      # slots left after target + extra
    rng = np.random.default_rng(seed)
    if len(pool) <= n_sample:
        sampled = pool                        # not enough; take all
    else:
        idx = rng.choice(len(pool), size=n_sample, replace=False)
        sampled = [pool[i] for i in idx]

    combined = [target, extra] + sampled

    os.makedirs(out_dir, exist_ok=True)
    base_two = os.path.basename(two_seq_msa).replace(".a3m", "")
    base_swap = os.path.basename(swapped_msa).replace(".a3m", "")
    out_path = os.path.join(out_dir, f"{base_swap}_with_2seq.a3m")

    with open(out_path, "w") as f:
        for hdr, seq in combined:
            f.write(f">{hdr}\n{seq}\n")

    print(f"  wrote {out_path}: {len(combined)} rows "
          f"(target + extra + {len(sampled)} sampled)")
    return out_path


def main():
    total = 512          # inclusive of target
    out_dir = './original_list/2seq_msas/2seq_with_swap/'

    two_seq_dir = './original_list/2seq_msas/2seq/'
    two_seq_msas = glob.glob(f'{two_seq_dir}*.a3m')
    for two_seq_msa in two_seq_msas:
        pdb = two_seq_msa.split('/')[-1].split('_')[0]
        print(two_seq_msa, '->', pdb)
        swapped_msas = glob.glob(
            f'./original_list/original_junk_msas_af_sampling/{pdb}/*/*.a3m')
        for swapped_msa in swapped_msas:
            combine_msa(two_seq_msa, swapped_msa, f'{out_dir}{pdb}/',
                        total=total, seed=0)
            #sys.exit()


if __name__ == "__main__":
    main()
