"""Compare all 34 experimental targets with tmtools, including across size groups.

Run: python create_tmscore_matrix.py
Smoke test: python create_tmscore_matrix.py --targets 1BOR 1GB1 1DE3 --output-dir /tmp/tm_smoke

Select the first coordinate model and the mapped chain (A for this analysis).
Use only amino-acid residues with C-alpha coordinates, keeping modified amino
acids explicitly identified by PDB MODRES records. Missing coordinates are not fabricated. Observed
sequences must be subsequences of the target FASTAs; incompatible structures
fail validation rather than silently aligning an unrelated chain/full protein.

tm_align takes two (N, 3) arrays and their sequences, not file paths. Each pair
is aligned once in deterministic target order and returns two TM-scores,
normalized by its two observed coordinate lengths. These are NOT normalized by
FASTA length, aligned length, or the average length. The directed matrix's
[row i, column j] entry is normalized by row i's coordinate length. Its diagonal
is NA. The symmetric mean and max matrices combine the two returned scores;
these combinations are explicit summaries, not new TM-align normalizations.
They do not imply a second reverse-order alignment was run.

Statistics include each unordered, non-self pair exactly once. All, within-class,
and between-class summaries are reported separately for mean and max scores.
No fold-distinctness conclusion or universal cutoff is inferred automatically.
The pair table retains both raw scores. Older tmtools releases lack RMSD; in
that case it is NA, not reconstructed from an unrelated transformation metric.

Outputs include labeled CSV matrices, pairs.csv, targets.csv (coordinate
coverage), summary.json and provenance.json (hashes, versions and conventions).
Existing results in tmscores_data are never used as computational inputs.
Execution is sequential, with numerical-library thread counts bounded to one.
On Biowulf run inside sinteractive/sbatch, not on the login node.
"""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import socket
import sys
import tempfile
from datetime import datetime, timezone
from itertools import combinations

# This script is deliberately serial. Set these before importing NumPy/BioPython.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import PDBParser
from tmtools import tm_align
from create_table_S1 import BASE_DIR, EXPECTED_COUNTS, read_csv, read_fasta, unique_lookup


def checksum(path):
    """Hash local inputs for reproducibility; no network access is required."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_coordinates(path, chain_id, target_sequence):
    """Read one chain once, rejecting unknown residues and incompatible sequences."""
    model = next(PDBParser(PERMISSIVE=False, QUIET=True).get_structure(path.stem, str(path)).get_models(), None)
    if model is None or chain_id not in model:
        raise ValueError(f"{path}: chain {chain_id!r} absent from first model")
    modified = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith('MODRES') and line[16] == chain_id:
                modified[(int(line[18:22]), line[22], line[12:15].strip())] = line[24:27].strip()
    coords, sequence, residue_ids = [], [], []
    no_ca = 0
    for residue in model[chain_id]:
        parent = modified.get((residue.id[1], residue.id[2], residue.resname))
        if residue.id[0] != " " and parent is None:
            continue  # Exclude water, ions and non-protein ligands (including calcium).
        if "CA" not in residue:
            no_ca += 1
            continue
        amino_acid = protein_letters_3to1_extended.get((parent or residue.resname).upper())
        if amino_acid is None:
            raise ValueError(f"{path}: unsupported amino acid {residue.resname} at {residue.id}")
        # Biopython selects the highest-occupancy alternate location by default.
        coords.append(residue["CA"].coord)
        sequence.append(amino_acid)
        residue_ids.append(f"{residue.id[1]}{residue.id[2].strip()}")
    observed = "".join(sequence)
    if len(observed) < 3:
        raise ValueError(f"{path}: at least three C-alpha coordinates are required")
    iterator = iter(target_sequence)
    if not all(amino_acid in iterator for amino_acid in observed):
        raise ValueError(f"{path}: observed chain {chain_id} sequence is not a subsequence of target FASTA")
    array = np.ascontiguousarray(coords, dtype=np.float64)
    if array.shape != (len(observed), 3) or not np.isfinite(array).all():
        raise ValueError(f"{path}: invalid coordinate array")
    return array, observed, residue_ids, no_ca


def get_tmscore_matrix(targets):
    """Return row-normalized, mean, max matrices and one record per unique pair.

    targets are preloaded dictionaries containing pdb_id, chain_id, length_class,
    coords and observed_sequence. Matrices contain NaN on the diagonal. No
    structure is compared to itself, and no pair is silently skipped on failure.
    """
    n = len(targets)
    directed = np.full((n, n), np.nan)
    mean = directed.copy()
    maximum = directed.copy()
    pairs = []
    for i, j in combinations(range(n), 2):
        first, second = targets[i], targets[j]
        result = tm_align(first["coords"], second["coords"],
                          first["observed_sequence"], second["observed_sequence"])
        a, b = float(result.tm_norm_chain1), float(result.tm_norm_chain2)
        if any(not math.isfinite(score) or not 0 <= score <= 1 for score in (a, b)):
            raise ValueError(f"Invalid TM-score for {first['pdb_id']} / {second['pdb_id']}")
        directed[i, j], directed[j, i] = a, b
        mean[i, j] = mean[j, i] = (a + b) / 2
        maximum[i, j] = maximum[j, i] = max(a, b)
        rmsd = getattr(result, "rmsd", None)
        pairs.append({
            "pdb_id_1": first["pdb_id"], "chain_1": first["chain_id"],
            "pdb_id_2": second["pdb_id"], "chain_2": second["chain_id"],
            "length_class_1": first["length_class"], "length_class_2": second["length_class"],
            "n_ca_1": len(first["coords"]), "n_ca_2": len(second["coords"]),
            "tm_norm_1": a, "tm_norm_2": b, "tm_mean": (a + b) / 2, "tm_max": max(a, b),
            "rmsd_angstrom": float(rmsd) if rmsd is not None else "NA",
        })
    return directed, mean, maximum, pairs


def summarize(pairs):
    """Summarize unordered pairs, keeping normalization choices explicit."""
    groups = {"all": pairs}
    for pair in pairs:
        sizes = sorted((pair["length_class_1"], pair["length_class_2"]))
        label = f"{sizes[0]}_vs_{sizes[1]}"
        groups.setdefault(label, []).append(pair)
    result = {}
    for label, records in groups.items():
        result[label] = {"n_pairs": len(records)}
        for metric in ("tm_mean", "tm_max"):
            values = np.array([p[metric] for p in records])
            top = records[int(np.argmax(values))]
            result[label][metric] = {
                "min": float(values.min()), "mean": float(values.mean()),
                "median": float(np.median(values)), "max": float(values.max()),
                "std": float(values.std(ddof=1)) if len(values) > 1 else None,
                "std_ddof": 1,
                "max_pair": [top["pdb_id_1"], top["pdb_id_2"]],
            }
    return result


def atomic_write(path, writer):
    """Replace a complete output file; never expose a partly written CSV/JSON."""
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=path.parent,
                                         prefix=f".{path.name}.", delete=False) as handle:
            temporary = Path(handle.name)
            writer(handle)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def save_csv(path, fields, rows):
    def write(handle):
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    atomic_write(path, write)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdb-dir", type=Path, default=BASE_DIR / "original_pdbs")
    parser.add_argument("--chain-map", type=Path, default=BASE_DIR / "target_chains.csv")
    for size in EXPECTED_COUNTS:
        parser.add_argument(f"--fasta-{size}", type=Path, default=BASE_DIR / f"targets_{size}.fa")
    parser.add_argument("--targets", nargs="+", help="Optional subset for smoke tests; at least two distinct IDs")
    parser.add_argument("--output-dir", type=Path, default=BASE_DIR / "outputs/tmscores")
    args = parser.parse_args(argv)
    try:
        if socket.gethostname().split('.')[0].startswith('biowulf') and not os.environ.get('SLURM_JOB_ID'):
            raise ValueError("Run on a Biowulf compute allocation (sinteractive/sbatch), not the login node")
        mappings = unique_lookup(read_csv(args.chain_map, ["pdb_id", "chain_id"]), args.chain_map)
        pdbs = {}
        for path in args.pdb_dir.iterdir():
            if path.is_file() and path.suffix.lower() == '.pdb':
                key = path.stem.upper()
                if key in pdbs:
                    raise ValueError(f"Duplicate PDB filenames for {key}")
                pdbs[key] = path
        inputs, records, seen = [args.chain_map, Path(__file__).resolve()], [], set()
        for size, count in EXPECTED_COUNTS.items():
            path = getattr(args, f"fasta_{size}")
            inputs.append(path)
            fasta = list(read_fasta(path))
            if len(fasta) != count:
                raise ValueError(f"{path}: expected {count} targets, found {len(fasta)}")
            for pdb_id, header_chain, sequence in fasta:
                if pdb_id in seen:
                    raise ValueError(f"Duplicate FASTA target {pdb_id}")
                seen.add(pdb_id)
                records.append((size, pdb_id, header_chain, sequence))
        requested = [pid.upper() for pid in args.targets] if args.targets else list(seen)
        if len(requested) < 2 or len(set(requested)) != len(requested) or set(requested) - seen:
            raise ValueError("--targets must contain at least two distinct, known FASTA IDs")
        targets, coverage = [], []
        # Match Table S1 ordering, independently of the order of --targets.
        for size, pdb_id, header_chain, sequence in sorted(records, key=lambda r: (r[0], len(r[3]), r[1])):
            if pdb_id not in requested:
                continue
            chain = mappings.get(pdb_id, {}).get('chain_id', '')
            if not chain or (header_chain and header_chain != chain):
                raise ValueError(f"{pdb_id}: missing or conflicting chain assignment; regenerate target_chains.csv")
            if pdb_id not in pdbs:
                raise ValueError(f"Missing PDB for {pdb_id}")
            coords, observed, residue_ids, no_ca = load_coordinates(pdbs[pdb_id], chain, sequence)
            inputs.append(pdbs[pdb_id])
            targets.append(dict(pdb_id=pdb_id, chain_id=chain, length_class=size,
                                coords=coords, observed_sequence=observed))
            coverage.append(dict(pdb_id=pdb_id, chain_id=chain, length_class=size,
                                 fasta_length=len(sequence), n_ca=len(coords),
                                 coordinate_coverage=len(coords)/len(sequence), residues_without_ca=no_ca,
                                 observed_sequence=observed, residue_ids=';'.join(residue_ids),
                                 pdb_file=str(pdbs[pdb_id].resolve())))
        print(f"Aligning {len(targets)} targets: {len(targets)*(len(targets)-1)//2} unique non-self pairs", flush=True)
        directed, mean, maximum, pairs = get_tmscore_matrix(targets)
        summary = summarize(pairs)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        ids = [t['pdb_id'] for t in targets]
        for name, matrix in [('row_normalized', directed), ('mean', mean), ('max', maximum)]:
            rows = [dict(pdb_id=pid, **{other: 'NA' if i == j else float(matrix[i,j])
                                      for j, other in enumerate(ids)}) for i, pid in enumerate(ids)]
            save_csv(args.output_dir / f'tmscore_{name}.csv', ['pdb_id', *ids], rows)
        save_csv(args.output_dir / 'pairs.csv', list(pairs[0]), pairs)
        save_csv(args.output_dir / 'targets.csv', list(coverage[0]), coverage)
        atomic_write(args.output_dir / 'summary.json', lambda f: json.dump(summary, f, indent=2, allow_nan=False))
        provenance = dict(
            created_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
            packages={name: importlib.metadata.version(name) for name in ['tmtools', 'biopython', 'numpy']},
            target_order=ids, n_pairs=len(pairs), model='first coordinate model',
            normalization='observed amino-acid C-alpha counts; both tmtools scores retained',
            diagonal='NA; self-alignments never computed or included in summaries',
            matrix_convention='row_normalized[i,j] normalized by row i; one alignment per unordered pair',
            symmetric_summaries='mean and max of the two returned normalization scores',
            command=sys.argv if argv is None else [str(Path(__file__)), *argv],
            slurm_job_id=os.environ.get('SLURM_JOB_ID'), threads=1,
            inputs={str(p.resolve()): checksum(p) for p in inputs},
        )
        atomic_write(args.output_dir / 'provenance.json', lambda f: json.dump(provenance, f, indent=2, allow_nan=False))
        print(json.dumps(summary['all'], indent=2))
        print(f"Saved results to {args.output_dir}")
    except (OSError, ValueError, csv.Error) as error:
        parser.exit(1, f"Error: {error}\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
