"""Small, dependency-free input/provenance helpers for the three-input pilot."""
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path

ALPHABET = set("ARNDCQEGHILKMFPSTWYVXBZUO-")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def atomic_bytes(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".partial-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_json(path, value):
    atomic_bytes(path, (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode())


def read_fasta(path):
    """Remove A3M insertion characters; never uppercase insertions or drop columns."""
    rows = []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            rows.append([line[1:], ""])
        else:
            require(bool(rows), f"Sequence before FASTA header: {path}")
            rows[-1][1] += re.sub(r"[a-z]|\.|\*", "", line)
    require(bool(rows), f"Empty FASTA: {path}")
    length = len(rows[0][1])
    require(length > 0, f"Empty query: {path}")
    require(all(len(s) == length for _, s in rows), f"Nonrectangular alignment: {path}")
    require(all(set(s) <= ALPHABET for _, s in rows), f"Unsupported residues: {path}")
    require("-" not in rows[0][1], f"Pilot requires an ungapped query: {path}")
    return rows


def fasta_bytes(rows):
    return "".join(f">{header}\n{seq}\n" for header, seq in rows).encode()


def infer_donor_offset(swapped, native):
    """Register columns from shared homolog rows, independently of contact scores."""
    reference = {}
    for header, seq in native[1:]:
        if header in reference:
            require(reference[header] == seq, "Conflicting duplicate donor header")
        reference[header] = seq
    shared = {h: (s, reference[h]) for h, s in swapped[1:] if h in reference}
    require(len(shared) >= 3, "Need at least three shared homolog identifiers")
    ls, ln = len(swapped[0][1]), len(native[0][1])
    candidates = []
    for offset in range(-12, 13):
        lo, hi = max(0, -offset), min(ls, ln - offset)
        if hi - lo < 0.8 * min(ls, ln):
            continue
        matches = total = 0
        for s, n in shared.values():
            for a, b in zip(s[lo:hi], n[lo + offset:hi + offset]):
                if a != "-" or b != "-":
                    total += 1
                    matches += a == b
        candidates.append((matches / total if total else 0, offset, lo, hi))
    candidates.sort(reverse=True)
    best = candidates[0]
    require(best[0] >= 0.99, "Donor registration has <99% identity")
    require(sum(abs(c[0] - best[0]) < 1e-12 for c in candidates) == 1,
            "Ambiguous donor registration")
    return dict(offset=best[1], start=best[2], stop=best[3], identity=best[0],
                shared_homologs=len(shared), convention="donor_column = swap_column + offset")
