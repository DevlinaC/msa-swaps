"""Analysis functions preserved from msa_transformer_jaccard_v3.py.

Dataset discovery and fixed cluster paths are replaced by compute_contacts.py.
See docs/file_sources.csv for original source checksums.
"""
import numpy as np


def crop_center(m: np.ndarray, n: int) -> np.ndarray:
    """Crop symmetric matrix m to n x n, removing extras from both termini.
    Odd differences take one extra residue off the front (N-terminus)."""
    d = m.shape[0] - n
    if d == 0:
        return m
    front = (d + 1) // 2  # ceil -> 1 extra off the front when d is odd
    back = d // 2         # floor
    return m[front:m.shape[0] - back, front:m.shape[0] - back]


def top(m: np.ndarray, mask: np.ndarray, k: int) -> np.ndarray:
    thr = np.sort(m[mask])[-k]
    return (m >= thr) & mask


def tolerant_jaccard(matrix_1: np.ndarray, matrix_2: np.ndarray) -> float:
    """Crop both matrices to shared size and compute the tolerant Jaccard
    index (top-L vs top-1.5L agreement, upper triangle, |i-j| >= 6)."""
    n = min(matrix_1.shape[0], matrix_2.shape[0])
    matrix_1 = crop_center(matrix_1, n)
    matrix_2 = crop_center(matrix_2, n)
    assert matrix_1.shape == matrix_2.shape == (n, n)

    mask = np.triu(np.ones((n, n), bool), k=6)
    L = n

    A, B = top(matrix_1, mask, L), top(matrix_2, mask, L)
    Awide = top(matrix_1, mask, int(1.5 * L))
    Bwide = top(matrix_2, mask, int(1.5 * L))

    # agreed = top-L in one and top-1.5L in the other
    agreed = (A & Bwide) | (B & Awide)
    return agreed.sum() / (A | B).sum()
