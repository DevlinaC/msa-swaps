"""Small boundary tests for the scientific selection and missing-data rules."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import analyze_gremlin_separation as analysis


class SeparationTests(unittest.TestCase):
    def test_inclusive_separation_and_score(self):
        frame = pd.DataFrame({"i": [0]*5, "j": [5, 6, 12, 24, 25],
                              "zscore": [2, 1, 2, 2, .999]})
        mask = np.ones(26, dtype=bool)
        for sep, expected in ((6, {(0,6),(0,12),(0,24)}),
                              (12, {(0,12),(0,24)}), (24, {(0,24)})):
            self.assertEqual(analysis.select_contacts(frame, mask, sep)[0], expected)

    def test_strict_gap_cutoff_and_multiline_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "alignment.fas"
            path.write_text(">a\nA-\nA\n>b\nAAA\n>c\nAAA\n>d\nA-A\n")
            mask, n = analysis.gap_mask(path)
            np.testing.assert_array_equal(mask, [True, False, True])
            self.assertEqual(n, 4)
            path.write_text(">a\nA-A\n>b\nAAA\n>c\nAAA\n>d\nAAA\n")
            self.assertFalse(analysis.gap_mask(path)[0][1])  # exactly 25% excluded
            path.write_text(">a\nAA\n>b\nAAA\n")
            with self.assertRaises(ValueError):
                analysis.gap_mask(path)

    def test_top100_is_applied_after_separation(self):
        frame = pd.DataFrame({"i": [0]*102, "j": range(1,103),
                              "zscore": np.arange(103,1,-1)})
        selected, eligible, _ = analysis.select_contacts(frame, np.ones(103,dtype=bool), 6)
        self.assertEqual(eligible, 97)
        self.assertEqual(selected, {(0,j) for j in range(6,103)})
        selected, eligible, _ = analysis.select_contacts(frame, np.ones(103,dtype=bool), 1)
        self.assertEqual(len(selected),100)
        self.assertNotIn((0,101),selected)

    def test_offset_and_empty_policies(self):
        a,b = {(0,6),(3,15)}, {(1,7),(4,16)}
        self.assertEqual(analysis.compare(a,b,False)["jaccard"],0)
        match = analysis.compare(a,b,True)
        self.assertEqual((match["jaccard"],match["offset"]),(1,1))
        empty = analysis.compare(set(),set(),True)
        self.assertTrue(np.isnan(empty["jaccard"]))
        self.assertEqual(empty["historical_jaccard"],0)
        self.assertEqual(empty["offset"],0)
        self.assertEqual(analysis.compare(a,set(),True)["jaccard"],0)

    def test_bad_indices_fail_instead_of_skipping_mask(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "scores.csv"
            path.write_text("i,j,zscore\n0,5,1\n")
            with self.assertRaises(ValueError):
                analysis.read_scores(path,np.ones(5,dtype=bool))


if __name__ == "__main__":
    unittest.main()
