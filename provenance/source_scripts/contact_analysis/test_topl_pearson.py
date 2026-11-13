"""Boundary checks for common-coordinate score correlation and selection."""
import unittest
import numpy as np
import pandas as pd
import analyze_topl_pearson as p


class PearsonTests(unittest.TestCase):
    def test_union_retains_pair_identity_and_ties(self):
        x=np.array([9.,8.,8.,1.]); y=np.array([1.,8.,8.,9.])
        result,selected=p.top_union(x,y,2)
        self.assertEqual(result["n_top_swapped"],3)
        self.assertEqual(result["n_top_native"],3)
        self.assertEqual(result["n_union"],4)
        self.assertLess(result["pearson_r"],0)
        np.testing.assert_array_equal(selected[0],[True]*4)

    def test_undefined_and_positive_affine_scaling(self):
        self.assertEqual(p.top_union(np.array([2.]),np.array([3.]),2)[0]["status"],"insufficient_eligible_for_L")
        self.assertEqual(p.correlation([1,1,1],[1,2,3])[1],"constant_vector")
        self.assertEqual(p.correlation([1],[2])[1],"too_few_pairs")
        r,status,_=p.correlation([1,3,7],[4,8,16])
        self.assertEqual(status,"ok");self.assertAlmostEqual(r,1)

    def test_offset_join_excludes_unavailable_pairs(self):
        a=pd.DataFrame({"i":[0,1,2],"j":[6,7,8],"zscore":[4.,3.,2.]})
        b=pd.DataFrame({"i":[1,3,4],"j":[7,9,10],"zscore":[8.,4.,1.]})
        ij,x,y,length=p.gremlin_pairs(a,b,10,12,1,6)
        np.testing.assert_array_equal(ij,[[1,7],[3,9]])
        np.testing.assert_array_equal(x,[4,2]);np.testing.assert_array_equal(y,[8,4])
        self.assertEqual(length,10)
        self.assertEqual(len(p.gremlin_pairs(a,b,10,12,1,7)[0]),0)

    def test_odd_center_crop_and_inclusive_separation(self):
        a=np.arange(81.).reshape(9,9);b=np.arange(64.).reshape(8,8)
        ij,x,y,length,ca,cb=p.matrix_pairs(a,b,6)
        self.assertEqual((length,ca,cb),(8,1,0))
        np.testing.assert_array_equal(ij,[[0,6],[0,7],[1,7]])
        np.testing.assert_array_equal(x,[a[1,7],a[1,8],a[2,8]])


if __name__=="__main__":
    unittest.main()
