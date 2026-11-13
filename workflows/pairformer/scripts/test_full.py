import tempfile
import unittest
from pathlib import Path
import numpy as np
from prepare_full import select_indices, mapping
from analyze_full import compare, replay
from run_full import save_matrix
from run_inference import validated_checkpoint


class FullTests(unittest.TestCase):
    def test_manuscript_depth_cap_and_nested_selection(self):
        rng=np.random.default_rng(8)
        rows=[[str(i), ''.join(x)] for i,x in enumerate(rng.choice(list('ACGT'),size=(600,8))) ]
        selected=select_indices(rows)
        self.assertEqual(len(selected),512)
        self.assertEqual(selected[0],0)
        self.assertTrue(set(select_indices(rows,256)).issubset(selected))

    def test_query_order_and_ties(self):
        rows=[['q','AAAA'],['a','AAAA'],['b','CCCC'],['c','CCAA'],['d','AACC']]
        # After opposite sequences are selected, remaining candidates tie.
        self.assertEqual(select_indices(rows,3).tolist(),[0,1,2])
        self.assertEqual(select_indices(rows,8).tolist(),list(range(5)))

    def test_ambiguous_mapping_remains_undefined(self):
        a=[['q','AAAAAAAAA']]+[[str(i),'AAAAAAAAA'] for i in range(3)]
        b=[['q','AAAAAAAAAAA']]+[[str(i),'AAAAAAAAAAA'] for i in range(3)]
        self.assertEqual(mapping(a,b)['status'],'undefined')

    def test_independent_replay_all_length_classes(self):
        rng=np.random.default_rng(42)
        for L in (56,100,150):
            a=rng.random((L,L)).astype('float32'); a=(a+a.T)/2
            b=rng.random((L,L)).astype('float32'); b=(b+b.T)/2
            for sep in (6,12,24):
                replay(a,b,L,sep,compare(a,b,L,sep))

    def test_checkpoint_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); save_matrix(p,np.eye(5,dtype='float32'),{'run_fingerprint':'test'})
            self.assertTrue(validated_checkpoint(p,'test',5))


if __name__=='__main__': unittest.main()
