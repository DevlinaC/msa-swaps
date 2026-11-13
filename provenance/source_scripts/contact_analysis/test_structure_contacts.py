import unittest
import numpy as np
import tempfile
from pathlib import Path
from analyze_structure_contacts import calculate, infer_offset, structure


class Contacts(unittest.TestCase):
    def test_glycine_CA_and_missing_CB(self):
        lines=[]
        for serial,(atom,res,index,x) in enumerate([('CA','ALA',1,20.),('CB','ALA',1,0.),('CA','GLY',2,8.),('CA','ALA',3,12.)],1):
            lines.append(f'ATOM  {serial:5d} {atom:^4s} {res:3s} A{index:4d}    {x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{90.:6.2f}           C\n')
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'test.pdb';path.write_text(''.join(lines)+'END\n')
            distances,details=structure(path,'AGA')
        self.assertEqual(distances[0,1],8.)
        self.assertTrue(np.isnan(distances[0,2]))
        self.assertEqual(details['missing_contact_atoms'],1)

    def test_hand_calculated_overlap_and_ties(self):
        d=np.full((8,8),20.)
        d[0,6]=8.; d[0,7]=7.; d[1,7]=9.
        ij=np.array([[0,6],[0,7],[1,7]])
        r,m=calculate(ij,np.array([3.,2.,2.]),d,0,6,2,'topL')
        self.assertEqual(r['n_selected'],3)
        self.assertAlmostEqual(r['precision'],2/3)
        self.assertAlmostEqual(r['recall'],1)
        self.assertAlmostEqual(r['jaccard'],2/3)

    def test_missing_atoms_and_undefined_L(self):
        d=np.full((8,8),np.nan)
        r,_=calculate(np.array([[0,7]]),np.array([4.]),d,0,6,1,'topL')
        self.assertEqual(r['status'],'insufficient_eligible_for_L')
        self.assertTrue(np.isnan(r['precision']))

    def test_offset_boundary_and_separation(self):
        d=np.full((9,9),7.)
        r,_=calculate(np.array([[0,6],[1,7],[2,8],[0,5]]),np.array([1.,2.,4.,9.]),d,1,6,2,'topL')
        self.assertEqual(r['n_eligible'],2)
        self.assertEqual(r['precision'],1)

    def test_infer_crop_without_contacts(self):
        seq='ACDEFGHIKLMNPQRSTVWYACDEFGHIKL'
        native=[('query',seq)]+[(str(i),seq) for i in range(3)]
        swap=[('query',seq[2:-2])]+[(str(i),seq[2:-2]) for i in range(3)]
        r=infer_offset(swap,native)
        self.assertEqual(r['status'],'ok'); self.assertEqual(r['offset'],2)

    def test_ambiguous_registration_excluded(self):
        rows=[('query','A'*30)]+[(str(i),'A'*30) for i in range(3)]
        self.assertNotEqual(infer_offset(rows,rows)['status'],'ok')


if __name__=='__main__':
    unittest.main()
