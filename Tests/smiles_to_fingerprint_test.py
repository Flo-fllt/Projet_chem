# fichier: smiles_to_fingerprint_test.py

import unittest
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# Function to test
def smiles_to_fingerprint(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((1,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

class TestSmilesToFingerprint(unittest.TestCase):
    def test_valid_smiles(self):
        arr = smiles_to_fingerprint("CCO")  # ethanol
        self.assertEqual(arr.shape[0], 2048)
        self.assertTrue((arr == 0).sum() < 2048)  # at least one bit activated

    def test_invalid_smiles(self):
        with self.assertRaises(ValueError):
            smiles_to_fingerprint("not_a_smiles")

if __name__ == '__main__':
    unittest.main()
