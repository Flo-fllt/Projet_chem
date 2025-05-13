import unittest

# Import Function to test
from Package_retrosynth.Package_functions.Interface_functions import smiles_to_fingerprint

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
