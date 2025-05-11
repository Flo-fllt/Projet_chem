import unittest

from Package_functions.Model_training_functions import smiles_to_fingerprints

class TestSmilesToFingerprints(unittest.TestCase):
    def test_valid_reaction(self):
        rxn = "CCO.CN>>CC=O.CO"
        r_fps, p_fps = smiles_to_fingerprints(rxn)
        self.assertEqual(len(r_fps), 2)
        self.assertEqual(len(p_fps), 2)
        for fp in r_fps + p_fps:
            self.assertEqual(fp.shape[0], 2048)
            self.assertTrue((fp == 0).sum() < 2048)

    def test_invalid_molecule_in_reaction(self):
        rxn = "CCO.invalid>>CC=O"
        r_fps, p_fps = smiles_to_fingerprints(rxn)
        self.assertEqual(len(r_fps), 1)  # CCO is valid, "invalid" is skipped
        self.assertEqual(len(p_fps), 1)

    def test_invalid_reaction_format(self):
        rxn = "CCO.CN>CC=O"  # missing one >
        r_fps, p_fps = smiles_to_fingerprints(rxn)
        self.assertEqual(r_fps, [])
        self.assertEqual(p_fps, [])

    def test_empty_input(self):
        r_fps, p_fps = smiles_to_fingerprints("")
        self.assertEqual((r_fps, p_fps), ([], []))

if __name__ == "__main__":
    unittest.main()
