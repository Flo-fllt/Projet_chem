import unittest
import pandas as pd

#Import function under test 
from RetroChem.Package_functions.Model_training_functions import prepare_fingerprints_for_training

class TestPrepareFingerprints(unittest.TestCase):

    def test_valid_dataframe(self):
        df = pd.DataFrame({
            "RxnSmilesClean": ["CCO.CN>>CC=O", "CCC>>CC=C"],
            "TemplateHash": ["template_1", "template_2"]
        })
        X, y = prepare_fingerprints_for_training(df)

        expected_total = 3 + 2

        self.assertEqual(X.shape[0], expected_total)
        self.assertEqual(y.shape[0], expected_total)
        self.assertTrue(all(len(fp) == 2048 for fp in X), "All fingerprints must be of length 2048.")

    def test_with_invalid_smiles(self):
        df = pd.DataFrame({
            "RxnSmilesClean": ["INVALID>>CC=O", "CCC>>CC=C"],
            "TemplateHash": ["bad", "template"]
        })
        X, y = prepare_fingerprints_for_training(df)

        expected_total = 2

        self.assertEqual(X.shape[0], expected_total)
        self.assertEqual(y.shape[0], expected_total)

    def test_empty_dataframe(self):
        df = pd.DataFrame(columns=["RxnSmilesClean", "TemplateHash"])
        X, y = prepare_fingerprints_for_training(df)

        self.assertEqual(X.shape[0], 0)
        self.assertEqual(y.shape[0], 0)

    def test_partial_invalid_reaction(self):
        df = pd.DataFrame({
            "RxnSmilesClean": ["CCO.INVALID>>CC=O", "CN>>INVALID"],
            "TemplateHash": ["t1", "t2"]
        })
        X, y = prepare_fingerprints_for_training(df)

        expected_total = 2

        self.assertEqual(X.shape[0], expected_total)
        self.assertEqual(y.shape[0], expected_total)

if __name__ == "__main__":
    unittest.main()

