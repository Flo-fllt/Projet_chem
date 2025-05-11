# file: prepare_fingerprints_for_traning_test.py

import unittest
import pandas as pd
import numpy as np

# Mocked dependency
def smiles_to_fingerprints(smiles):
    if "INVALID" in smiles:
        return [], []
    # return 2 fake reactants and 1 fake product fingerprint
    fake_fp = np.zeros((2048,), dtype=int)
    fake_fp[0] = 1
    return [fake_fp, fake_fp], [fake_fp]

# Function under test
def prepare_fingerprints_for_training(df):
    X = []
    y = []

    print("Début du traitement des données.")
    
    for idx, (smiles, target) in enumerate(zip(df['RxnSmilesClean'], df['TemplateHash'])):
        if idx < 5:
            print(f"Index {idx} - SMILES: {smiles} | Target: {target}")
        
        reactants_fps, products_fps = smiles_to_fingerprints(smiles)
        
        if reactants_fps and products_fps:
            X.extend(reactants_fps)
            y.extend([target] * len(reactants_fps))
            X.extend(products_fps)
            y.extend([target] * len(products_fps))
        else:
            print(f"Skipping reaction at index {idx}: {smiles}")

    X = np.array(X)
    y = np.array(y)

    print(f"Fingerprint preparation finished. Total examples: {X.shape[0]}")
    if X.shape[0] != y.shape[0]:
        raise ValueError(f"Mismatch: X {X.shape[0]} vs y {y.shape[0]}")

    return X, y

class TestPrepareFingerprints(unittest.TestCase):
    def test_valid_dataframe(self):
        df = pd.DataFrame({
            "RxnSmilesClean": ["CCO.CN>>CC=O", "CCC>>CC=C"],
            "TemplateHash": ["template_1", "template_2"]
        })
        X, y = prepare_fingerprints_for_training(df)
        # Each entry returns 2 reactant + 1 product = 3 fingerprints
        self.assertEqual(X.shape[0], 6)
        self.assertEqual(y.shape[0], 6)
        self.assertTrue(all(len(fp) == 2048 for fp in X))

    def test_with_invalid_smiles(self):
        df = pd.DataFrame({
            "RxnSmilesClean": ["INVALID>>CC=O", "CCC>>CC=C"],
            "TemplateHash": ["bad", "template"]
        })
        X, y = prepare_fingerprints_for_training(df)
        # Only one valid entry (CCC>>CC=C): 2 reactant + 1 product = 3
        self.assertEqual(X.shape[0], 3)
        self.assertEqual(y.shape[0], 3)

    def test_empty_dataframe(self):
        df = pd.DataFrame(columns=["RxnSmilesClean", "TemplateHash"])
        X, y = prepare_fingerprints_for_training(df)
        self.assertEqual(X.shape[0], 0)
        self.assertEqual(y.shape[0], 0)

if __name__ == "__main__":
    unittest.main()
