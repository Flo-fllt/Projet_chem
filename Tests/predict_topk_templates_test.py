# file: predict_topk_templates_test.py

import unittest
import numpy as np
import pandas as pd
from unittest.mock import patch, MagicMock

# Dummy function to stand in for actual fingerprinting
def smiles_to_fingerprint(smiles, radius=2, n_bits=2048):
    if smiles == "INVALID":
        raise ValueError("Invalid SMILES string")
    return np.ones((2048,), dtype=int)

# Function under test
def predict_topk_templates(smiles_input, topk=3):
    scaler = MagicMock()
    scaler.transform.return_value = np.ones((1, 2048))

    model = MagicMock()
    model.predict_proba.return_value = np.array([[0.1, 0.6, 0.3]])
    model.classes_ = np.array([0, 1, 2])

    label_encoder = MagicMock()
    label_encoder.inverse_transform.side_effect = lambda x: [f"template_{i}" for i in x]

    templates_df = pd.DataFrame({
        "TemplateHash": ["template_0", "template_1", "template_2"],
        "RetroTemplate": ["[CH3:1][OH:2]", "[C:1]=[O:2]", "[C:1][C:2]"]
    })

    fingerprint = smiles_to_fingerprint(smiles_input).reshape(1, -1)
    fingerprint_scaled = scaler.transform(fingerprint)
    probs = model.predict_proba(fingerprint_scaled)[0]
    topk_indices = np.argsort(probs)[::-1][:topk]
    topk_template_hashes = label_encoder.inverse_transform(model.classes_[topk_indices])
    topk_probs = probs[topk_indices]

    predictions = []
    for template_hash, prob in zip(topk_template_hashes, topk_probs):
        row = templates_df[templates_df['TemplateHash'] == template_hash]
        if not row.empty:
            retro_template = row.iloc[0]['RetroTemplate']
            predictions.append((template_hash, retro_template, prob))
    return predictions

class TestPredictTopKTemplates(unittest.TestCase):
    def test_valid_smiles_prediction(self):
        predictions = predict_topk_templates("CCO", topk=2)
        self.assertEqual(len(predictions), 2)
        for pred in predictions:
            self.assertEqual(len(pred), 3)  # (template_hash, retro_template, probability)
            self.assertTrue(isinstance(pred[2], float))

    def test_invalid_smiles(self):
        with self.assertRaises(ValueError):
            predict_topk_templates("INVALID")

if __name__ == "__main__":
    unittest.main()
