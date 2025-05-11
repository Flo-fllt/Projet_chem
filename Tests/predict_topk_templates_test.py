import unittest

# import function under test
from Package_functions.Interface_functions import predict_topk_templates

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
