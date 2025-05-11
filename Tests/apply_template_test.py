import unittest

# Import function under test
from Package_functions.Interface_functions import apply_template

class TestApplyTemplate(unittest.TestCase):
    def test_valid_reaction(self):
        # Simple ester hydrolysis template
        template = "[C:1](=O)[O:2][C:3]>>[C:1](=O)[O:2].[C:3][OH]"
        input_smiles = "CC(=O)OC"
        result = apply_template(template, input_smiles)
        self.assertIsInstance(result, list)
        self.assertGreater(len(result), 0)
        self.assertTrue(all(isinstance(prod, list) for prod in result))

    def test_invalid_smiles(self):
        result = apply_template("[C:1]>>[C:1]", "not_a_smiles")
        self.assertEqual(result, [])

    def test_invalid_template(self):
        result = apply_template("invalid_template", "CCO")
        self.assertEqual(result, [])

if __name__ == "__main__":
    unittest.main()
