# file: apply_template_test.py

import unittest
from rdkit import Chem
from rdkit.Chem import rdChemReactions

# Function under test
def apply_template(template_smarts, smiles_input):
    mol = Chem.MolFromSmiles(smiles_input)
    if mol is None:
        return []
    try:
        rxn = rdChemReactions.ReactionFromSmarts(template_smarts)
        products = rxn.RunReactants((mol,))
        product_smiles = []
        for prod_set in products:
            prod_list = [Chem.MolToSmiles(p) for p in prod_set if p is not None]
            if prod_list:
                product_smiles.append(prod_list)
        return product_smiles
    except:
        return []

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
