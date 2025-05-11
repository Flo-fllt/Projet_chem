import unittest

# Import Function under test
from Package_functions.Model_training_functions import split_rxn_smiles

class TestSplitRxnSmiles(unittest.TestCase):
    def test_valid_reaction(self):
        rxn = "CCO.CN>>CC=O.CO"
        reactants, products = split_rxn_smiles(rxn)
        self.assertEqual(reactants, ["CCO", "CN"])
        self.assertEqual(products, ["CC=O", "CO"])

    def test_single_reactant_product(self):
        rxn = "CCO>>CC=O"
        reactants, products = split_rxn_smiles(rxn)
        self.assertEqual(reactants, ["CCO"])
        self.assertEqual(products, ["CC=O"])

    def test_invalid_format_missing_arrow(self):
        rxn = "CCO.CN>CC=O.CO"
        reactants, products = split_rxn_smiles(rxn)
        self.assertEqual(reactants, [])
        self.assertEqual(products, [])

    def test_empty_string(self):
        reactants, products = split_rxn_smiles("")
        self.assertEqual((reactants, products), ([], []))

    def test_extra_dots(self):
        rxn = "CCO.>>.CC=O"
        reactants, products = split_rxn_smiles(rxn)
        self.assertEqual(reactants, ["CCO", ""])
        self.assertEqual(products, ["", "CC=O"])

if __name__ == "__main__":
    unittest.main()
