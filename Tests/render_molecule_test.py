# file: render_molecule_test.py

import unittest
from rdkit import Chem
from PIL import Image

# Mocked image-rendering function
def mol_to_high_quality_image(mol, size=(800, 800)):
    return Image.new("RGB", size, color="white")

# Function under test
def render_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return mol_to_high_quality_image(mol)
    return None

class TestRenderMolecule(unittest.TestCase):
    def test_valid_smiles(self):
        result = render_molecule("CCO")
        self.assertIsInstance(result, Image.Image)
        self.assertEqual(result.size, (800, 800))

    def test_invalid_smiles(self):
        result = render_molecule("not_a_smiles")
        self.assertIsNone(result)

if __name__ == "__main__":
    unittest.main()
