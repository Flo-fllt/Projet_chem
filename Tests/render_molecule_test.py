import unittest
from PIL import Image

from Package_functions.Interface_functions import render_molecule

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
