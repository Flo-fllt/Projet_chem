import unittest
from rdkit import Chem
from PIL import Image

# Import function under test
from Package_functions.Interface_functions import mol_to_high_quality_image

class TestMolToHighQualityImage(unittest.TestCase):
    def test_image_output(self):
        mol = Chem.MolFromSmiles("CCO")  # ethanol
        img = mol_to_high_quality_image(mol)
        self.assertIsInstance(img, Image.Image)
        self.assertEqual(img.size, (800, 800))

    def test_invalid_molecule(self):
        # If mol is None, there should be an error
        with self.assertRaises(Exception):
            mol_to_high_quality_image(None)

if __name__ == '__main__':
    unittest.main()
