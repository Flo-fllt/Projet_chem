# file: mol_to_high_quality_image_test.py

import unittest
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
from io import BytesIO

# Function under test
def mol_to_high_quality_image(mol, size=(800, 800)):
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.0
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    return Image.open(BytesIO(png))

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
