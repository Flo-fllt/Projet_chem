import unittest
from PIL import Image

# Import Function under test
from Package_retrosynth.Package_functions.Interface_functions import image_to_base64

class TestImageToBase64(unittest.TestCase):
    def test_valid_image(self):
        img = Image.new("RGB", (100, 100), color="blue")
        html = image_to_base64(img, width=300)
        self.assertIsInstance(html, str)
        self.assertIn('<img src="data:image/png;base64,', html)
        self.assertIn('width="300px"', html)
        self.assertTrue(html.strip().endswith("/>") or html.strip().endswith(">"))

    def test_invalid_input(self):
        with self.assertRaises(AttributeError):
            image_to_base64("not_an_image")

if __name__ == "__main__":
    unittest.main()
