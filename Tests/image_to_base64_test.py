# file: test_image_to_base64.py

import unittest
from PIL import Image
from io import BytesIO
import base64

# Function under test
def image_to_base64(img, width=350):
    buffer = BytesIO()
    img.save(buffer, format="PNG")
    encoded = base64.b64encode(buffer.getvalue()).decode("utf-8")
    return f'''
    <div style="display:inline-block; margin: 0 10px;">
        <img src="data:image/png;base64,{encoded}" width="{width}px" />
    </div>
    '''

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
