import unittest
from PIL import Image

#import function under test 
from Package_functions.Interface_functions import st_scaled_image

class TestStScaledImage(unittest.TestCase):
    def test_valid_image_display(self):
        img = Image.new("RGB", (100, 100), color="red")
        try:
            st_scaled_image(img, width_display_px=150)
        except Exception as e:
            self.fail(f"st_scaled_image raised an exception with a valid image: {e}")

    def test_invalid_image_input(self):
        with self.assertRaises(AttributeError):  # e.g., 'str' has no 'save' method
            st_scaled_image("not_an_image")

if __name__ == "__main__":
    unittest.main()
