# file: st_scaled_image_test.py

import unittest
from PIL import Image
from io import BytesIO
import base64
import streamlit as st

# Function under test
def st_scaled_image(image, width_display_px=200):
    buffered = BytesIO()
    image.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode()
    html = f"""
    <div style="display:inline-block;">
        <img src="data:image/png;base64,{img_str}" style="width:{width_display_px}px; height:auto;" />
    </div>
    """
    st.markdown(html, unsafe_allow_html=True)

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
