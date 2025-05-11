# file: render_reaction_scheme_test.py

import unittest
from PIL import Image

# Dummy helpers (replace with real imports in your codebase)
def render_molecule(smi):
    return Image.new("RGB", (100, 100), color="white")

def image_to_base64(img, width=250):
    return f'<img src="data:image/png;base64,DUMMY" width="{width}"/>'

# Function to test
def render_reaction_scheme(smiles_chain):
    html_parts = []

    for i, smi_group in enumerate(smiles_chain):
        mol_imgs = [render_molecule(smi) for smi in smi_group]
        img_htmls = [image_to_base64(img, width=250) for img in mol_imgs if img]
        html_parts.append(" + ".join(img_htmls))
        if i < len(smiles_chain) - 1:
            html_parts.append('<span style="font-size: 28px; font-weight: bold; margin: 0 12px;">→</span>')

    full_html = "".join(html_parts)

    return f'''
    <div style="display: flex; flex-wrap: nowrap; align-items: center; overflow-x: auto;">
        {full_html}
    </div>
    '''

class TestRenderReactionScheme(unittest.TestCase):
    def test_output_contains_html(self):
        smiles_chain = [["CCO"], ["CC=O"]]
        html = render_reaction_scheme(smiles_chain)
        self.assertIsInstance(html, str)
        self.assertIn("→", html)
        self.assertIn("<img", html)
        self.assertTrue(len(html.strip()) > 0)

    def test_single_step(self):
        smiles_chain = [["CCO"]]
        html = render_reaction_scheme(smiles_chain)
        self.assertNotIn("→", html)

    def test_empty_input(self):
        html = render_reaction_scheme([])
        self.assertIn("flex", html)  # basic container still rendered

if __name__ == "__main__":
    unittest.main()
