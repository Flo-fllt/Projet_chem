import unittest

# Import function to test
from Package_retrosynth.Package_functions.Interface_functions import render_reaction_scheme

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
