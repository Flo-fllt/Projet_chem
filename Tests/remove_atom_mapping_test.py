# file: remove_atom_mapping_test.py

import unittest
import re

# Function under test
def remove_atom_mapping(smiles):
    """Remove :number mapping from SMILES."""
    return re.sub(r":\d+", "", smiles)

class TestRemoveAtomMapping(unittest.TestCase):
    def test_removes_atom_mappings(self):
        self.assertEqual(remove_atom_mapping("[CH3:1][CH2:2][OH:3]"), "[CH3][CH2][OH]")
        self.assertEqual(remove_atom_mapping("[C:12](=[O:4])[O-:7]"), "[C](=[O])[O-]")

    def test_smiles_without_mappings(self):
        self.assertEqual(remove_atom_mapping("CC(=O)O"), "CC(=O)O")

    def test_empty_input(self):
        self.assertEqual(remove_atom_mapping(""), "")

    def test_nested_or_complex_smiles(self):
        self.assertEqual(
            remove_atom_mapping("[C@@H:1]([O:2])[N:3](C)[C:4]"),
            "[C@@H]([O])[N](C)[C]"
        )

if __name__ == "__main__":
    unittest.main()
