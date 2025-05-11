# file: smiles_to_fingerprints_test.py

import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdchem, DataStructs

# Dependencies 
def split_rxn_smiles(rxn_smiles):
    parts = rxn_smiles.split(">>")
    if len(parts) != 2:
        return [], []
    return parts[0].split("."), parts[1].split(".")

def remove_atom_mapping(smiles):
    import re
    return re.sub(r":\d+", "", smiles)

# Function under test
def smiles_to_fingerprints(smiles):
    try:
        reactants_smiles, products_smiles = split_rxn_smiles(smiles)

        def mols_to_fps(smiles_list):
            fps = []
            n_bits = 2048
            for s in smiles_list:
                cleaned_s = remove_atom_mapping(s)
                mol = Chem.MolFromSmiles(cleaned_s)
                if mol is None:
                    print(f"Molécule invalide : {s}")
                    continue
                fp = AllChem.GetMorganFingerprint(mol, radius=3)
                arr = np.zeros((n_bits,), dtype=int)
                if isinstance(fp, DataStructs.UIntSparseIntVect):
                    on_bits = list(fp.GetNonzeroElements().keys())
                    for bit in on_bits:
                        arr[bit % n_bits] = 1
                fps.append(arr)
            return fps

        reactants_fps = mols_to_fps(reactants_smiles)
        products_fps = mols_to_fps(products_smiles)

        return reactants_fps, products_fps

    except Exception as e:
        print(f"Erreur parsing SMILES: {smiles} -> {e}")
        return [], []

class TestSmilesToFingerprints(unittest.TestCase):
    def test_valid_reaction(self):
        rxn = "CCO.CN>>CC=O.CO"
        r_fps, p_fps = smiles_to_fingerprints(rxn)
        self.assertEqual(len(r_fps), 2)
        self.assertEqual(len(p_fps), 2)
        for fp in r_fps + p_fps:
            self.assertEqual(fp.shape[0], 2048)
            self.assertTrue((fp == 0).sum() < 2048)

    def test_invalid_molecule_in_reaction(self):
        rxn = "CCO.invalid>>CC=O"
        r_fps, p_fps = smiles_to_fingerprints(rxn)
        self.assertEqual(len(r_fps), 1)  # CCO is valid, "invalid" is skipped
        self.assertEqual(len(p_fps), 1)

    def test_invalid_reaction_format(self):
        rxn = "CCO.CN>CC=O"  # missing one >
        r_fps, p_fps = smiles_to_fingerprints(rxn)
        self.assertEqual(r_fps, [])
        self.assertEqual(p_fps, [])

    def test_empty_input(self):
        r_fps, p_fps = smiles_to_fingerprints("")
        self.assertEqual((r_fps, p_fps), ([], []))

if __name__ == "__main__":
    unittest.main()
