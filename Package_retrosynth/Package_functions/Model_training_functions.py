import re
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
import joblib
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedShuffleSplit

# --- 1. Remove atom mapping ---
def remove_atom_mapping(smiles):
    """Remove :number mapping from SMILES."""
    return re.sub(r":\d+", "", smiles)

# --- 2. Split SMILES into reactants and products ---
def split_rxn_smiles(rxn_smiles):
    try:
        parts = rxn_smiles.split(">>")
        if len(parts) != 2:
            return [], []
        reactants_raw, products_raw = parts
        reactants = reactants_raw.split(".")
        products = products_raw.split(".")
        return reactants, products
    except Exception as e:
        print(f"Erreur dans split_rxn_smiles: {rxn_smiles} -> {e}")
        return [], []
    
# --- 3. Convert SMILES to fingerprints (after cleaning) ---
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
                    continue  # <-- SKIP invalid molecules
                fp = AllChem.GetMorganFingerprint(mol, radius=3)
                arr = np.zeros((n_bits,), dtype=int)
                if isinstance(fp, Chem.DataStructs.UIntSparseIntVect):
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

# --- 4. Prepare training data X, y ---
def prepare_fingerprints_for_training(df):
    X = []
    y = []
    
    print("Start of data processing.")
    
    for idx, (smiles, target) in enumerate(zip(df['RxnSmilesClean'], df['TemplateHash'])):
        if idx < 5:  
            print(f"Index {idx} - SMILES: {smiles} | Target: {target}")
        
        reactants_fps, products_fps = smiles_to_fingerprints(smiles)
        
        if reactants_fps and products_fps:
            X.extend(reactants_fps)
            y.extend([target] * len(reactants_fps))
            X.extend(products_fps)
            y.extend([target] * len(products_fps))
        else:
            print(f"Skipping reaction at index {idx}: {smiles}")

    X = np.array(X)
    y = np.array(y)

    print(f"Fingerprint preparation finished. Total examples: {X.shape[0]}")
    if X.shape[0] != y.shape[0]:
        raise ValueError(f"Mismatch: X {X.shape[0]} vs y {y.shape[0]}")

    return X, y