import joblib
import numpy as np
import pandas as pd
import re
from sklearn.metrics import accuracy_score
from rdkit import Chem
from rdkit.Chem import AllChem

# --- 1. Load the saved model and scaler ---
model = joblib.load('mlp_classifier_model.pkl')
scaler = joblib.load('scaler.pkl')

# --- 2. Load your test data ---
file_path = r"C:\Users\flopi\projet_reaction_chimie\uspto50\uspto50\reaction_templates_50k_test.csv"
df = pd.read_csv(file_path, sep="\t")

# --- 3. Prepare the data (fixed version) ---

def remove_atom_mapping(smiles):
    """Remove :number atom mappings from SMILES."""
    return re.sub(r":\d+", "", smiles)

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
        print(f"Erreur split_rxn_smiles: {rxn_smiles} -> {e}")
        return [], []

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
                    continue  # skip invalid molecules
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

        return reactants_fps + products_fps

    except Exception as e:
        print(f"Erreur parsing SMILES: {smiles} -> {e}")
        return []

def prepare_fingerprints_for_testing(df):
    X_test = []
    y_test = []
    for idx, (smiles, target) in enumerate(zip(df['RxnSmilesClean'], df['TemplateHash'])):
        fps = smiles_to_fingerprints(smiles)
        if fps:
            X_test.extend(fps)
            y_test.extend([target] * len(fps))
        else:
            print(f"Skipping invalid SMILES at index {idx}: {smiles}")
    return np.array(X_test), np.array(y_test)

# --- 4. Prepare the clean X_test and y_test ---
X_test, y_test = prepare_fingerprints_for_testing(df)

print(f"X_test shape: {X_test.shape}, y_test shape: {y_test.shape}")

print(len(np.unique(y_test)))

# --- 5. Scale test data ---
X_test_scaled = scaler.transform(X_test)

# --- 6. Predict with the loaded model ---
y_pred = model.predict(X_test_scaled)

# --- 7. Evaluate ---
acc = accuracy_score(y_test, y_pred)
print(f"✅ Test Accuracy: {acc*100:.2f}%")
