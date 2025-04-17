import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score
import re
from rdkit import Chem
from rdkit.Chem import AllChem

# 1. Load the saved model and scaler
model = joblib.load('mlp_classifier_model.pkl')
scaler = joblib.load('scaler.pkl')


# 2. Load your test data (or reuse X, y from before if you don't have separate test set)
file_path = "/Users/giuliogarotti/Documents/GitHub/Projet_chem/uspto50/uspto50/reaction_templates_50k_test.csv"
df = pd.read_csv(file_path, sep="\t")

# Prepare the data again (same as training)
def fix_common_errors(smiles):
    smiles = re.sub(r'@', 'O', smiles)
    return smiles

def split_rxn_smiles(rxn_smiles):
    try:
        rxn_smiles = fix_common_errors(rxn_smiles) 
        parts = rxn_smiles.split(">>")
        if len(parts) != 2:
            return []
        reactants_raw, products_raw = parts
        reactants = reactants_raw.split(".")
        products = products_raw.split(".")
        return reactants, products
    except:
        return []

def smiles_to_fingerprints(smiles):
    try:
        reactants_smiles, products_smiles = split_rxn_smiles(smiles)
        def mols_to_fps(smiles_list):
            fps = []
            n_bits = 2048
            for s in smiles_list:
                mol = Chem.MolFromSmiles(s)
                if mol is None:
                    fps.append(np.zeros(n_bits))
                else:
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
    except:
        return []

def prepare_fingerprints_for_testing(df):
    X_test = []
    y_test = []
    for smiles, target in zip(df['RxnSmilesClean'], df['TemplateHash']):
        fps = smiles_to_fingerprints(smiles)
        if fps:
            X_test.extend(fps)
            y_test.extend([target] * len(fps))
    return np.array(X_test), np.array(y_test)

X_test, y_test = prepare_fingerprints_for_testing(df)

# 3. Scale test data using the SAME scaler
X_test_scaled = scaler.transform(X_test)

# 4. Predict with the loaded model
y_pred = model.predict(X_test_scaled)

# 5. Evaluate
acc = accuracy_score(y_test, y_pred)
print(f"Test Accuracy: {acc*100:.2f}%")
