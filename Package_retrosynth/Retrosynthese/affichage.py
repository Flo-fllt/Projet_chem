import streamlit as st
from streamlit_ketcher import st_ketcher
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, rdChemReactions
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
from io import BytesIO
import base64
import os
import pandas as pd

# --- High quality molecule drawing ---
def mol_to_high_quality_image(mol, size=(800, 800)):
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.0
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    return Image.open(BytesIO(png))

# --- Display scaled image (visually smaller but full resolution) ---
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


# --- Helper functions ---
def smiles_to_fingerprint(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((1,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def apply_template(template_smarts, smiles_input):
    mol = Chem.MolFromSmiles(smiles_input)
    if mol is None:
        return []
    try:
        rxn = rdChemReactions.ReactionFromSmarts(template_smarts)
        products = rxn.RunReactants((mol,))
        product_smiles = []
        for prod_set in products:
            prod_list = [Chem.MolToSmiles(p) for p in prod_set if p is not None]
            if prod_list:
                product_smiles.append(prod_list)
        return product_smiles
    except:
        return []

def predict_topk_templates(smiles_input, topk=50):
    scaler = joblib.load("scaler.pkl")
    model = joblib.load("mlp_classifier_model.pkl")
    label_encoder = joblib.load("label_encoder.pkl")
    # Get the directory of the current file (affichage.py)
    base_dir = os.path.dirname(__file__)

# Build path to Data/combined_data.csv
    csv_path = os.path.join(base_dir, "Data", "combined_data.csv")

# Load the CSV
    templates_df = pd.read_csv(csv_path, sep="\t")
    
    fingerprint = smiles_to_fingerprint(smiles_input).reshape(1, -1)
    fingerprint_scaled = scaler.transform(fingerprint)
    probs = model.predict_proba(fingerprint_scaled)[0]
    topk_indices = np.argsort(probs)[::-1][:topk]
    topk_template_hashes = label_encoder.inverse_transform(model.classes_[topk_indices])
    topk_probs = probs[topk_indices]

    predictions = []
    for template_hash, prob in zip(topk_template_hashes, topk_probs):
        row = templates_df[templates_df['TemplateHash'] == template_hash]
        if not row.empty:
            retro_template = row.iloc[0]['RetroTemplate']
            predictions.append((template_hash, retro_template, prob))
    return predictions