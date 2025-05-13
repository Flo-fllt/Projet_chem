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
import importlib.resources as pkg_resources
from pathlib import Path
import io
import RetroChem.Model as model_pkg
import RetroChem.Data as data_pkg


# --- Helper functions ---
def smiles_to_fingerprint(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((1,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def render_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return mol_to_high_quality_image(mol)
    return None

def image_to_base64(img, width=350):
    buffer = BytesIO()
    img.save(buffer, format="PNG")
    encoded = base64.b64encode(buffer.getvalue()).decode("utf-8")
    return f'''
    <div style="display:inline-block; margin: 0 10px;">
        <img src="data:image/png;base64,{encoded}" width="{width}px" />
    </div>
    '''

# --- High quality molecule drawing ---
def mol_to_high_quality_image(mol, size=(800, 800)):
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.0
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    return Image.open(BytesIO(png))

# --- Display scaled image ---
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

# --- Reaction scheme rendering ---
def render_reaction_scheme(smiles_chain):
    html_parts = []

    for i, smi_group in enumerate(smiles_chain):
        mol_imgs = [render_molecule(smi) for smi in smi_group]
        img_htmls = [image_to_base64(img, width=250) for img in mol_imgs if img]
        html_parts.append(" + ".join(img_htmls))
        if i < len(smiles_chain) - 1:
            html_parts.append('<span style="font-size: 28px; font-weight: bold; margin: 0 12px;">→</span>')

    full_html = "".join(html_parts)

    # Wrap the full scheme in a flex container to keep it on one line
    return f'''
    <div style="display: flex; flex-wrap: nowrap; align-items: center; overflow-x: auto;">
        {full_html}
    </div>
    '''


# --- Prediction functions ---
def predict_topk_templates(smiles_input, topk=50):

    # --- Load model components from Package_retrosynth/Model/ ---
    with pkg_resources.files(model_pkg).joinpath("scaler.pkl").open("rb") as f:
        scaler = joblib.load(f)
    with pkg_resources.files(model_pkg).joinpath("mlp_classifier_model.pkl").open("rb") as f:
        model = joblib.load(f)
    with pkg_resources.files(model_pkg).joinpath("label_encoder.pkl").open("rb") as f:
        label_encoder = joblib.load(f)


    # --- Load CSV from Package_retrosynth/Data/ ---
    with pkg_resources.files(data_pkg).joinpath("combined_data.csv").open("r", encoding="utf-8") as f:
        templates_df = pd.read_csv(f, sep="\t")
    
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
    