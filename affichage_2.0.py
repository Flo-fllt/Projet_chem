import streamlit as st
from streamlit_ketcher import st_ketcher
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem, rdChemReactions
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
from io import BytesIO

# --- High quality molecule drawing ---
def mol_to_high_quality_image(mol, size=(700, 700)):
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.0
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    return Image.open(BytesIO(png))

# --- Load model and data ---
scaler = joblib.load("scaler.pkl")
model = joblib.load("mlp_classifier_model.pkl")
label_encoder = joblib.load("label_encoder.pkl")
templates_df = pd.read_csv("/Users/giuliogarotti/Documents/GitHub/Projet_chem/uspto50/uspto50/combined_data.csv", sep="\t")

# --- Functions ---
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
    fingerprint = smiles_to_fingerprint(smiles_input)
    fingerprint = fingerprint.reshape(1, -1)
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

# --- Streamlit App ---
st.title("🧪 RetroSynthesis Prediction Tool")
with st.expander("1. Draw Molecule"):
    molecule = st.text_input("**Paste SMILES (optional)**")
    ketcher_smiles = st_ketcher(molecule, height=600)

final_smiles = ketcher_smiles or molecule
if final_smiles:
    st.success(f"✅ SMILES: `{final_smiles}`")

if st.button("Run Retrosynthesis") and final_smiles:
    try:
        st.info("🔍 Predicting templates and generating precursors...")
        topk_predictions = predict_topk_templates(final_smiles, topk=50)

        seen_reactants = set()
        successful_predictions = []

        for rank, (template_hash, retro_template, prob) in enumerate(topk_predictions, start=1):
            predicted_reactants = apply_template(retro_template, final_smiles)
            for prod_set in predicted_reactants:
                canon_prod_set = tuple(sorted(prod_set))
                if canon_prod_set not in seen_reactants:
                    seen_reactants.add(canon_prod_set)
                    successful_predictions.append((template_hash, retro_template, prob, prod_set))

        if successful_predictions:
            st.success("🎯 Successful Predictions:")
            for idx, (template_hash, smarts, prob, reactants) in enumerate(successful_predictions, 1):
                with st.expander(f"🔹 Prediction {idx} - {prob*100:.2f}% confidence"):
                    st.markdown(f"**Template Hash:** `{template_hash}`")
                    st.markdown(f"**Template SMARTS:** `{smarts}`")
                    st.markdown("**Predicted Reactants:**")
                    cols = st.columns(len(reactants))
                    for i, smi in enumerate(reactants):
                        mol = Chem.MolFromSmiles(smi)
                        if mol:
                            img = mol_to_high_quality_image(mol)
                            cols[i].image(img, caption=smi, use_container_width=True)
        else:
            st.error("❌ No valid templates produced any reactants.")

    except Exception as e:
        st.error(f"⚠️ Error: {e}")
