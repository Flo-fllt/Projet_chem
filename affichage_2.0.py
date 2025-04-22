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

# --- High quality molecule drawing ---
def mol_to_high_quality_image(mol, size=(300, 300)):
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.0
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    return Image.open(BytesIO(png))

# --- Load model and data ---
scaler = joblib.load('scaler.pkl')
model = joblib.load('mlp_classifier_model.pkl')
templates_df = pd.read_csv("/Users/giuliogarotti/Documents/GitHub/Projet_chem/uspto50/uspto50/reaction_templates_50k_train.csv", sep='\t')

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
    topk_template_hashes = model.classes_[topk_indices]
    topk_probs = probs[topk_indices]

    predictions = []
    for template_hash, prob in zip(topk_template_hashes, topk_probs):
        row = templates_df[templates_df['TemplateHash'] == template_hash]
        if not row.empty:
            retro_template = row.iloc[0]['RetroTemplate']
            predictions.append((template_hash, retro_template, prob))
        else:
            predictions.append((template_hash, None, prob))
    return predictions

# --- Streamlit App ---
st.title("🧪 RetroSynthesis Prediction Tool")
with st.expander("1. Draw Molecule"):
    molecule = st.text_input("**Paste SMILES (optional)**")
    ketcher_smiles = st_ketcher(molecule, height=600)

# Show the generated SMILES
final_smiles = ketcher_smiles or molecule
if final_smiles:
    st.success(f"✅ SMILES: `{final_smiles}`")
else:
    st.warning("Draw or enter a SMILES string to begin.")

# Run retrosynthesis when button is clicked
if st.button("Run Retrosynthesis") and final_smiles:
    try:
        st.info("🔍 Predicting templates and generating precursors...")
        topk_predictions = predict_topk_templates(final_smiles, topk=50)
        
        successful_predictions = []
        seen_reactants = set()  # Track unique reactant sets
        
        for rank, (template_hash, retro_template, prob) in enumerate(topk_predictions, start=1):
            if retro_template is None:
                continue

            predicted_products = apply_template(retro_template, final_smiles)

            if predicted_products:
                for prod_set in predicted_products:
                    canon_prod_set = tuple(sorted(prod_set))  # Canonical form
                    if canon_prod_set not in seen_reactants:
                        seen_reactants.add(canon_prod_set)
                        successful_predictions.append((template_hash, retro_template, prob, [prod_set]))

        if successful_predictions:
            st.success("🎯 Successful Retrosynthesis Predictions:")
            for idx, (template_hash, retro_template, prob, products) in enumerate(successful_predictions, start=1):
                with st.expander(f"🔹 Prediction {idx} - {prob*100:.2f}% confidence"):
                    st.markdown(f"**Template Hash:** `{template_hash}`")
                    st.markdown(f"**Template SMARTS:** `{retro_template}`")
                    st.markdown("**Predicted Reactants:**")
                    for prod_set in products:
                        st.write("**Reactant Set:**")
                        cols = st.columns(len(prod_set))
                        for i, smiles in enumerate(prod_set):
                            mol = Chem.MolFromSmiles(smiles)
                            if mol:
                                img = mol_to_high_quality_image(mol, size=(700, 700))
                                cols[i].image(img, caption=smiles, use_container_width=True)
        else:
            st.error("❌ No valid templates produced reactants for this molecule.")

    except Exception as e:
        st.error(f"⚠️ Error: {e}")
