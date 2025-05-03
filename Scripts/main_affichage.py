import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from io import BytesIO
from Retrosynthese.affichage import (mol_to_high_quality_image,st_scaled_image,apply_template,predict_topk_templates,render_reaction_scheme)


# --- Streamlit UI ---
with st.sidebar:
    st.image("images/logo.png", width=1000)

st.title("RetroSynthesis Prediction Tool")

# Molecule input section
with st.expander("1. Draw Molecule"):
    molecule = st.text_input("**Paste SMILES (optional)**")
    ketcher_smiles = st_ketcher(molecule, height=600)

final_smiles = ketcher_smiles or molecule
if final_smiles:
    st.success(f"✅ SMILES: {final_smiles}")

# Retrosynthesis logic
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

        total_prob = sum(prob for _, _, prob, _ in successful_predictions)
        normalized_predictions = [
            (template_hash, smarts, prob / total_prob, reactants)
            for (template_hash, smarts, prob, reactants) in successful_predictions
        ] if total_prob > 0 else []

        if normalized_predictions:
            st.markdown("### Retrosynthesis Predictions")
            for idx, (template_hash, smarts, norm_prob, reactants) in enumerate(normalized_predictions, 1):
                with st.expander(f" Prediction {idx} - {norm_prob * 100:.2f}% confidence"):

                    # Display first-level reactants
                    st.markdown("**Reactants:**")
                    cols = st.columns(len(reactants))
                    for i, smi in enumerate(reactants):
                        mol = Chem.MolFromSmiles(smi)
                        if mol:
                            img = mol_to_high_quality_image(mol)
                            with cols[i]:
                                st_scaled_image(img, width_display_px=300)

                    # Try second-level retrosynthesis if there's only one reactant
                    step2_reactants = None
                    if len(reactants) == 1:
                        st.markdown("**↪ Step 2 - Retrosynthesis:**")
                        second_predictions = predict_topk_templates(reactants[0], topk=50)
                        for t_hash, smarts2, p2 in second_predictions:
                            reactant_products = apply_template(smarts2, reactants[0])
                            if reactant_products:
                                step2_reactants = reactant_products[0]
                                subcols = st.columns(len(step2_reactants))
                                for j, smi2 in enumerate(step2_reactants):
                                    mol2 = Chem.MolFromSmiles(smi2)
                                    if mol2:
                                        img2 = mol_to_high_quality_image(mol2)
                                        with subcols[j]:
                                            st_scaled_image(img2, width_display_px=300)
                                break
                        else:
                            st.markdown("- No further retrosynthesis found.")

                    # --- Reaction scheme with images ---
                    if step2_reactants:
                        smiles_chain = [step2_reactants, reactants, [final_smiles]]
                    else:
                        smiles_chain = [reactants, [final_smiles]]

                    scheme_html = render_reaction_scheme(smiles_chain)
                    st.markdown("**Reaction Scheme:**", unsafe_allow_html=True)
                    st.markdown(scheme_html, unsafe_allow_html=True)

        else:
            st.error("❌ No valid templates produced any reactants.")
    except Exception as e:
        st.error(f"⚠️ Error: {e}")



