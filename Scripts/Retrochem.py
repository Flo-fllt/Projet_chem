import streamlit as st
import requests
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from Package_retrosynth.Package_functions.Interface_functions import (mol_to_high_quality_image, st_scaled_image, apply_template, predict_topk_templates, render_reaction_scheme)

# --- Sidebar ---
with st.sidebar:
    st.image("images/logo.png", width=1000)

st.title("RetroSynthesis Prediction Tool")

# --- Molecule input section ---
with st.expander("1. Input Molecule", expanded=True):
    
    # Search by name (from PubChem)
    mol_name = st.text_input("Search molecule by name", "")
    name_smiles = ""
    if mol_name:
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{mol_name}/property/CanonicalSMILES/JSON"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                name_smiles = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
                st.success(f"Found SMILES: {name_smiles}")
            else:
                st.error("❌ Molecule not found in PubChem library.")
        except Exception as e:
            st.error(f"Error fetching SMILES: {e}")

    # Draw or edit the molecule (using SMILES from name or blank)
    ketcher_smiles = st_ketcher(name_smiles, height=600)

# Final SMILES from Ketcher
final_smiles = ketcher_smiles
if final_smiles:
    st.success(f"✅ SMILES: {final_smiles}")


# --- Retrosynthesis logic ---
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
                    st.markdown("**Reactants:**")
                    cols = st.columns(len(reactants))
                    for i, smi in enumerate(reactants):
                        mol = Chem.MolFromSmiles(smi)
                        if mol:
                            img = mol_to_high_quality_image(mol)
                            with cols[i]:
                                st_scaled_image(img, width_display_px=300)

                    # Step 2 retrosynthesis 
                    step2_reactants = None
                    if len(reactants) == 1:
                        st.markdown("**↪ Step 2 - Retrosynthesis:**")
                        second_predictions = predict_topk_templates(reactants[0], topk=50)
                        for t_hash, smarts2, p2 in second_predictions:
                            reactant_products = apply_template(smarts2, reactants[0])
                            if not reactant_products or not reactant_products[0]:
                                continue

                            candidate_reactants = reactant_products[0]

                            # Canonicalize comparison to avoid format mismatch
                            canonical_final = Chem.MolToSmiles(Chem.MolFromSmiles(final_smiles))
                            canonical_candidates = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in candidate_reactants]

                            if canonical_final in canonical_candidates:
                                continue  

                            step2_reactants = candidate_reactants
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

                    # Reaction Scheme
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