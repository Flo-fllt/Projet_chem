import streamlit as st
from streamlit_ketcher import st_ketcher

# Define the function at the top level
def save_smiles_to_file(smiles, filename="smiles.txt"):
    with open(filename, "w") as f:
        f.write(smiles)

with st.expander("Draw Molecule"):
    molecule = st.text_input("**Smiles**") 
    ketcher_smiles = st_ketcher(molecule, height=600)
    
    # Add a button to convert to SMILES
    if st.button("Convert to SMILES"):
        if ketcher_smiles:
            save_smiles_to_file(ketcher_smiles)
            st.success(f"SMILES saved: {ketcher_smiles}")
        else:
            st.warning("Please draw a molecule first!")