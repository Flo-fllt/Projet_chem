import streamlit as st
from streamlit_ketcher import st_ketcher


with st.expander("Draw Molecule"):
    molecule = st.text_input("**Smiles**") 
    ketcher_smiles= st_ketcher(molecule,height=600)



