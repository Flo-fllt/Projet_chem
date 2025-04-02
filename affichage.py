import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors
from stmol import showmol
import py3Dmol
from pathlib import Path
import pandas as pd
import os
from streamlit_ketcher import st_ketcher
from rdkit.Chem import rdFingerprintGenerator
import numpy as np
import streamlit.components.v1 as components
from typing import Tuple, List


def visualize_3D(molstring):
    "Visualize the molecule in 3D using stmol"
    w, h = 400, 400
    xyzview = py3Dmol.view(width=w,height=w)
    xyzview.addModel(molstring,'mol')
    xyzview.setStyle({'sphere':{'colorscheme':'cyanCarbon', 'scale':0.25}, 'stick':{'colorscheme':'cyanCarbon'}})
    xyzview.zoomTo()
    xyzview.spin()
    xyzview.setBackgroundColor('white')
    showmol(xyzview, height = w,width=w)

def download_data():
    "Download the ChEMBL database"
    csv_file = Path(r"C:\Users\flopi\python-course\git\ppchem\chembl_drugs.csv")# modifier la source
    df = pd.read_csv(csv_file, sep= ";")
    df.head()
    return df


with st.expander("Draw Molecule"):
    molecule = st.text_input("**Smiles**") 
    ketcher_smiles= st_ketcher(molecule,height=600)

df = download_data()



