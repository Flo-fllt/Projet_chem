# Retrosynth

**Retrosynth** is a lightweight Python package designed for retrosynthetic analysis based on SMILES inputs.  
It includes tools for model training, template preprocessing, and visualization of reaction predictions.


## Features

- Clean SMILES strings by removing atom mappings
- Train machine learning models (e.g. MLPClassifier)
- Load and preprocess reaction template datasets
- Visualize top predicted templates (Streamlit-compatible)


## Installation

Clone this repository and install the package locally:

```bash or terminal
#Clone the repository
git clone https://github.com/Flo-fllt/Projet_chem.git

#naviguate to the package folder
cd Projet_chem/Package_retrosynth

#Install the package locally in editable mode
pip install -e .

#This will install the retrosynth package on your machine. You can now import and use its functions anywhere in your Python environment:
from retrosynthese.training import train_model
from retrosynthese.affichage import remove_atom_mapping

