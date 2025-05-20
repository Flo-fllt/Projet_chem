<img width="500" alt="logo" src="https://github.com/Flo-fllt/Projet_chem/blob/main/Images/logo.png?raw=true">

# -         RetroChem       -                                                                                                                                                                                   

## 📜 Package information
RetroChem, is a one-step retrosynthesis engine, it is also a pip installable python package. It includes tools for model training, template preprocessing, and visualization of reaction predictions. Developped as part of a project for the practical Programming in Chemistry at EPFL (2025). 

[![EPFL Course](https://img.shields.io/badge/EPFL-red?style=for-the-badge)](https://edu.epfl.ch/coursebook/en/practical-programming-in-chemistry-CH-200)

### 🪄 Features

- Simple to use retrosynthesis engine. 
- Draw molecules or use molecule name (eg. paracetamol).
- Gives clean and chemically correct molecules.
- Works well for known molecules (eg. GHB, EDTA and etc).
- Predicts possible reactants and reaction steps.
- Gives the predictions in order of confidence

### 👥 Contributors
- Giulio Matteo Garotti, second year chemical engineer at EPFL           [![GitHub](https://img.shields.io/badge/GitHub-Giulio--grt-181717.svg?style=flat&logo=github)](http://github.com/Giulio-grt)
- Florian Follet, second year chemist at EPFL                            [![GitHub](https://img.shields.io/badge/GitHub-Flo--fllt-181717.svg?style=flat&logo=github)](https://github.com/Flo-fllt)
- Noah Paganzzi, second year chemical engineer at EPFL                   [![GitHub](https://img.shields.io/badge/GitHub-Noah--Paga-181717.svg?style=flat&logo=github)](https://github.com/Noah-Paga)
- Jacques Maurice Grandjean, second year chemical engineer at EPFL       [![GitHub](https://img.shields.io/badge/GitHub-JacquesGrandjean-181717.svg?style=flat&logo=github)](https://github.com/JacquesGrandjean)

![Maintained](https://img.shields.io/badge/Maintained-Yes-green) ![Contributors](https://img.shields.io/badge/Contributors-4-blue) [![python3.10](https://img.shields.io/badge/Python-3.10-3776AB.svg?style=flat&logo=python&logoColor=orange)](https://www.python.org) [![LICENSE](https://img.shields.io/badge/License-MIT-purple.svg)](https://github.com/Flo-fllt/Projet_chem/blob/main/LICENSE.txt) ![coverage](https://img.shields.io/badge/coverage-96%25-brightgreen)


View contributors : [Contributors on GitHub](https://github.com/Flo-fllt/Projet_chem/graphs/contributors)

View commit activity: [Commit activity on GitHub](https://github.com/Flo-fllt/Projet_chem/graphs/commit-activity)

View code frequency : [Code frequency on GitHub](https://github.com/Flo-fllt/Projet_chem/graphs/code-frequency)

[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/ThomasCsson/MASSIVEChem)
[![Python](https://img.shields.io/badge/Python-FFD43B?style=for-the-badge&logo=python&logoColor=blue)](https://www.python.org/)
[![Jupyter](https://img.shields.io/badge/Jupyter-F37626.svg?&style=for-the-badge&logo=Jupyter&logoColor=purple)](https://jupyter.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-FF4B4B.svg?&style=for-the-badge&logo=Streamlit&logoColor=white)](https://streamlit.io/) 
[![Anaconda](https://img.shields.io/badge/Anaconda-44A833.svg?&style=for-the-badge&logo=anaconda&logoColor=white)](https://www.anaconda.com/)

[![Made_with_python](http://ForTheBadge.com/images/badges/made-with-python.svg)](https://www.python.org/)


### 🔍 Retrosynthesis, what is it?
Organic retrosynthesis is a problem-solving technique used in organic chemistry to design a synthetic route for a target molecule by breaking it down into simpler precursor structures. This process, known as retrosynthetic analysis, helps chemists plan the step-by-step synthesis of complex organic compounds by working backward from the desired product.

#### The approach involves two main stages:

Disconnection: The target molecule is “disconnected” at strategic bonds to identify simpler molecules that could be used as starting materials. These disconnections are guided by known chemical reactions and functional group transformations. Each disconnection leads to one or more synthons, which are idealised fragments representing the functional components of the molecule. For example, a carbon–carbon bond may be disconnected to yield a nucleophile and an electrophile.

Reagent Identification: Once the synthons are defined, they are translated into real, purchasable or synthesizable compounds known as synthetic equivalents. These are then used to build the molecule in the forward direction. This step requires the selection of appropriate reagents, protecting groups, and reaction conditions to ensure a feasible and efficient synthesis.

A retrosynthetic tree is often constructed to visualise all possible synthetic pathways, with the target molecule at the top and potential starting materials branching below. Chemists use both strategic bonds and known reaction patterns (like aldol condensations, Grignard additions, or Diels-Alder reactions) to inform these choices.

- Retrosynthesis is widely applied in:
- Pharmaceutical chemistry, for the development of active pharmaceutical ingredients (APIs).
- Natural product synthesis, for recreating complex bioactive compounds from simpler building blocks.
- Material science, for designing functional molecules with specific properties.
- Green chemistry, to plan more sustainable and efficient synthetic routes.

The logic of retrosynthesis also lends itself well to computational chemistry. Algorithms and machine learning tools are increasingly used to automate retrosynthetic analysis, generating synthetic routes based on reaction databases and predictive models.

The theoretical underpinning of retrosynthesis often involves mapping the retrosynthetic steps to known reaction mechanisms and using heuristics based on reactivity and selectivity. The goal is to find the shortest, most cost-effective, and most reliable pathway to the target compound.

Let us now walk through how to apply retrosynthetic analysis to a real molecule using this framework!
## 🕹️ How to install it

Firstly it is advised to create a CONDA environment:
```bash or terminal
#Open bash or terminal
#Name the environment as you wish and specify python 3.10
conda create -n env.name python==3.10

#Activate your environment
conda activate env.name

#Alternative way of activating it
source activate env.name
````
Clone this repository and install the python package locally:

```bash or terminal
#Clone the repository
git clone https://github.com/Flo-fllt/RetroChem.git

#Naviguate to the RetroChem folder
cd RetroChem

#Install the package locally in editable mode, make sure to activate your environment before doing so
pip install retrochem

#This will install the retrosynth package on your machine. You can now run the program with:
retrochem
````
The streamlit web page for RetroChem will be opened on your default browser.
## 💡 Requirements
````
python>=3.10
streamlit
pandas
numpy
joblib
rdkit
scikit-learn
matplotlib
pillow
streamlit-ketcher
````

## 💻 Guide
Once the RetroChem web page is open the retrosynthesis engine is ready to run!

<img width="800" alt="interface guide 1" src="https://github.com/Flo-fllt/RetroChem/blob/main/Images/Interface%20guide%201.png?raw=true">

If the program has a difficult time finding reactants, try another molecule and especially known molecules.

<img width="800" alt="interface guide 2" src="https://github.com/Flo-fllt/RetroChem/blob/main/Images/Interface%20guide%202.png?raw=true">
 
<img width="800" alt="interface guide 3" src="https://github.com/Flo-fllt/RetroChem/blob/main/Images/Interface%20guide%203.png?raw=true">

<img width="800" alt="interface guide 4" src="https://github.com/Flo-fllt/RetroChem/blob/main/Images/Interface%20guide%204.png?raw=true">

## 🛟 Need help?

If for some reason the program does not work or an issue occurs during the search of retrosynthesis templates there are a few steps you can take to ensure everything is in order.

First make sure that you are in the correct environment where you installed the package: 
```` bash or terminal
#Show what current environnement you are in
which python
````

If not: 
```` bash or terminal
#Activate the environnement where you downloaded the retrochem package
conda activate env.name
````

Then navigate to the top-level RetroChem repository folder: 
```` bash or terminal
cd retrochem

# Confirm your current location
pwd

# It should end with: /your/path/to/retrochem
````
❗️ Warning: If your path ends in /retrochem/retrochem, you are one level too deep, go back using:
```` bash or terminal
cd ..
````

Thirdly check that you have the latest version of RetroChem
```` bash or terminal
#Shows you which version of retrochem you are using
pip show retrochem

#Updates the package
pip install --upgrade retrochem

#If the issue persists try uninstalling the package, then installing it again with the latest version specified
pip uninstall retrochem
pip install retrochem==x.x.x #change the x.x.x to the current version (check PyPi page below)
````
You can compare the version you have downloaded to that of the newest available version, which can be found on the PiPy page: 

[![RetroChem.PYPI](https://img.shields.io/badge/pypi-3775A9?style=for-the-badge&logo=pypi&logoColor=white)](https://pypi.org/project/RetroChem/)

If for some reason the issue is not resolved, it may derive from a pip issue, in this case checking for a pip update may solve it:
````
#For a virtual environment
pip install --upgrade pip

#For Linux/MacOS
python3 -m pip install --upgrade pip
````
For questions on how to set-up/download the necesseties (ANACONDA and etc) we welcome you to take a look at this:

[![GitHub](https://img.shields.io/badge/GitHub-PPCHEM-black?logo=github)](https://github.com/schwallergroup/practical-programming-in-chemistry-exercises/tree/main/Lecture01)

## 🤖 How to retrain our model

To retrain the retrosynthesis prediction model, use the notebook located in:
````
retrain_model/Train_model.ipynb
````
By default, this notebook is configured to work with the original USPTO-50K dataset.

1. Retrain with the default dataset
   - Simply run the two cell blocks in the notebook:
   - The first one combines the three original datasets (train, valid, and test) into one file.
   - The second cell block loads the combined file, retrains the model, and saves:
````
retrain_model/new_mlp_classifier_model.pkl
retrain_model/new_scaler.pkl
retrain_model/new_label_encoder.pkl
````
2. Retrain with your own data (not recommended)
   - Skip the first cell block
   - Update the file path in the second cell block:
````  
# Replace this path with your own
combined_file_path = os.path.join(data_dir, "your_data.csv")
````
   - Ensure your CSV includes the following required columns:
````
RxnSmilesClean, PseudoHash, RetroTemplate, TemplateHash
````
Disclaimer: Retraining with custom data is not recommended unless you are familiar with the dataset structure and preprocessing requirements. The model was trained with curated USPTO-50K data, and using inconsistent or improperly formatted data may lead to poor performance or errors.

## 📚 Want more information?

Here is the link to the report of the project, logging detailed information regarding the model we used and how it was trained along with the general context of this project.

[![Jupyter Notebook](https://img.shields.io/badge/Jupyter_Notebook-orange.svg)](https://github.com/Flo-fllt/RetroChem/blob/main/Notebook/project_report.ipynb)

## 🕺Your turn!

With all this information you should be more than ready to give RetroChem a try and don't forget to let us know what you think!
