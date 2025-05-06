<img width="500" alt="logo" src="https://github.com/Flo-fllt/Projet_chem/blob/main/Images/logo.png?raw=true">

[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/ThomasCsson/MASSIVEChem)
[![Python](https://img.shields.io/badge/Python-FFD43B?style=for-the-badge&logo=python&logoColor=blue)](https://www.python.org/)
[![Jupyter](https://img.shields.io/badge/Jupyter-F37626.svg?&style=for-the-badge&logo=Jupyter&logoColor=purple)](https://jupyter.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-FF4B4B.svg?&style=for-the-badge&logo=Streamlit&logoColor=white)](https://streamlit.io/) 
# -         RetroChem       -                                                                                                                                                                                   
#### Poject in practical programming in chemistry

## 📦Package information
RetroChem, is a one-step retrosynthesis engine, which essentially allows for a given molecule a prediction of the reactants needed to produce it. Developped as part of a project for the practical Programming in Chemistry course at EPFL (2025). 
### 👥 Contributors
- Giulio Matteo Garotti, second year chemical engineer at EPFL [![GitHub](https://img.shields.io/badge/GitHub-Giulio--grt-181717.svg?style=flat&logo=github)](http://github.com/Giulio-grt)
- Florian Follet, second year chemist at EPFL [![GitHub](https://img.shields.io/badge/GitHub-Flo--fllt-181717.svg?style=flat&logo=github)](https://github.com/Flo-fllt)
- Noah Paganzzi, second year chemical engineer at EPFL [![GitHub](https://img.shields.io/badge/GitHub-Noah--Paga-181717.svg?style=flat&logo=github)](https://github.com/Noah-Paga)
- Jacques Maurice Grandjean, second year chemical engineer at EPFL [![GitHub](https://img.shields.io/badge/GitHub-JacquesGrandjean-181717.svg?style=flat&logo=github)](https://github.com/JacquesGrandjean)

![Maintained](https://img.shields.io/badge/Maintained-Yes-green) ![Contributors](https://img.shields.io/badge/Contributors-4-blue) [![python3.10](https://img.shields.io/badge/Python-3.10-3776AB.svg?style=flat&logo=python&logoColor=orange)](https://www.python.org) [![LICENSE](https://img.shields.io/badge/License-MIT-purple.svg)](https://github.com/Flo-fllt/Projet_chem/blob/main/LICENSE.txt)

View contributors : [Contributors on GitHub](https://github.com/Flo-fllt/Projet_chem/graphs/contributors)

View commit activity: [Commit activity on GitHub](https://github.com/Flo-fllt/Projet_chem/graphs/commit-activity)

View code frequency : [Code frequency on GitHub](https://github.com/Flo-fllt/Projet_chem/graphs/code-frequency)
### 🧪 Retrosynthesis, what is it?
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
## 🛠️ How to install it
