<img width="500" alt="logo" src="https://github.com/Flo-fllt/Projet_chem/blob/main/Images/logo.png?raw=true">

[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/ThomasCsson/MASSIVEChem)
[![Python](https://img.shields.io/badge/Python-FFD43B?style=for-the-badge&logo=python&logoColor=blue)](https://www.python.org/)
[![Jupyter](https://img.shields.io/badge/Jupyter-F37626.svg?&style=for-the-badge&logo=Jupyter&logoColor=purple)](https://jupyter.org/)

# -         Retro-chem       -

## Packages 

RetroChem, is a one-step retrosynthesis engine, which essentially allows for a given molecule a prediction of the reactants needed to produce it. Developped as part of a project for the practical Programming in Chemistry course at EPFL (2025). 

Contributors:
- Giulio Matteo Garotti, second year chemical engineer at EPFL [![GitHub](https://img.shields.io/badge/GitHub-Giulio--grt-181717.svg?style=flat&logo=github)](http://github.com/Giulio-grt)
- Florian Follet, second year chemist at EPFL [![GitHub](https://img.shields.io/badge/GitHub-Flo--fllt-181717.svg?style=flat&logo=github)](https://github.com/Flo-fllt)
- Noah Paganzzi, second year chemical engineer at EPFL [![GitHub](https://img.shields.io/badge/GitHub-Noah--Paga-181717.svg?style=flat&logo=github)](https://github.com/Noah-Paga)
- Jacques Maurice Grandjean, second year chemical engineer at EPFL [![GitHub](https://img.shields.io/badge/GitHub-JacquesGrandjean-181717.svg?style=flat&logo=github)](https://github.com/JacquesGrandjean)

  import requests

# Remplace avec ton nom d'utilisateur et le repo
owner = "Flo-fllt"
repo = "Projet_chem"

# Si repo privé, ajoute un token GitHub perso avec les droits repo
headers = {
    "Authorization": "token ghp_0epWU6m791ql7WoxxMLg1jnEsOLvmtK4aUzWn"
}

url = f"https://api.github.com/repos/{owner}/{repo}/contributors"

response = requests.get(url, headers=headers)

if response.status_code == 200:
    contributors = response.json()
    for contributor in contributors:
        print(f"{contributor['login']} ({contributor['contributions']} contributions) - {contributor['html_url']}")
else:
    print("Erreur d'accès :", response.status_code)

## 📊 GitHub Stats

- 👥 [Contributors](https://github.com/Flo-fllt/Projet_chem/graphs/contributors)
- 📈 [Commit Activity](https://github.com/Flo-fllt/Projet_chem/graphs/commit-activity)
- 📊 [Code Frequency](https://github.com/Flo-fllt/Projet_chem/graphs/code-frequency)
- 🔁 [Pull Requests](https://github.com/Flo-fllt/Projet_chem/pulls)
