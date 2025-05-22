from setuptools import setup, find_packages
from pathlib import Path

# Lire le contenu de ton README
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

setup(
    name="RetroChem",
    version="0.2.0",   
    author="Jacques Grandjean, Noah Paganuzzi, Florian Follet, Giulio Garotti",
    description="Retrosynthesis, visualization and machine learning tools",
    long_description=long_description,
    long_description_content_type="text/markdown", 
    packages=find_packages(include=["RetroChem", "Retrochem.*"]),
    py_modules=["retrochem_launcher"],  # pour la commande CLI
    include_package_data=True,          # pour inclure les fichiers non-Python
    install_requires=[
        "streamlit",
        "pandas",
        "numpy==2.2.6",
        "joblib",
        "rdkit",
        "scikit-learn",
        "matplotlib",
        "pillow",
        "streamlit-ketcher"
    ],
    entry_points={
        "console_scripts": [
            "retrochem=retrochem_launcher:main"
        ]
    },
    python_requires=">=3.10",
)
