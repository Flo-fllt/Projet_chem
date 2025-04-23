from setuptools import setup, find_packages

setup(
    name="retrosynth",
    version="0.1.0",
    author="Ton Nom",
    description="Package pour la rétrosynthèse, visualisation et entraînement de modèles",
    packages=find_packages(),
    install_requires=[
        "streamlit",
        "pandas",
        "numpy",
        "joblib",
        "rdkit",
        "scikit-learn",
        "matplotlib",
        "pillow", 
        "streamlit-ketcher"
    ],
    python_requires=">=3.10"
)
