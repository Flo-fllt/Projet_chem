from setuptools import setup, find_packages

setup(
    name="RetroChem",
    version="0.1.0",
    author="Jacques Grandjean, Noah Paganuzzi, Florian Follet, Giulio Garotti",
    description="Outil de rétrosynthèse, visualisation et machine learning",
    packages=find_packages(include=["Package_retrosynth", "Package_retrosynth.*"]),
    py_modules=["retrochem_launcher"],  # pour la commande CLI
    include_package_data=True,          # pour inclure les fichiers non-Python
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
    entry_points={
        "console_scripts": [
            "retrochem=retrochem_launcher:main"
        ]
    },
    python_requires=">=3.10",
)
