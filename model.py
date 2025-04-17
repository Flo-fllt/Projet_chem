import joblib
import pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import re

# Chargement des données
file_path = "/Users/giuliogarotti/Documents/GitHub/Projet_chem/uspto50/uspto50/reaction_templates_50k_test.csv"
df = pd.read_csv(file_path, sep="\t")

# Sélectionner les colonnes nécessaires
colonnes_necessaires = ['RxnSmilesClean', 'PseudoHash', 'RetroTemplate', 'TemplateHash']
df = df[colonnes_necessaires]
df_unique = df.drop_duplicates(subset=['PseudoHash'], keep='first')

# Fonction pour corriger les erreurs communes dans les SMILES
def fix_common_errors(smiles):
    smiles = re.sub(r'@', 'O', smiles)
    return smiles

# Fonction pour séparer les réactifs et les produits à partir de la réaction SMILES
def split_rxn_smiles(rxn_smiles):
    try:
        rxn_smiles = fix_common_errors(rxn_smiles) 
        parts = rxn_smiles.split(">>")  # Séparer réactifs et produits
        if len(parts) != 2:
            return []  # Retourner une liste vide si le format est incorrect

        reactants_raw, products_raw = parts
        reactants = reactants_raw.split(".")  # Séparer les réactifs s'il y en a plusieurs
        products = products_raw.split(".")  # Séparer les produits s'il y en a plusieurs

        return reactants, products

    except Exception as e:
        print(f"Erreur dans split_rxn_smiles: {rxn_smiles} -> {e}")
        return []

# Fonction pour convertir les SMILES en empreintes
def smiles_to_fingerprints(smiles):
    try:
        reactants_smiles, products_smiles = split_rxn_smiles(smiles)  # Séparer réactifs et produits

        # Fonction pour générer les empreintes à partir des SMILES
        def mols_to_fps(smiles_list):
            fps = []
            n_bits = 2048  # Taille de l'empreinte
            for s in smiles_list:
                mol = Chem.MolFromSmiles(s)  # Convertir le SMILES en molécule
                if mol is None:
                    print(f"Molécule invalide : {s}")
                    fps.append(np.zeros(n_bits))  # Retourner un vecteur de zéros pour les erreurs
                else:
                    fp = AllChem.GetMorganFingerprint(mol, radius=3)
                    arr = np.zeros((n_bits,), dtype=int)
                    if isinstance(fp, Chem.DataStructs.UIntSparseIntVect):
                        on_bits = list(fp.GetNonzeroElements().keys())
                        for bit in on_bits:
                            arr[bit % n_bits] = 1
                    fps.append(arr)
            return fps

        # Générer les empreintes pour les réactifs et produits
        reactants_fps = mols_to_fps(reactants_smiles)
        products_fps = mols_to_fps(products_smiles)

        max_len = 2048
        reactants_fps = [np.pad(fp, (0, max_len - len(fp))) if len(fp) < max_len else fp for fp in reactants_fps]
        products_fps = [np.pad(fp, (0, max_len - len(fp))) if len(fp) < max_len else fp for fp in products_fps]
        
        return reactants_fps, products_fps

    except Exception as e:
        print(f"Erreur lors du parsing de SMILES: {smiles}, erreur: {e}")
        return [], []  # Retourner des listes vides en cas d'erreur
    
def prepare_fingerprints_for_training(df):
    X = []
    y = []

    # Traiter chaque ligne du DataFrame
    for smiles, target in zip(df['RxnSmilesClean'], df['TemplateHash']):
        reactants_fps, products_fps = smiles_to_fingerprints(smiles)
        
        # Vérifier si les réactifs et produits ont été correctement générés
        if reactants_fps and products_fps:
            # Ajouter les empreintes des réactifs à X et dupliquer la cible pour chaque réactif
            X.extend(reactants_fps)
            y.extend([target] * len(reactants_fps))  # Associer la même cible pour chaque réactif
            
            # Ajouter les empreintes des produits à X et dupliquer la cible pour chaque produit
            X.extend(products_fps)
            y.extend([target] * len(products_fps))  # Associer la même cible pour chaque produit

    # Conversion en tableau NumPy
    X = np.array(X)
    y = np.array(y)

    # Vérification de la correspondance des tailles
    if X.shape[0] != y.shape[0]:
        raise ValueError(f"Le nombre d'exemples dans X ({X.shape[0]}) ne correspond pas à celui de y ({y.shape[0]}).")

    return X, y


# Utilisation de la fonction de préparation des empreintes
X, y = prepare_fingerprints_for_training(df)

# Entraîner ton modèle
from sklearn.neural_network import MLPClassifier

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

model = MLPClassifier(hidden_layer_sizes=(64, 64), max_iter=1000, random_state=42)
model.fit(X, y)  # Entraîner le modèle

print("Modèle entraîné avec succès")

# Sauvegarder le modèle et le scaler
joblib.dump(model, 'mlp_classifier_model.pkl')
joblib.dump(scaler, 'scaler.pkl')
print("Modèle sauvegardé dans 'mlp_classifier_model.pkl'")
print("Scaler sauvegardé dans 'scaler.pkl'")