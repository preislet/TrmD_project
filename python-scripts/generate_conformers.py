from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import PDBWriter
import os
import subprocess


def generate_protonation_states(smiles, pH=7.4):
    """
    Generate protonation states for a given SMILES string using Open Babel.

    Parameters:
        smiles (str): SMILES representation of the molecule.
        pH (float): Target pH for protonation.

    Returns:
        list: List of SMILES strings representing different protonation states.
    """
    try:
        # Write the input SMILES to a temporary file
        with open("input.smi", "w") as f:
            f.write(smiles)

        # Use Open Babel to generate protonation states
        command = [
            "obabel",
            "input.smi",
            "-osmi",
            "--protonate",
            f"--pH={pH}",
            "-O",
            "output.smi",
        ]
        subprocess.run(command, check=True)

        # Read the output SMILES
        with open("output.smi", "r") as f:
            protonated_smiles = [line.split()[0] for line in f]

        # Clean up temporary files
        os.remove("input.smi")
        os.remove("output.smi")

        return protonated_smiles

    except Exception as e:
        print(f"Error during protonation: {e}")
        return [smiles]  # Return original SMILES if error occurs

def generate_conformers(smiles, output_dir, num_conformers=5):
    """
    Generuje 3D konformery pro zadanou molekulu a ukládá je do PDB formátu.

    Parameters:
        smiles (str): SMILES reprezentace molekuly.
        output_dir (str): Složka, kam se uloží PDB soubory.
        num_conformers (int): Počet 3D konformerů, které se mají generovat.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Načti molekulu z SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Neplatný SMILES řetězec: {smiles}")
        return
    
    # Přidej vodíky
    mol = Chem.AddHs(mol)
    
    # Generuj 3D konformery
    AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, randomSeed=42)
    
    # Ulož každý konformer jako PDB
    for conf_id in range(mol.GetNumConformers()):
        pdb_file = os.path.join(output_dir, f"mol_conf_{conf_id}.pdb")
        with PDBWriter(pdb_file) as writer:
            writer.write(mol, confId=conf_id)
    
    print(f"Všechny konformery byly uloženy do složky: {output_dir}")

# Použití
if __name__ == "__main__":
    # Vstupní SMILES pro molekulu
    smiles = "CCO"  # Etanol (příklad)
    
    # Výstupní složka
    output_dir = "output_conformers"
    
    # Spusť generování
    generate_conformers(smiles, output_dir, num_conformers=5)

