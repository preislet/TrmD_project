import os
import requests
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser(description="Download a single molecule from ChEMBL by its ChEMBL ID.")
parser.add_argument("--chembl_id", type=str, required=True, help="The ChEMBL ID of the molecule to download.")
parser.add_argument("--output_format", type=str, default="sdf", choices=["sdf", "mol", "smiles"], help="The format to save the molecule.")
parser.add_argument("--output_folder", type=str, default="downloaded_molecule", help="Folder to save the downloaded molecule.")
parser.add_argument("--preprocess", action="store_true", default=False, help="Whether to preprocess the molecule using RDKit.")

def download_molecule(chembl_id: str, output_format: str, output_folder: str, preprocess: bool) -> None:
    """
    Downloads a single molecule from ChEMBL and optionally preprocesses it.

    Args:
        chembl_id (str): The ChEMBL ID of the molecule.
        output_format (str): The format to save the molecule (sdf, mol, or smiles).
        output_folder (str): Folder to save the downloaded molecule.
        preprocess (bool): Whether to preprocess the molecule using RDKit.
    """
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Construct the URL for downloading the molecule
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.{output_format}"
    file_path = os.path.join(output_folder, f"{chembl_id}.{output_format}")
    
    print(f"Downloading molecule {chembl_id} in {output_format} format...")

    # Download the molecule
    response = requests.get(url)
    if response.status_code == 200:
        with open(file_path, "wb") as file:
            file.write(response.content)
        print(f"Downloaded molecule saved to: {file_path}")

        # Preprocess if required
        if preprocess and output_format in ["mol", "sdf"]:
            preprocess_molecule(file_path)
    else:
        print(f"Failed to download molecule {chembl_id}. Status Code: {response.status_code}")

def preprocess_molecule(file_path: str) -> None:
    """
    Preprocesses a molecule file using RDKit (e.g., conformer generation).

    Args:
        file_path (str): Path to the molecule file to preprocess.
    """
    print(f"Preprocessing molecule at {file_path}...")
    mol = Chem.MolFromMolFile(file_path)
    if mol:
        AllChem.EmbedMolecule(mol)
        processed_path = file_path.replace(".mol", "_preprocessed.mol").replace(".sdf", "_preprocessed.sdf")
        Chem.MolToMolFile(mol, processed_path)
        print(f"Preprocessed molecule saved to: {processed_path}")
    else:
        print(f"Failed to preprocess molecule at {file_path}. RDKit could not parse the file.")

if __name__ == "__main__":
    args = parser.parse_args()
    download_molecule(args.chembl_id, args.output_format, args.output_folder, args.preprocess)
