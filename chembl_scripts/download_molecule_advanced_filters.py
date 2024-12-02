import os
import requests
import argparse
import pandas as pd
from typing import List, Dict
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem
from rdkit.Chem import AllChem

# Argument parser setup
parser = argparse.ArgumentParser(description="Download molecules from ChEMBL with advanced filtering and preprocessing.")
parser.add_argument("--max_results", type=int, default=1000, help="Maximum number of molecules to download.")
parser.add_argument("--output_format", type=str, default="sdf", choices=["sdf", "mol", "smiles"], help="Output format for molecule structures.")
parser.add_argument("--processors", type=int, default=4, help="Number of processors to use for downloading.")
parser.add_argument("--folder_name", type=str, default="molecule_structures", help="Base folder to store molecule files.")
parser.add_argument("--max_phase", type=int, default=4, help="Maximum phase for filtering molecules (e.g., 4 for FDA-approved drugs).")
parser.add_argument("--min_molecular_weight", type=float, help="Minimum molecular weight for filtering.")
parser.add_argument("--max_molecular_weight", type=float, help="Maximum molecular weight for filtering.")
parser.add_argument("--min_logp", type=float, help="Minimum logP for filtering.")
parser.add_argument("--max_logp", type=float, help="Maximum logP for filtering.")
parser.add_argument("--min_hba", type=int, help="Minimum hydrogen bond acceptors.")
parser.add_argument("--max_hbd", type=int, help="Maximum hydrogen bond donors.")
parser.add_argument("--max_rotatable_bonds", type=int, help="Maximum number of rotatable bonds.")
parser.add_argument("--min_psa", type=float, help="Minimum polar surface area for filtering.")
parser.add_argument("--max_psa", type=float, help="Maximum polar surface area for filtering.")
parser.add_argument("--max_lipinski_violations", type=int, help="Maximum Lipinski Rule of 5 violations allowed.")
parser.add_argument("--min_np_likeness", type=float, help="Minimum natural product likeness score.")
parser.add_argument("--max_np_likeness", type=float, help="Maximum natural product likeness score.")
parser.add_argument("--molecular_species", type=str, help="Filter by molecular species (e.g., NEUTRAL, ACID, BASE).")
parser.add_argument("--preprocess", action="store_true", default=False, help="Preprocess molecules using RDKit.")

def create_folders(base_folder: str) -> Dict[str, str]:
    """Creates folders for storing original and preprocessed files."""
    original_folder = os.path.join(base_folder, "original")
    preprocessed_folder = os.path.join(base_folder, "preprocessed")
    os.makedirs(original_folder, exist_ok=True)
    os.makedirs(preprocessed_folder, exist_ok=True)
    print(f"Created folders: {original_folder}, {preprocessed_folder}")
    return {"original": original_folder, "preprocessed": preprocessed_folder}

def preprocess_molecule(file_path: str, output_folder: str) -> None:
    """Preprocess a molecule file using RDKit."""
    mol = Chem.MolFromMolFile(file_path)
    if mol:
        AllChem.EmbedMolecule(mol)
        processed_path = os.path.join(output_folder, os.path.basename(file_path).replace(".mol", "_preprocessed.mol").replace(".sdf", "_preprocessed.sdf"))
        Chem.MolToMolFile(mol, processed_path)
        print(f"Preprocessed and saved to: {processed_path}")

def download_structure(chembl_id: str, folders: Dict[str, str], output_format: str, preprocess: bool) -> None:
    """Downloads a molecule structure in the specified format."""
    original_folder = folders["original"]
    preprocessed_folder = folders["preprocessed"]
    structure_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.{output_format}"
    file_path = os.path.join(original_folder, f"{chembl_id}.{output_format}")
    
    if os.path.exists(file_path):
        print(f"Skipping already downloaded: {chembl_id}")
        return
    
    response = requests.get(structure_url)
    if response.status_code == 200:
        with open(file_path, "wb") as file:
            file.write(response.content)
        print(f"Downloaded: {chembl_id}.{output_format}")
        if preprocess and output_format in ["mol", "sdf"]:
            preprocess_molecule(file_path, preprocessed_folder)
    else:
        print(f"Failed to download {chembl_id}: {response.status_code}")

def fetch_and_download_molecules(args: argparse.Namespace) -> None:
    """Fetches molecule IDs, applies filters, downloads molecules, and saves metadata."""
    base_url: str = "https://www.ebi.ac.uk/chembl/api/data/molecule"
    params: dict = {
        "format": "json",
        "limit": 100,
        "offset": 0,
        "max_phase": args.max_phase
    }

    folders: dict = create_folders(args.folder_name)
    metadata: List[Dict] = []
    chembl_ids: List[str] = []
    downloaded_count: int = 0

    # Fetch molecule IDs in batches
    while downloaded_count < args.max_results:
        response = requests.get(base_url, params=params)
        if response.status_code != 200:
            print(f"Failed to fetch molecule data: {response.status_code}")
            break

        data = response.json()
        molecules = data.get("molecules", [])
        if not molecules:
            break

        for molecule in molecules:
            chembl_id: str = molecule.get("molecule_chembl_id", "")
            molecule_properties = molecule.get("molecule_properties", {})

            molecular_weight: float = float(molecule_properties.get("full_mwt", 0))
            logp: float = float(molecule_properties.get("alogp", 0))
            hba: int = int(molecule_properties.get("hba", 0))
            hbd: int = int(molecule_properties.get("hbd", 0))
            rotatable_bonds: int = int(molecule_properties.get("rtb", 0))
            psa: float = float(molecule_properties.get("psa", 0))
            np_likeness: float = float(molecule_properties.get("np_likeness_score", 0))
            lipinski_violations: int = int(molecule_properties.get("num_lipinski_ro5_violations", 0))
            species: str = molecule_properties.get("molecular_species", "")

            # Apply filters
            if args.min_molecular_weight and molecular_weight < args.min_molecular_weight: continue
            if args.max_molecular_weight and molecular_weight > args.max_molecular_weight: continue
            if args.min_logp and logp < args.min_logp: continue
            if args.max_logp and logp > args.max_logp: continue
            if args.min_hba and hba < args.min_hba: continue
            if args.max_hbd and hbd > args.max_hbd: continue
            if args.max_rotatable_bonds and rotatable_bonds > args.max_rotatable_bonds:continue
            if args.min_psa and psa < args.min_psa: continue
            if args.max_psa and psa > args.max_psa: continue
            if args.max_lipinski_violations and lipinski_violations > args.max_lipinski_violations: continue
            if args.min_np_likeness and np_likeness < args.min_np_likeness: continue
            if args.max_np_likeness and np_likeness > args.max_np_likeness: continue
            if args.molecular_species and species.lower() != args.molecular_species.lower(): continue

            chembl_ids.append(chembl_id)
            metadata.append({
                "chembl_id": chembl_id,
                "molecular_weight": molecular_weight,
                "logp": logp,
                "hba": hba,
                "hbd": hbd,
                "rotatable_bonds": rotatable_bonds,
                "psa": psa,
                "np_likeness": np_likeness,
                "lipinski_violations": lipinski_violations,
                "species": species
            })

            downloaded_count += 1
            if downloaded_count >= args.max_results:
                break

        params["offset"] += len(molecules)

    print(f"Fetched {len(chembl_ids)} molecule IDs. Starting downloads...")

    # Save metadata to CSV
    metadata_file = os.path.join(args.folder_name, "molecule_metadata.csv")
    pd.DataFrame(metadata).to_csv(metadata_file, index=False)
    print(f"Metadata saved to {metadata_file}")

    # Download structures in parallel
    with ThreadPoolExecutor(max_workers=args.processors) as executor:
        executor.map(
            lambda cid: download_structure(cid, folders, args.output_format, args.preprocess),
            chembl_ids
        )

    print(f"Downloaded {len(chembl_ids)} molecules.")

if __name__ == "__main__":
    args = parser.parse_args()
    fetch_and_download_molecules(args)
