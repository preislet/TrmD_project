import os
import requests
import argparse
from typing import List
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem
from rdkit.Chem import AllChem

# Argument parser setup
parser = argparse.ArgumentParser(description="Quickly download molecules from ChEMBL with filtering and preprocessing.")
parser.add_argument("--max_results", type=int, default=500, help="Maximum number of molecules to download.")
parser.add_argument("--output_format", type=str, default="sdf", choices=["sdf", "mol", "smiles"], help="Output format for molecule structures.")
parser.add_argument("--processors", type=int, default=4, help="Number of processors to use for downloading.")
parser.add_argument("--folder_name", type=str, default="molecule_structures", help="Base folder to save molecule files.")
parser.add_argument("--max_phase", type=int, default=4, help="Maximum approval phase for filtering molecules (e.g., 4 for FDA-approved drugs).")
parser.add_argument("--preprocess", action="store_true", default=False, help="Preprocess molecules using RDKit.")
parser.add_argument("--max_attempts", type=int, default=5, help="Maximum number of download attempts per molecule.")
parser.add_argument("--verbose", action="store_true", default=False, help="Print verbose output.")
args = parser.parse_args()

def create_folders(base_folder: str) -> dict:
    """Creates folders for storing original and preprocessed files.

    Args:
        base_folder (str): The base folder for storing molecule files.

    Returns:
        dict: A dictionary containing paths for original and preprocessed folders.
    """
    original_folder = os.path.join(base_folder, "original")
    preprocessed_folder = os.path.join(base_folder, "preprocessed")
    os.makedirs(original_folder, exist_ok=True)
    os.makedirs(preprocessed_folder, exist_ok=True)
    if args.verbose: print(f"Created folders: {original_folder}, {preprocessed_folder}")
    # make log file
    log_file = os.path.join(base_folder, "download_log.txt")
    return {"original": original_folder, "preprocessed": preprocessed_folder, "log": log_file}

def preprocess_molecule(file_path: str, output_folder: str) -> None:
    """Preprocess a molecule file using RDKit.

    Args:
        file_path (str): Path to the molecule file to preprocess.
        output_folder (str): Folder to save the preprocessed molecule file.
    """
    mol = Chem.MolFromMolFile(file_path)
    if mol:
        AllChem.EmbedMolecule(mol)
        processed_path = os.path.join(output_folder, os.path.basename(file_path).replace(".mol", "_preprocessed.mol").replace(".sdf", "_preprocessed.sdf"))
        Chem.MolToMolFile(mol, processed_path)
        print(f"Preprocessed and saved to: {processed_path}")

def download_structure(chembl_id: str, folders: dict, output_format: str, preprocess: bool, retries:int = args.max_attempts) -> None:
    """Downloads a molecule structure in the specified format.

    Args:
        chembl_id (str): The ChEMBL ID of the molecule to download.
        folders (dict): Dictionary containing paths for original and preprocessed folders.
        output_format (str): The format of the molecule structure file.
        preprocess (bool): Whether to preprocess the molecule after downloading.
    """
    original_folder = folders["original"]
    preprocessed_folder = folders["preprocessed"]
    log_file = folders["log"]
    structure_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.{output_format}"
    file_path = os.path.join(original_folder, f"{chembl_id}.{output_format}")
    
    if os.path.exists(file_path):
        if args.verbose: print(f"Skipping already downloaded: {chembl_id}")
        return

    for _ in range(retries):
        response = requests.get(structure_url)
        if response.status_code == 200:
            with open(file_path, "wb") as file:
                file.write(response.content)
            if args.verbose: print(f"Downloaded: {chembl_id}.{output_format}")
            if preprocess and output_format in ["mol", "sdf"]:
                preprocess_molecule(file_path, preprocessed_folder)
            return
        else:
            if args.verbose: 
                print(f"Failed to download {chembl_id}. Status Code: {response.status_code}")
                print(f"Retrying download for {chembl_id}...")

    if args.verbose: print(f"Failed to download {chembl_id} after {retries} retries.")
    with open(log_file, "a") as log:
                log.write(f"Failed to download {chembl_id} | https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}\n")

def fetch_and_download_molecules(args: argparse.Namespace) -> None:
    """Fetches molecule IDs and downloads them in parallel.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    base_url: str = "https://www.ebi.ac.uk/chembl/api/data/molecule"
    params: dict = {
        "format": "json",
        "limit": 500,
        "offset": 0,
        "max_phase": args.max_phase
    }

    folders: dict = create_folders(args.folder_name)
    chembl_ids: List[str] = []
    downloaded_count: int = 0

    # Fetch molecule IDs in batches
    while downloaded_count < args.max_results:
        response = requests.get(base_url, params=params)
        if response.status_code != 200:
            if args.verbose: print(f"Failed to fetch molecule data: {response.status_code}")
            # log file
            with open(folders["log"], "a") as log:
                log.write(f"Failed to fetch molecule data: {response.status_code}\n")
            break

        data = response.json()
        molecules = data.get("molecules", [])
        if not molecules:
            break

        for molecule in molecules:
            chembl_id: str = molecule.get("molecule_chembl_id", "")
            if not chembl_id:
                continue
            chembl_ids.append(chembl_id)
            downloaded_count += 1
            if downloaded_count >= args.max_results:
                break

        params["offset"] += len(molecules)
        if args.verbose: print(f"Fetched {len(chembl_ids)} molecule IDs so far...")

    if args.verbose: print(f"Fetched {len(chembl_ids)} molecule IDs. Starting downloads...")

    # Download structures in parallel
    with ThreadPoolExecutor(max_workers=args.processors) as executor:
        executor.map(
            lambda cid: download_structure(cid, folders, args.output_format, args.preprocess),
            chembl_ids
        )

    if args.verbose: print(f"Downloaded {len(chembl_ids)} molecules.")

if __name__ == "__main__":
    fetch_and_download_molecules(args)
