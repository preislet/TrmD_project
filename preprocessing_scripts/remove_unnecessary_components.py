import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

# Argument parser for input and output folder paths
parser = argparse.ArgumentParser(description="Process SDF files to retain the largest connected structure.")
parser.add_argument("input_folder", type=str, help="Path to the folder containing input SDF files.")
parser.add_argument("output_folder", type=str, help="Path to the folder to save processed SDF files.")

def remove_components(input_folder, output_folder):
    """
    Processes all SDF files in the input folder:
      1. Retains the largest connected molecular structure (if multiple components exist).
      2. Optimizes the structure to generate new conformations independent of deleted structures.
      3. Saves the processed molecule to the output folder.

    Args:
        input_folder (str): Path to the folder containing input SDF files.
        output_folder (str): Path to the folder to save processed SDF files.
    """
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    for sdf_file in os.listdir(input_folder):
        if sdf_file.endswith(".sdf"):
            input_path = os.path.join(input_folder, sdf_file)
            output_path = os.path.join(output_folder, sdf_file)

            # Read the molecule from the SDF file
            supplier = Chem.SDMolSupplier(input_path, removeHs=False)
            molecule = next(supplier)

            if molecule is None:
                print(f"Error: Could not read molecule from {sdf_file}")
                continue

            # Find connected components
            components = Chem.rdmolops.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)

            if len(components) > 1:
                # Select the largest connected component
                largest_component = max(components, key=lambda mol: mol.GetNumAtoms())
                molecule = largest_component
                # Optimize the structure (if embedded)
                AllChem.UFFOptimizeMolecule(molecule)

            # Write the processed molecule to the output SDF
            writer = Chem.SDWriter(output_path)
            writer.write(molecule)
            writer.close()

            print(f"Processed {sdf_file}: Saved to {output_path}")

if __name__ == "__main__":    
    args = parser.parse_args()
    remove_components(args.input_folder, args.output_folder)
