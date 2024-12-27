import sys
import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser(description="Preprocess ligands for docking.")
parser.add_argument("input_dir", help="Directory containing input SDF files.")
parser.add_argument("output_dir", help="Directory to store processed SDF files.")
parser.add_argument("--optimize", action="store_true", default=False, help="Perform a geometry optimization on each ligand.")

def preprocess_ligand(input_sdf, output_sdf, optimize=True):
    # Load the ligand
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    mols = [m for m in supplier if m is not None]
    if len(mols) == 0:
        raise ValueError(f"No valid molecules found in {input_sdf}.")
    mol = mols[0]

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Compute Gasteiger charges (optional)
    AllChem.ComputeGasteigerCharges(mol)

    # Optional geometry optimization
    if optimize:
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        # Or: AllChem.MMFFOptimizeMolecule(mol, maxIters=200)

    # Write to output SDF
    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()


def main(args):
    input_dir = args.input_dir
    output_dir = args.output_dir
    optimize = args.optimize

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Process all SDF files in the input directory
    for filename in os.listdir(input_dir):
        if filename.lower().endswith(".sdf"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)
            try:
                preprocess_ligand(input_path, output_path, optimize=optimize)
                print(f"Processed: {input_path} -> {output_path}")
            except Exception as e:
                print(f"Error processing {input_path}: {e}", file=sys.stderr)


if __name__ == "__main__":
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
