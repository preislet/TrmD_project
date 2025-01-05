import sys
import os
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

parser = argparse.ArgumentParser(description="Preprocess ligands for docking.")
parser.add_argument("input_dir", help="Directory containing input SDF files.")
parser.add_argument("output_dir", help="Directory to store processed SDF files.")
parser.add_argument("--generate_conformers", action="store_true", default=False, help="Generate conformers for each ligand.")
parser.add_argument("--num_conformers", type=int, default=10, help="Number of conformers to generate.")
parser.add_argument("--optimize", action="store_true", default=False, help="Perform a geometry optimization on each ligand.")
parser.add_argument("--single_file_for_all_conformers", action="store_true", default=False, help="Write all conformers to a single file.")
parser.add_argument("--mol_weight_limit", type=float, default=1000, help="Maximum molecular weight of ligands to process.")


def _generate_conformers(mol, num_conformers, max_attempts=1000):
    # Embed multiple conformers
    params = AllChem.ETKDGv3()  # Use appropriate parameters for your RDKit version
    params.maxAttempts = max_attempts
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)
    if not cids:
        raise ValueError("Failed to generate conformers.")
    # Optimize each conformer
    for cid in cids:
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=200)
        except Exception as e:
            print(f"Conformer {cid} optimization failed: {e}", file=sys.stderr)
    return cids


def preprocess_ligand(input_sdf, output_sdf, optimize=False, generate_conformers=False, num_conformers=0, multiconf_file=False, mol_weight_limit=1000):
    try:
        # Load the ligand
        supplier = Chem.SDMolSupplier(input_sdf, removeHs=True)
        mols = [m for m in supplier if m is not None]
        if len(mols) == 0:
            raise ValueError(f"No valid molecules found in {input_sdf}.")
        mol = mols[0]

        # Check molecular weight
        mol_weight = Descriptors.MolWt(mol)
        if mol_weight > mol_weight_limit:
            print(f"Skipping ligand with molecular weight {mol_weight} > {mol_weight_limit}.")
            return
        print(f"Molecular weight: {mol_weight}")

        # Clean-up before proceeding
        try:
            Chem.SanitizeMol(mol)
        except Chem.MolSanitizeException as e:
            print(f"Sanitization failed for {input_sdf}: {e}")
            return
        print("SANITIZED")

        # Add hydrogens
        mol = Chem.AddHs(mol, addCoords=True)
        print("HYDROGENS ADDED")

        # Compute Gasteiger charges (optional)
        AllChem.ComputeGasteigerCharges(mol)
        print("CHARGES COMPUTED")

        # Optional geometry optimization
        if optimize:
            if mol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                # Or: AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except Exception as e:
                print(f"Geometry optimization failed for {input_sdf}: {e}", file=sys.stderr)

        # Optional conformer generation
        if generate_conformers:
            if num_conformers <= 0:
                raise ValueError("Number of conformers must be a positive integer.")
            cids = _generate_conformers(mol, num_conformers=num_conformers)
            print(f"Generated {len(cids)} conformers for {input_sdf}")

        # Write to output SDF
        if multiconf_file == True:
            print(f"Writing all conformers to a single file: {output_sdf}")
            with Chem.SDWriter(output_sdf) as writer:
                print(mol.GetNumConformers())
                for conf_id in range(mol.GetNumConformers()):
                    writer.write(mol, confId=conf_id)
        else:
            # Write all conformers in separate files
            print(f"Writing conformers to separate files: {output_sdf}_confX.sdf")
            for conf_id in range(mol.GetNumConformers()):
                output_path = output_sdf.replace(".sdf", f"_conf{conf_id}.sdf")
                with Chem.SDWriter(output_path) as writer:
                    writer.write(mol, confId=conf_id)

    except Exception as e:
        raise ValueError(f"Error processing {input_sdf}: {e}")


def main(args):
    input_dir = args.input_dir
    output_dir = args.output_dir
    optimize = args.optimize
    generate_conformers = args.generate_conformers
    num_conformers = args.num_conformers
    multiconf_file = args.single_file_for_all_conformers
    mol_weight_limit = args.mol_weight_limit

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Process all SDF files in the input directory
    for filename in os.listdir(input_dir):
        if filename.lower().endswith(".sdf"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)
            try:
                print(f"Processing: {input_path}")
                preprocess_ligand(input_path, output_path, optimize, generate_conformers, num_conformers, multiconf_file, mol_weight_limit)
                print(f"Processed ligand successfully!")
            except Exception as e:
                print(f"Error processing {input_path}: {e}", file=sys.stderr)
                continue


if __name__ == "__main__":
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
