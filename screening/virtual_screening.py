#!/usr/bin/env python3

import argparse
import subprocess
import os
import glob
import time
import csv

# Argument parser
parser = argparse.ArgumentParser(description="Uni-Dock Docking Workflow with Batch Processing")

# Input arguments
parser.add_argument("--receptor", required=True, help="Rigid part of the receptor (PDBQT or PDB)")
parser.add_argument("--flex", help="Flexible side chains, if any (PDBQT or PDB)")
parser.add_argument("--ligand", help="Single ligand (PDBQT)")
parser.add_argument("--ligand_index", help="File containing paths to ligands (PDBQT or SDF)")
parser.add_argument("--batch", help="Batch ligand (PDBQT)")
parser.add_argument("--gpu_batch", help="GPU batch ligand (PDBQT or SDF)")
parser.add_argument("--scoring", choices=["ad4", "vina", "vinardo"], default="vina", help="Scoring function")

# Search space arguments
parser.add_argument("--maps", help="Affinity maps for AD4.2 or vina scoring function")
parser.add_argument("--center_x", type=float, required=True, help="X coordinate of the center (Angstrom)")
parser.add_argument("--center_y", type=float, required=True, help="Y coordinate of the center (Angstrom)")
parser.add_argument("--center_z", type=float, required=True, help="Z coordinate of the center (Angstrom)")
parser.add_argument("--size_x", type=float, default=22.0, help="Size in X dimension (Angstrom)")
parser.add_argument("--size_y", type=float, default=22.0, help="Size in Y dimension (Angstrom)")
parser.add_argument("--size_z", type=float, default=22.0, help="Size in Z dimension (Angstrom)")
parser.add_argument("--autobox", action="store_true", help="Set map dimensions based on input ligand(s)")

# Output arguments
parser.add_argument("--out", help="Output models (PDBQT)")
parser.add_argument("--dir", default="./result", help="Output directory for batch mode")
parser.add_argument("--write_maps", help="Output filename (directory + prefix name) for maps")
parser.add_argument("--csv_output", default="docking_results.csv", help="Name of the output CSV file") # NEW

# Miscellaneous arguments
parser.add_argument("--cpu", type=int, default=0, help="Number of CPUs to use")
parser.add_argument("--seed", type=int, default=0, help="Explicit random seed")
parser.add_argument("--exhaustiveness", type=int, default=8, help="Exhaustiveness of the global search")
parser.add_argument("--max_evals", type=int, default=0, help="Number of evaluations in each MC run")
parser.add_argument("--num_modes", type=int, default=9, help="Maximum number of binding modes to generate")
parser.add_argument("--min_rmsd", type=float, default=1.0, help="Minimum RMSD between output poses")
parser.add_argument("--energy_range", type=float, default=3.0, help="Maximum energy difference (kcal/mol)")
parser.add_argument("--spacing", type=float, default=0.375, help="Grid spacing (Angstrom)")
parser.add_argument("--verbosity", type=int, default=1, help="Verbosity (0=no output, 1=normal, 2=verbose)")
parser.add_argument("--max_step", type=int, default=0, help="Maximum number of steps in each MC run")
parser.add_argument("--refine_step", type=int, default=5, help="Number of steps in refinement")
parser.add_argument("--max_gpu_memory", type=int, default=0, help="Maximum GPU memory to use (MB)")
parser.add_argument("--search_mode", choices=["fast", "balance", "detail"], help="Search mode of Uni-Dock")

# Configuration file
parser.add_argument("--config", help="Path to a configuration file")

# Help and information
parser.add_argument("--version", action="store_true", help="Display program version")

# Batch processing
parser.add_argument("--batch_size", "-bs", type=int, default=0, help="Number of ligands to process in each batch") # NEW


def split_ligands(ligand_index, batch_size):
    """Split the ligand files into smaller batches."""
    with open(ligand_index, "r") as f:
        ligands = f.read().splitlines()
    for i in range(0, len(ligands), batch_size):
        yield ligands[i:i + batch_size]


def write_ligand_file(ligands, batch_number):
    """Write ligands for the current batch to a temporary file."""
    batch_file = f"batch_{batch_number}.txt"
    with open(batch_file, "w") as f:
        f.write("\n".join(ligands))
    return batch_file


def build_command(batch_file, args):
    """Build the Uni-Dock command for the current batch."""
    command = [
        "unidock",
        f"--receptor {args.receptor}",
        f"--ligand_index {batch_file}",
        f"--scoring {args.scoring}",
        f"--center_x {args.center_x}",
        f"--center_y {args.center_y}",
        f"--center_z {args.center_z}",
        f"--size_x {args.size_x}",
        f"--size_y {args.size_y}",
        f"--size_z {args.size_z}",
        f"--exhaustiveness {args.exhaustiveness}",
        f"--num_modes {args.num_modes}",
        f"--refine_step {args.refine_step}",
        f"--dir {args.dir}"
    ]
    if args.flex:
        command.append(f"--flex {args.flex}")
    if args.maps:
        command.append(f"--maps {args.maps}")
    if args.autobox:
        command.append("--autobox")
    if args.out:
        command.append(f"--out {args.out}")
    if args.write_maps:
        command.append(f"--write_maps {args.write_maps}")
    if args.cpu:
        command.append(f"--cpu {args.cpu}")
    if args.seed:
        command.append(f"--seed {args.seed}")
    if args.max_evals:
        command.append(f"--max_evals {args.max_evals}")
    if args.min_rmsd:
        command.append(f"--min_rmsd {args.min_rmsd}")
    if args.energy_range:
        command.append(f"--energy_range {args.energy_range}")
    if args.spacing:
        command.append(f"--spacing {args.spacing}")
    if args.verbosity:
        command.append(f"--verbosity {args.verbosity}")
    if args.max_step:
        command.append(f"--max_step {args.max_step}")
    if args.max_gpu_memory:
        command.append(f"--max_gpu_memory {args.max_gpu_memory}")
    if args.search_mode:
        command.append(f"--search_mode {args.search_mode}")
    return " ".join(command)


def run_docking(batch_file, batch_number, args):
    """Run docking for a specific batch."""
    os.makedirs(args.dir, exist_ok=True)
    command = build_command(batch_file, args)
    print(f"Running batch {batch_number}: {command}")
    start_time = time.time()
    subprocess.run(command, shell=True, check=True)
    elapsed_time = time.time() - start_time
    print(f"Batch {batch_number} completed in {elapsed_time:.2f} seconds.")


def extract_scores(result_dir, csv_output):
    """Extract scores from the result files and write them to a CSV."""
    result_files = glob.glob(f"{result_dir}/*.pdbqt")
    with open(csv_output, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Ligand", "Score"])
        for result_file in result_files:
            with open(result_file, "r") as f:
                for line in f:
                    if line.startswith("REMARK VINA RESULT:"):
                        score = float(line.split()[3])
                        ligand_name = os.path.basename(result_file)
                        csvwriter.writerow([ligand_name, score])
                        break


def main(args):
    if not args.ligand_index or not os.path.exists(args.ligand_index):
        raise FileNotFoundError("Ligand index file not found.")
    batch_size = args.batch_size or 1  # Default to processing one ligand per batch if not specified
    batch_generator = split_ligands(args.ligand_index, batch_size)
    for batch_number, ligands in enumerate(batch_generator, start=1):
        batch_file = write_ligand_file(ligands, batch_number)
        try:
            run_docking(batch_file, batch_number, args)
        finally:
            os.remove(batch_file)
    print("Extracting scores and generating CSV...")
    extract_scores(args.dir, args.csv_output)
    print(f"Results saved to {args.csv_output}")


if __name__ == "__main__":
    args = parser.parse_args([] if "__file__" not in globals() else None)
    main(args)
