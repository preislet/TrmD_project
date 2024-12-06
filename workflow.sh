#!/bin/sh
# MeetEU PROJECT WORKFLOW

# Exit on any command failure
set -e

# Configuration
TARGETS="4yvg"                              # Comma-separated list of target PDB-IDs
LIGANDS_FOLDER="./ligands"                  # Directory for ligand files
TARGETS_FOLDER="./targets"                  # Directory for target files
DOCKED_STRUCTURES="$OUTPUT_DIR/ligand_target_pairs.txt"
CONFIG_FILE="./params/unidock_config.json"  # Configuration file for docking
OUTPUT_DIR="./output"                       # Directory for output files

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Log start of workflow
echo "Starting MeetEU project workflow..."


# Step 1: Fetch ligands
# IF NOT ALREADY FETCHED
echo "Step 1: Fetching ligands from ChEMBL database..."
python3 download_molecules.py [OPTIONS]
echo "Ligands fetched successfully inside '$LIGANDS_FOLDER'."


# Step 2: Prepare ligands
echo "Step 2: Refining ligands in '$LIGANDS_FOLDER'..."
bash prepare_ligands.sh "$LIGANDS_FOLDER"


# Step 3: Fetch targets
# IF NOT ALREADY FETCHED
echo "Step 3: Fetching '$TARGETS'..."
bash get_targets.sh "$TARGETS" > "$PATH_TO_TARGETS"


# Step 4: Prepare targets
echo "Step 4: Refining targets in '$PATH_TO_TARGETS'..."
python3 prepare_targets.py [OPTIONS]


# Step 5: Perform virtual screening with docking
echo "Step 5: Screening viable candidates with docking..."
python3 unidock_screening.py [OPTIONS] > "$DOCKED_STRUCTURES"


# Step 6: Generate additional ligands with a generative model
echo "Step 7: Generating additional ligands with a generative model..."
python3 generate_ligands.py "$LIGAND_TARGET_PAIRS" > "$OUTPUT_DIR/generated_ligands.sdf"


# Step 7: Optimize ligand-target interactions
echo "Step 6: Optimizing ligand-target interactions..."
bash optimize_interactions.sh "$LIGAND_TARGET_PAIRS"


# Step 8: Calculate final scores for all possible ligands
echo "Step 8: Calculating final scores of ligand-target pairs..."
python3 generate_scores.py "$LIGAND_TARGET_PAIRS" > "$OUTPUT_DIR/scores.txt"


# Step 9: Generate final report
echo "Step 9: Generating final report..."
bash generate_report.sh "$LIGAND_TARGET_PAIRS" "$OUTPUT_DIR"


# Log completion
echo "Workflow completed successfully. Results are saved in '$OUTPUT_DIR'."
