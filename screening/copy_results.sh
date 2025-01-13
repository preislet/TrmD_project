#!/bin/bash

# Define the remote path, local destination, and number of top ligands to download
REMOTE_PATH="metacentrum:/storage/praha1/home/tobiasma/meet-eu/trmd-project"
LOCAL_PATH="./adomet_pocket"
TOP_N=6  # Number of top ligands to process

# Temporary file to store top ligand names
TOP_LIG_NAMES=$(mktemp)

# Function to handle file transfer
transfer_files() {
    local molecule="$1"

    # Extract the part before the first "_"
    if [[ "$molecule" == *"_"* ]]; then
        molecule="${molecule%%_*}"
    fi

    local molecule_dir="${LOCAL_PATH}/${molecule}"

    # Check if the directory exists, and skip the molecule if it does
    if [[ -d "${molecule_dir}" ]]; then
        echo "Directory ${molecule_dir} already exists. Skipping molecule ${molecule}."
        return
    fi

    # Create necessary directories
    mkdir -p "${molecule_dir}/out"

    # Define remote paths for different files
    scp "${REMOTE_PATH}/result/${molecule}[!0-9]*.sdf" "${molecule_dir}/out/" 2>/dev/null
    scp "${REMOTE_PATH}/molecule_structures/original/${molecule}.sdf" "${molecule_dir}/" 2>/dev/null
    scp "${REMOTE_PATH}/molecule_structures/processed_single_components/${molecule}[!0-9]*.sdf" "${molecule_dir}/" 2>/dev/null
    scp "${REMOTE_PATH}/molecule_structures/preprocessed/${molecule}[!0-9]*.sdf" "${molecule_dir}/" 2>/dev/null
}

# Extract top ligand names based on docking results
if [ -f "${LOCAL_PATH}/docking_results.csv" ]; then
    cat "${LOCAL_PATH}/docking_results.csv" | sort -t, -k2,2 -n | head -n "${TOP_N}" | cut -d, -f1 > "${TOP_LIG_NAMES}"
else
    echo "Error: docking_results.csv not found in ${LOCAL_PATH}"
    exit 1
fi

# Loop through the molecule names and copy each file
while read -r molecule; do
    if [ -n "${molecule}" ]; then
        transfer_files "${molecule}"
    else
        echo "Warning: Skipping empty molecule name."
    fi
done < "${TOP_LIG_NAMES}"

# Remove temporary file
rm -f "${TOP_LIG_NAMES}"

echo "File transfer completed successfully."

