#!/bin/bash

# Define the remote path and local destination
REMOTE_PATH="metacentrum:/storage/praha1/home/tobiasma/meet-eu/trmd-project"
LOCAL_PATH="."

# Loop through the molecule names and copy each file
while read -r MOLECULE; do
    mkdir -p "${MOLECULE}/out"
    scp "${REMOTE_PATH}/result/${MOLECULE}[!0-9]*sdf" "${LOCAL_PATH}/${MOLECULE}/out/"
    scp "${REMOTE_PATH}/molecule_structures/original/${MOLECULE}[!0-9]*sdf" "${LOCAL_PATH}/${MOLECULE}/"
    scp "${REMOTE_PATH}/molecule_structures/preprocessed_w_h_separate/${MOLECULE}[!0-9]*sdf" "${LOCAL_PATH}/${MOLECULE}/"
    scp "${REMOTE_PATH}/molecule_structures/preprocessed/${MOLECULE}[!0-9]*sdf" "${LOCAL_PATH}/${MOLECULE}/"
done < top_lig_names.txt

