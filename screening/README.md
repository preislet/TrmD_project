# Uni-Dock: High-Throughput Docking Workflow

This project provides a Python script for batch-processing ligand docking using the Uni-Dock tool. It supports flexible scoring functions, advanced batch management, and automated output handling, making it suitable for high-throughput virtual screening workflows.

---

## Features

- **Batch Processing**: Efficiently split large ligand files into smaller batches for parallel docking.
- **Flexible Scoring**: Choose from `ad4`, `vina`, or `vinardo` scoring functions.
- **Comprehensive Input Support**: Accepts rigid receptors, flexible side chains, single ligands, or ligand lists in `PDBQT` or `SDF` format.
- **Grid-Based Docking**: Define custom search space dimensions or use automatic ligand-based box generation.
- **Results Extraction**: Extract docking scores and save results to a CSV file.
- **Advanced Configuration**: Fine-tune docking parameters such as exhaustiveness, binding modes, and search modes.
- **GPU Support**: Enable GPU-accelerated docking for batch ligands.

---

## Requirements

### Uni-Dock Software

Ensure that the Uni-Dock software is installed and available in your system's PATH.

---

## Usage

Run the script from the command line to perform docking operations. The script supports various options to customize receptor setup, ligand handling, and docking parameters.

### Basic Command
```bash
python3 script.py [OPTIONS]
```
### Command-Line Arguments

#### Input Arguments
| Argument        | Type   | Default | Description                                                |
|-----------------|--------|---------|------------------------------------------------------------|
| `--receptor`    | `str`  | None    | Rigid part of the receptor (`PDBQT` or `PDB`).             |
| `--flex`        | `str`  | None    | Flexible side chains (`PDBQT` or `PDB`).                  |
| `--ligand`      | `str`  | None    | Single ligand file (`PDBQT`).                             |
| `--ligand_index`| `str`  | None    | File containing paths to multiple ligands (`PDBQT`, `SDF`).|
| `--batch`       | `str`  | None    | Batch ligand file (`PDBQT`).                              |
| `--gpu_batch`   | `str`  | None    | GPU-accelerated batch ligand file (`PDBQT`, `SDF`).       |

#### Docking Parameters
| Argument            | Type   | Default  | Description                                           |
|---------------------|--------|----------|-------------------------------------------------------|
| `--scoring`         | `str`  | `vina`   | Scoring function (`ad4`, `vina`, `vinardo`).          |
| `--exhaustiveness`  | `int`  | `8`      | Search exhaustiveness.                                |
| `--num_modes`       | `int`  | `9`      | Maximum number of binding modes to generate.          |
| `--min_rmsd`        | `float`| `1.0`    | Minimum RMSD between output poses.                   |
| `--energy_range`    | `float`| `3.0`    | Maximum energy difference (kcal/mol).                |

#### Search Space Parameters
| Argument        | Type   | Default | Description                                              |
|-----------------|--------|---------|----------------------------------------------------------|
| `--center_x`    | `float`| None    | X coordinate of the center (Angstrom).                  |
| `--center_y`    | `float`| None    | Y coordinate of the center (Angstrom).                  |
| `--center_z`    | `float`| None    | Z coordinate of the center (Angstrom).                  |
| `--size_x`      | `float`| `22.0`  | Size in X dimension (Angstrom).                         |
| `--size_y`      | `float`| `22.0`  | Size in Y dimension (Angstrom).                         |
| `--size_z`      | `float`| `22.0`  | Size in Z dimension (Angstrom).                         |
| `--autobox`     | `bool` | `False` | Automatically set grid dimensions based on input ligand.|

#### Output Arguments
| Argument       | Type   | Default            | Description                                              |
|----------------|--------|--------------------|----------------------------------------------------------|
| `--out`        | `str`  | None               | Output file for docking results (`PDBQT`).               |
| `--dir`        | `str`  | `./result`         | Output directory for batch docking.                      |
| `--csv_output` | `str`  | `docking_results.csv` | Name of the CSV file for extracted docking scores.        |

#### Miscellaneous
| Argument          | Type   | Default | Description                                             |
|-------------------|--------|---------|---------------------------------------------------------|
| `--cpu`           | `int`  | `0`     | Number of CPUs to use.                                  |
| `--max_gpu_memory`| `int`  | `0`     | Maximum GPU memory to use (MB).                         |
| `--verbosity`     | `int`  | `1`     | Verbosity level (`0=no output`, `1=normal`, `2=verbose`).|

---

### Example Usage

#### Single Ligand Docking
```bash
python3 script.py --receptor receptor.pdbqt --ligand ligand.pdbqt --center_x 10.0 --center_y 10.0 --center_z 10.0
```
#### Batch Docking
```bash
python3 script.py --receptor receptor.pdbqt --ligand_index ligands.txt --center_x 10.0 --center_y 10.0 --center_z 10.0 --batch_size 100 --csv_output results.csv
```
#### GPU-Accelerated Batch Docking
```bash
python3 script.py --receptor receptor.pdbqt --gpu_batch ligands.sdf --center_x 15.0 --center_y 15.0 --center_z 15.0
```
---

## Folder Structure

The results are organized in the specified output directory:

result/
├── ligand1.pdbqt
├── ligand2.pdbqt
└── docking_results.csv

---

## Advanced Features

- **Batch Management**: Efficiently split large ligand files for processing in smaller chunks.
- **Automated Scoring Extraction**: Generate CSV files summarizing docking scores for downstream analysis.
- **Flexible Output Options**: Save docking results and grid maps for visualization.

