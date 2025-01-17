# Uni-Dock: High-Throughput Docking Workflow

This project provides a Python script for batch-processing ligand docking using the Uni-Dock tool (GPU-accelerated molecular docking tool). It supports flexible scoring functions, advanced batch management, and automated output handling, making it suitable for high-throughput virtual screening workflows.

---

## Requirements

### Uni-Dock Software

Ensure that the Uni-Dock software is installed and available in your system's PATH or in the path of your container. You can follow the installation steps [here](https://github.com/dptech-corp/Uni-Dock/tree/main/unidock#installation) or you can use the Docker image provided in our repository. One thing to note is that you have to have sufficient NVIDIA drivers.

---

## Usage

Run the script from the command line to perform docking operations and save the results in a CSV table. The script supports various options to customize receptor setup, ligand handling, and docking parameters.

### Basic Command

The basic command should include a receptor file, ligand file or a txt file of paths to ligand files (on a single line, space-separated like this: `./ligands/lig1.pdbqt ./ligands/lig2.pdbqt ...`) and the information about the location of where to dock.

```bash
python3 run_unidock.py --receptor <target_file> --ligand|ligand_index <ligand_file|ligands.txt> \
     -cx <center_x> -cy <center_y> -cz <center_z> [OPTIONS]
```

```bash
python3 run_unidock.py --receptor <receptor.pdbqt> \
     --ligand_index ligands.txt \
     --search_mode balance \
     --scoring vina \
     -cx <center_x> -cy <center_y> -cz <center_z> \
     -sx <size_x> -sy <size_y> -sz <size_z> \
     --num_modes 1 \
     --dir <save dir>
```

### Command-Line Arguments

#### Input Arguments
| Argument                      | Type   | Default | Description                                                |
|-------------------------------|--------|---------|------------------------------------------------------------|
| `--receptor`, `-r`            | `str`  | None    | Rigid part of the receptor (`PDBQT` or `PDB`).             |
| `--flex`, `-f`                | `str`  | None    | Flexible side chains (`PDBQT` or `PDB`).                  |
| `--ligand`, `-l`              | `str`  | None    | Single ligand file (`PDBQT`).                             |
| `--ligand_index`, `-li`       | `str`  | None    | File containing paths to multiple ligands (`PDBQT`, `SDF`).|
| `--batch`, `-b`               | `str`  | None    | Batch ligand file (`PDBQT`).                              |
| `--gpu_batch`, `-gb`          | `str`  | None    | GPU-accelerated batch ligand file (`PDBQT`, `SDF`).       |

#### Docking Parameters
| Argument                      | Type    | Default  | Description                                           |
|-------------------------------|---------|----------|-------------------------------------------------------|
| `--scoring`, `-s`             | `str`   | `vina`   | Scoring function (`ad4`, `vina`, `vinardo`).          |
| `--exhaustiveness`, `-ex`     | `int`   | `8`      | Search exhaustiveness.                                |
| `--num_modes`, `-nm`          | `int`   | `9`      | Maximum number of binding modes to generate.          |
| `--min_rmsd`, `-mr`           | `float` | `1.0`    | Minimum RMSD between output poses.                   |
| `--energy_range`, `-er`       | `float` | `3.0`    | Maximum energy difference (kcal/mol).                |
| `--max_evals`, `-me`          | `int`   | `0`      | Number of evaluations in each Monte Carlo (MC) run.  |
| `--max_step`, `-ms`           | `int`   | `0`      | Maximum number of steps in each MC run.              |
| `--refine_step`, `-rs`        | `int`   | `5`      | Number of steps in refinement.                       |
| `--search_mode`, `-sm`        | `str`   | None     | Search mode (`fast`, `balance`, `detail`).           |

#### Search Space Parameters
| Argument                      | Type    | Default | Description                                              |
|-------------------------------|---------|---------|----------------------------------------------------------|
| `--center_x`, `-cx`           | `float` | None    | X coordinate of the center (Angstrom).                  |
| `--center_y`, `-cy`           | `float` | None    | Y coordinate of the center (Angstrom).                  |
| `--center_z`, `-cz`           | `float` | None    | Z coordinate of the center (Angstrom).                  |
| `--size_x`, `-sx`             | `float` | `22.0`  | Size in X dimension (Angstrom).                         |
| `--size_y`, `-sy`             | `float` | `22.0`  | Size in Y dimension (Angstrom).                         |
| `--size_z`, `-sz`             | `float` | `22.0`  | Size in Z dimension (Angstrom).                         |
| `--autobox`                   | `bool`  | `False` | Automatically set grid dimensions based on input ligand.|

#### Output Arguments
| Argument                      | Type   | Default            | Description                                              |
|-------------------------------|--------|--------------------|----------------------------------------------------------|
| `--out`, `-o`                 | `str`  | None               | Output file for docking results (`PDBQT`).               |
| `--dir`, `-d`                 | `str`  | `./result`         | Output directory for batch docking.                      |
| `--csv_output`, `-csv`        | `str`  | `docking_results.csv` | Name of the CSV file for extracted docking scores.      |
| `--write_maps`, `-wm`         | `str`  | None               | Output filename (directory + prefix name) for maps.      |

#### Miscellaneous
| Argument                      | Type    | Default | Description                                             |
|-------------------------------|---------|---------|---------------------------------------------------------|
| `--cpu`                       | `int`   | `0`     | Number of CPUs to use.                                  |
| `--max_gpu_memory`, `-mgm`    | `int`   | `0`     | Maximum GPU memory to use (MB).                         |
| `--spacing`, `-sp`            | `float` | `0.375` | Grid spacing (Angstrom).                                |
| `--seed`                      | `int`   | `42`    | Explicit random seed.                                   |
| `--verbosity`, `-v`           | `int`   | `1`     | Verbosity level (`0=no output`, `1=normal`, `2=verbose`).|
| `--batch_size`, `-bs`         | `int`   | `0`     | Number of ligands to process in each batch.             |
| `--config`, `-c`              | `str`   | None    | Path to a configuration file.                          |
| `--version`                   | `bool`  | `False` | Display program version.                                |

### Example Usage

#### Single Ligand Docking
```bash
python3 run_unidock.py -r receptor.pdbqt --ligand ligand.pdbqt -cx 10.0 -cy 10.0 -cz 10.0
```
#### Batch Docking
```bash
python3 run_unidock.py -r receptor.pdbqt --ligand_index ligands.txt -cx 10.0 -cy 10.0 -cz 10.0 --batch_size 100 --csv_output results.csv
```
---

## Folder Structure

The results are organized in the specified output directory:

```plaintext
result/
├── docking_results.csv
└── <target_name>/
    ├── lig1_out.pdbqt
    ├── lig2_out.pdbqt
    └── ...
```

---

## Features

- **Batch Processing**: Efficiently split large ligand files into smaller batches for parallel docking.
- **Flexible Scoring**: Choose from `ad4`, `vina`, or `vinardo` scoring functions.
- **Comprehensive Input Support**: Accepts rigid receptors, flexible side chains, single ligands, or ligand lists in `PDBQT` or `SDF` format.
- **Grid-Based Docking**: Define custom search space dimensions or use automatic ligand-based box generation.
- **Results Extraction**: Extract docking scores and save results to a CSV file.
- **Advanced Configuration**: Fine-tune docking parameters such as exhaustiveness, binding modes, and search modes.
- **GPU Support**: Enable GPU-accelerated docking for batch ligands.

