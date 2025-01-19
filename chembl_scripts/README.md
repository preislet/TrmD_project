# ChEMBL Molecule Downloader

This project provides a Python script for downloading molecules from the ChEMBL database with features such as filtering by approval phase (`max_phase`), preprocessing using RDKit, and organizing outputs into separate folders for original and preprocessed files.

---

## Features

- **Download Molecules**: Fetch and save molecules in formats such as `sdf`, `mol`, or `smiles`.
- **Approval Phase Filtering**: Filter molecules based on `max_phase` (e.g., FDA-approved drugs).
- **Preprocessing**: Optionally preprocess molecules using RDKit to generate conformers.
- **Parallel Downloads**: Speed up downloading by leveraging multithreading.
- **Organized Output**: Automatically create separate folders for original and preprocessed molecules.

---

## Requirements

### Python Libraries

The script requires the following Python libraries:

- `requests`
- `rdkit`

You can install the dependencies using:

```bash
pip install -r requirements.txt
```

## Usage

Run the script from the command line with various options to control molecule selection, output format, parallelization, and preprocessing.

```bash
python3 download_molecules.py [OPTIONS]
```

### Command-Line Arguments

| Argument         | Type    | Default               | Description                                                 |
|------------------|---------|-----------------------|-------------------------------------------------------------|
| `--max_results`  | `int`   | `1000`                 | Maximum number of molecules to download.                    |
| `--output_format`| `str`   | `sdf`                | Format of the molecule files (`sdf`, `mol`, or `smiles`).    |
| `--processors`   | `int`   | `4`                  | Number of threads to use for parallel downloads.             |
| `--folder_name`  | `str`   | `molecule_structures`| Base folder to store downloaded and preprocessed files.     |
| `--max_phase`    | `int`   | `4`                  | Maximum approval phase (e.g., 4 for FDA-approved drugs).    |
| `--preprocess`   | `bool`  | `False`              | Enable preprocessing of molecules using RDKit.              |
| `--verbose`      | `bool`  | `False`              | Activate verbose mode                                       |


## Folder Structure

The script automatically organizes downloaded files into original and preprocessed folders within the specified folder_name.

``` plaintext
molecule_structures/
│
├── original/
│   ├── CHEMBL1.sdf
│   ├── CHEMBL2.sdf
│   └── ...
│
└── preprocessed/
    ├── CHEMBL1_preprocessed.sdf
    ├── CHEMBL2_preprocessed.sdf
    └── ...
```

# Advanced ChEMBL Molecule Downloader

This Python script provides a powerful tool to download and preprocess molecules from the ChEMBL database, with support for advanced filtering based on molecular properties.

---

## Features

- **Download Molecules**: Retrieve molecule data from ChEMBL in `sdf`, `mol`, or `smiles` formats.
- **Advanced Filtering**: Apply filters for molecular weight, logP, PSA, Lipinski Rule violations, and more.
- **Metadata Export**: Save detailed metadata for downloaded molecules in a CSV file.
- **Preprocessing**: Use RDKit to preprocess downloaded molecules, including conformer generation.
- **Parallel Downloads**: Speed up downloading with multithreading.

---

## Requirements

### Python Libraries

Install the required libraries:

```bash
pip install -r requirements.txt
```

## Usage

Run the script from the command line with various options to control molecule selection, output format, parallelization, and preprocessing.

```bash
python3 download_molecules_advanced_filter.py [OPTIONS]
```

### Command-Line Arguments

| Argument                   | Type    | Default               | Description                                                  |
|----------------------------|---------|-----------------------|--------------------------------------------------------------|
| `--max_results`            | `int`   | `1000`                  | Maximum number of molecules to download.                     |
| `--output_format`          | `str`   | `sdf`                 | Format of the molecule files (`sdf`, `mol`, `smiles`).        |
| `--processors`             | `int`   | `4`                   | Number of threads for parallel downloads.                    |
| `--folder_name`            | `str`   | `molecule_structures` | Folder to store original and preprocessed molecules.         |
| `--max_phase`              | `int`   | `4`                   | Maximum approval phase (e.g., 4 for FDA-approved drugs).      |
| `--min_molecular_weight`   | `float` | `None`                | Minimum molecular weight for filtering.                      |
| `--max_molecular_weight`   | `float` | `None`                | Maximum molecular weight for filtering.                      |
| `--min_logp`               | `float` | `None`                | Minimum logP for filtering.                                  |
| `--max_logp`               | `float` | `None`                | Maximum logP for filtering.                                  |
| `--min_hba`                | `int`   | `None`                | Minimum hydrogen bond acceptors.                             |
| `--max_hbd`                | `int`   | `None`                | Maximum hydrogen bond donors.                                |
| `--max_rotatable_bonds`    | `int`   | `None`                | Maximum number of rotatable bonds.                           |
| `--min_psa`                | `float` | `None`                | Minimum polar surface area (PSA) for filtering.              |
| `--max_psa`                | `float` | `None`                | Maximum polar surface area (PSA) for filtering.              |
| `--max_lipinski_violations`| `int`   | `None`                | Maximum Lipinski Rule of 5 violations allowed.               |
| `--min_np_likeness`        | `float` | `None`                | Minimum natural product likeness score.                      |
| `--max_np_likeness`        | `float` | `None`                | Maximum natural product likeness score.                      |
| `--molecular_species`      | `str`   | `None`                | Filter by molecular species (`NEUTRAL`, `ACID`, `BASE`).      |
| `--preprocess`             | `bool`  | `False`               | Enable preprocessing of molecules using RDKit.               |
| `--verbose`                | `bool`  | `False`               | Activate verbose mode                                       |

### Example

```bash
python3 download_molecules_advanced_filter.py --max_results 5000 --min_molecular_weight 250 --max_molecular_weight 400 --min_logp 0 --max_logp 5 --min_psa 70 --max_psa 120 --max_rotatable_bonds 2
```
