# ChEMBL Molecule Downloader

This project provides a Python script for downloading molecules from the ChEMBL database with features such as filtering by approval phase (`max_phase`), preprocessing using RDKit, and organizing outputs into separate folders for original and preprocessed files.

---

## Features

- **Download Molecules**: Fetch and save molecules in formats such as `sdf`, `mol`, or `smiles`.
- **Approval Phase Filtering**: Filter molecules based on `min_phase` (e.g., FDA-approved drugs).
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
| `--max_results`  | `int`   | `10`                 | Maximum number of molecules to download.                    |
| `--output_format`| `str`   | `sdf`                | Format of the molecule files (`sdf`, `mol`, or `smiles`).    |
| `--processors`   | `int`   | `4`                  | Number of threads to use for parallel downloads.             |
| `--folder_name`  | `str`   | `molecule_structures`| Base folder to store downloaded and preprocessed files.     |
| `--min_phase`    | `int`   | `4`                  | Minimum approval phase for filtering molecules.             |
| `--preprocess`   | `bool`  | `False`              | Enable preprocessing of molecules using RDKit.              |

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
