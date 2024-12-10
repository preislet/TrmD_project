## Scripts

### `generate_conformers.py`

**Purpose**:  
Generates 3D conformers for a given SMILES string.

**Inputs**:
- SMILES representation of a molecule.
- Number of conformers to generate.

**Outputs**:  
A directory containing the generated conformers in PDB format.


### `prepare_targets.py`

**Purpose**:  
Prepares a target (protein) for docking by:
- Removing water molecules.
- Adding missing hydrogen atoms.

**Inputs**:
- PDB file of the target protein.

**Outputs**:  
A cleaned-up PDB file ready for docking.

**Usage**:
```bash
python generate_conformers.py
python prepare_targets.py


