#!/usr/bin/env python3

import sys
from Bio.PDB import PDBParser, PDBIO, Select

class PocketSelect(Select):
    """
    A custom Biopython 'Select' subclass to filter atoms
    based on whether they fall within a specified 3D box.
    """
    def __init__(self, x_min, x_max, y_min, y_max, z_min, z_max):
        super().__init__()
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max

    def accept_atom(self, atom):
        x, y, z = atom.get_coord()
        if (self.x_min <= x <= self.x_max and
            self.y_min <= y <= self.y_max and
            self.z_min <= z <= self.z_max):
            return True
        return False

def main():
    if len(sys.argv) != 8:
        print(f"Usage: {sys.argv[0]} <PDB_FILE> <x_center> <y_center> <z_center> <size_x> <size_y> <size_z>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    x_center = float(sys.argv[2])
    y_center = float(sys.argv[3])
    z_center = float(sys.argv[4])
    size_x   = float(sys.argv[5])
    size_y   = float(sys.argv[6])
    size_z   = float(sys.argv[7])

    # Calculate bounding box
    x_min = x_center - size_x / 2.0
    x_max = x_center + size_x / 2.0
    y_min = y_center - size_y / 2.0
    y_max = y_center + size_y / 2.0
    z_min = z_center - size_z / 2.0
    z_max = z_center + size_z / 2.0

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input_structure", pdb_file)

    # Create PDBIO object and save filtered structure
    io = PDBIO()
    io.set_structure(structure)
    io.save("pocket_cut.pdb", select=PocketSelect(x_min, x_max, y_min, y_max, z_min, z_max))

    print("Cut-out pocket file saved to: pocket_cut.pdb")

if __name__ == "__main__":
    main()
