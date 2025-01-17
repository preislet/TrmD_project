from pymol import cmd
import os

def prepare_target(pdb_file, output_dir, clean_water=True, add_hydrogens=True):
    """
    Připraví protein na docking: odstraní vodu, přidá vodíky a uloží strukturu.

    Parameters:
        pdb_file (str): Cesta k PDB souboru targetu.
        output_dir (str): Složka, kam se uloží připravený target.
        clean_water (bool): Odstranit molekuly vody (HOH).
        add_hydrogens (bool): Přidat vodíky k atomům.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Načti protein do PyMOL
    object_name = os.path.basename(pdb_file).replace(".pdb", "")
    cmd.load(pdb_file, object_name)

    # Volitelně odstranění molekul vody
    if clean_water:
        cmd.remove("resn HOH")
        print(f"Molekuly vody byly odstraněny pro {object_name}.")

    # Přidání vodíků
    if add_hydrogens:
        cmd.h_add("all")
        print(f"Atomové vodíky byly přidány pro {object_name}.")

    # Uložení vyčištěného PDB souboru
    output_file = os.path.join(output_dir, os.path.basename(pdb_file).replace(".pdb", "_prepared.pdb"))
    cmd.save(output_file, object_name)
    print(f"Připravený protein uložen jako: {output_file}")

    # Ukonči PyMOL relaci
    cmd.reinitialize()

# Použití
if __name__ == "__main__":
    # Cesty k PDB souborům (příklad více struktur)
    pdb_files = ["/Structure/TrmD/4yvg_manipulation/raw_4yvg.pdb"]

    # Výstupní složka
    output_dir = "prepared_targets"

    # Spusť přípravu pro každý target
    for pdb_file in pdb_files:
        prepare_target(pdb_file, output_dir, clean_water=True, add_hydrogens=True)

