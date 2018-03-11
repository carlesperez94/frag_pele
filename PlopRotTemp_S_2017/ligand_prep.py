import os
import sys
from schrodinger import structure as st
import main as plop


def convert_mae(ligands):
    """
       Desciption: From each structure retrieve
       a .mae file of the ligand in the receptor.
       Output:
            structure_mae: ligand
            res = residue
    """

    
    for structure in st.StructureReader(ligands):
        for residue in structure.residue:
            res = residue.pdbres.strip()
        str_name = "{}".format(res)
        try:
            structure.write(str_name + ".mae")
        except ValueError:
            str_name = "{}".format(res)
        finally:
            structure_mae = "{}.mae".format(str_name)
            structure.write(structure_mae)
    return structure_mae

def create_template(pdb):
   mae_file = convert_mae(pdb)
   plop.main(mae_file)
   os.remove(mae_file)

