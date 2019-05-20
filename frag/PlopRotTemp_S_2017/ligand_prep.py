import os
from schrodinger import structure as st
import main as plop
import subprocess
import argparse


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

def create_template(pdb, gridres):
   mae_file = convert_mae(pdb)
   plop.main(mae_file, out_temp="DataLocal/Templates/OPLS2005/HeteroAtoms/templates_generated", out_rot="DataLocal/LigandRotamerLibs/",
             gridres=gridres)
   os.remove(mae_file)

def arg_parse():
  parser = argparse.ArgumentParser()
  parser.add_argument("pdb", type=str, help="ligand file to templatize")
  parser.add_argument("gridres", type=str, help="Degrees of rotation.")
  args = parser.parse_args()
  return args.pdb, args.gridres

if __name__ == '__main__':
    pdb, gridres = arg_parse()
    create_template(pdb, gridres)
