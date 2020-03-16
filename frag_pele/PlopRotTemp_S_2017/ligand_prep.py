import os
import sys
from schrodinger import structure as st
import subprocess
import argparse
# Local import
frag_pele_dirname = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(frag_pele_dirname)
import frag_pele.PlopRotTemp_S_2017.main as plop


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


def create_template(pdb, gridres, out_temp, out_lib):
   mae_file = convert_mae(pdb)
   plop.main(mae_file, out_temp=out_temp, out_rot=out_lib, gridres=gridres)
   os.remove(mae_file)


def arg_parse():
  parser = argparse.ArgumentParser()

  parser.add_argument("pdb", type=str, help="ligand file to templatize")
  parser.add_argument("gridres", type=str, help="Degrees of rotation.")
  parser.add_argument("out_temp", type=str, help="Output template path.")
  parser.add_argument("out_lib", type=str, help="Output rotamer library path.")

  args = parser.parse_args()

  return args.pdb, args.gridres, args.out_temp, args.out_lib


if __name__ == '__main__':
    pdb, gridres, outtemp, outlib = arg_parse()
    create_template(pdb, gridres, outtemp, outlib)
