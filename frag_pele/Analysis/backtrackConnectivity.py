import argparse
import os
import sys
import itertools
import lib_prep.pdb_modifier as pm
from AdaptivePELE.utilities import utilities


def parse_arguments():
    """
        Parse the command-line options

        :returns: str, str, str --  path to file to backtrack,
            output path where to write the files, name of the files
    """
    desc = "Adds the connectivity information to a trajectory file that does not have it.\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("pathway", type=str, help="Trajectory file with the backtracked pathway.")
    parser.add_argument("pdb_with_connects", type=str, help="PDB of the final structure prepared. Must contain"
                                                            " CONNECTS section. Remember to not delete water molecules!")
    parser.add_argument("--ligand_chain", type=str, default="L", help="Chain of the ligand.")
    parser.add_argument("-o", type=str, default=None, help="Output file with the new trajectory.")

    args = parser.parse_args()
    return args.pathway, args.pdb_with_connects, args.ligand_chain, args.o


def main(pathway, pdb_with_connects, ligand_chain="L", outfilename=None):
    pdb = pm.PDB(pdb_with_connects, chain=ligand_chain)
    ligand = "".join(pdb.get_atoms_of_chain())
    connects = "".join(pdb.read_conect())
    ligand_with_connects = ligand + connects
    # Temporary file to extract the indexes of the ligand
    with open("ligand_with_connects_tmp.pdb", "w") as out:
        out.write(ligand_with_connects)
    pdb = pm.PDB(in_pdb="ligand_with_connects_tmp.pdb", chain=ligand_chain)
    names_dictionary = pdb.get_names_dictionary()
    os.remove("ligand_with_connects_tmp.pdb")
    pathway_snapshots = utilities.getSnapshots(pathway)
    new_pathway = []
    sys.stderr.write("Rebuilding indexes...\n")
    for snap in pathway_snapshots:
        with open("snap_tmp.pdb", "w") as out:
            out.write(snap)
        pdb = pm.PDB(in_pdb="snap_tmp.pdb", chain=ligand_chain)
        os.remove("snap_tmp.pdb")
        ligand = pdb.get_atoms_of_chain()
        complex_no_ligand = sorted(list(set(pdb.read_atoms_section()) ^ set(ligand)))
        new_ligand = []
        for index, name in names_dictionary.items():
            for line in ligand:
                pdb_name = pm.get_atom_pdb_name_from_line(line).strip()
                if name == pdb_name:
                    new_line = pm.set_index_to_line(line, index)
                    new_ligand.append(new_line)
        new_ligand_with_connects = "".join(new_ligand) + connects
        new_complex = "".join(complex_no_ligand) + new_ligand_with_connects
        new_pathway.append(new_complex)

    sys.stderr.write("Writing connects...\n")
    if not outfilename:
        outfilename = "{}_connected.pdb".format(pathway.split(".pdb")[0])
    if os.path.exists(outfilename):
        os.remove(outfilename)
    for n, model in enumerate(new_pathway):
        with open(outfilename, "a") as output:
            output.write("MODEL {}\n".format(n))
            output.write(model)
            output.write("ENDMDL\n")


if __name__ == "__main__":
    pathway, pdb_with_connects, chain, out = parse_arguments()
    main(pathway, pdb_with_connects, chain, out)

