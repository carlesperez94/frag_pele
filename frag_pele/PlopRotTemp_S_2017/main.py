"""
$Revision: 2.0.0 $
For ligands:
Reads in a maestro mae file and makes a "rotamer enabled" template and the rotamer library to accompany it.  
This consists of finding the backbone core that results in the least number of child bonds rotated with any rotatable bond rotation. 
Reads in the rotatable bonds from a macromodel atomtyping (should be easily modifyiable to read them from stdin).  
Hetgrp_ffgen is used for atomtyping and to determine the geometry in the template fromt the mae.  The mae must only have the residue to atomtype in it. 

Builds rotamer libraries for arbitrary ligand molecules by building and combining rotamer libraries.  There are two general algorithms that are implemented.  This first is using macromodel (or some other external tool) to sample the degrees of freedom and converting the resulting ensemble into a rotamer library.  The second is identifying the rotatable bonds, again using macromodel, and assigning rotamer libraries to these bonds.  For most bonds a simple freely rotatable library ( 0,10,20,30...360), but in the case of closed rings special libraries are built using macromodel sampling.  These component rotamer libraries are then arranged into groups for use in PLOP.  Each group consists of a single tree rooted at the central core.  This core can either be used chosen, or will be chosen based on an algorithm that minimizes the number of bond lengths from the farthest leeaf to the trunk.  Any built rotamer libraries are placed in the current directory and a file called <resname>.rot.assign is also written.  This tells PLOP how to assemble the full combinatoric library that will be used in sidehchain prediction/monte carlo.     

For unnatural amino acids:
Requires: 1) a maestro mae file of an unnatural amino acid with no NMA or ACE termini; the N-H and C=0 are to be left as they are found in a peptide
Outputs: 1) a re-ordered template file (old one is backed up in FILENMAE.hetgrp_ffgen)
         2) a PLOP nonstandard residue specification for pasting into a PLOP control file, both to stdout and to <maefile>_torsions.txt


Options:
   
   --core            Set a core atom 
   --mtor            Maximum number of torsion per group
   --n               Maximum number of entries in the rotamer file
   --ext_charges     Import charges from file e.g. (Inclusion of quantum charges)
   --clean           Clean Intermiadiate files

  
Mae file should be properly atomtyped and a residue name is needed

Most common problem:  As part of this procedure the pdb atom names are often renamed when they are not unique. Alsothis procedure works best if the ligand is minimized.  For these two reasons an atomtyped, minimzed version of the input ligand is written to (input).PlopRotTemp.pdb.  If at all possible, use the ligand structure and atom names from this file in any subsequent plop runs.    

examples:
Build a rotamer library for the following ligand
$SCHRODINGER/utilities/python PlopRotTemp.py 3ert_lig.mae


Make libraries for rotatable bonds in ligand.mae up to a maximum of 4 rotatable bonds in each library
$SCHRODINGER/utilities/python PlopRotTemp.py ligand.mae --mtor=4

For a given ligand named LIG the following files will be created:
ligz               - Template file for use in PLOP, its zmatrix matches the libraries created
LIG.rot.assign     - Summary of all libraries build or used for this ligand read into plop with the command
                     "rot assign all"

---------------------------------------------------------------------------------------------------------



All jobs run on the localhost


USAGE: "$SCHRODINGER/utilities/python main.py [file.mae] --options valueoption"

HELP: $SCHRODINGER/utilities/python main.py --help

"""


import argparse
import sys
import os
import re
import shutil
import PlopRotTemp as pl
from template.templateBuilder import TemplateBuilder



def main(mae_file, max_tors=-1, nrot=1000, user_core_atom=-1, mae_charges = False, clean = False, out_temp=".",
         out_rot=".", gridres="10.0"):

    #Defaults
    nrot = nrot
    max_tors = max_tors
    user_core_atom = user_core_atom
    template_file = ""
    debug = 0  # 1 means don't run exteral commands (assumes output is already there)
    conf_file = 'none';
    output_template_file = ""
    gridres = gridres
    nsamp = 10000
    Ecut = 100
    use_rings = 0
    do_init_min = 0
    max_dist_eq = 0.25
    user_tors = []
    back_tors = []
    back_algorithm = "none"
    back_conf_file = ""
    hetgrp_opt = ""
    OPLS = "2005"
    user_fixed_bonds = []
    files2clean = []
    use_mult_lib = 1
    run_conf = 0
    algorithm = "MCMM"
    gridres_oh = gridres
    unnat_res = 0  # for old-style PLOP nonstandard side chain
    resno = 1  #for old-style PLOP nonstandard side chain
    chain = 'A'  #for old-style PLOP nonstandard side chain
    grow = 0
    tree = 0  # for old-style PLOP nonstandard ligand tree-style torsion reordering
    R_group_root_atom_name = 'None'  # which atom do you want to start sampling at?

    if (unnat_res == 1):
        init_min = 0  #the input mae file is for a peptide and will not have a suitable Lewis structure
        if (template_file == ""):
            print("Cannot use unnatural residues without pre-made template files!")
            sys.exit(-1)
        use_mult_lib = 1  # so a dummy conformational search is performed just to see which bonds are rotatble
        use_rings = 0  #for now; I'm not sure that ring torsions will follow the tree pattern appropriately, this would need testing
        #For now, just try low energy ring conformations


    atomnames = pl.find_names_in_mae(mae_file)
    
    pl.check_repite_names(atomnames)


    root = pl.get_root_path(mae_file)

    print("\n")
    print("INPUT")
    print("mae_file {}".format(mae_file))
    print("root {}".format(root))
    print("OPLS {}".format(OPLS))
    print("hetgrp options '{}'".format(hetgrp_opt))
    print("User template file '{}'".format(template_file))
    print("User output template file '{}'".format(output_template_file))
    print("\n")


    #Build a template file

    print("TEMPLATE GENERATION")

    resname=pl.find_resnames_in_mae(mae_file)
    template_output = os.path.join(out_temp, "{}z".format(resname[0].lower()))

    template_builder = TemplateBuilder(mae_file, template_output)

    [template_file, output_template_file, mae_file_hetgrp_ffgen, files, resname] = \
        template_builder.build_template(mae_charges)

    print(output_template_file)



    

    print("\n")
    if (do_init_min == 1):
        mcu_mini = mu.ComUtil(ffld='opls2005', serial=True, solv=True, nant=False, demx=True)
        mcu_mini.SOLV[2] = 1  # water
        mini_root = root + "_mini"
        com_file = mcu_mini.mini(mae_file_hetgrp_ffgen, mini_root + '.com')
        print('\nMINIMIZATION\nRunning minimization: {0} -> {1} -out.mae\n'.format(mae_file_hetgrp_ffgen, mini_root))
        if (not debug):
            cmd = mcu_mini.getLaunchCommand(com_file)
            job = jc.launch_job(cmd)
            job.wait()
            files2clean.append(mini_root + '-out.mae')
            files2clean.append(mini_root + '.log')
    #       files2clean.append(mini_root + '-out.tmp')
    #        files2clean.append(mini_root + '.com')
        mae_min_file = mini_root + "-out.mae"
    else:
        print('\nSkipping Minimization\n ')
        mae_min_file = mae_file_hetgrp_ffgen


    print("\n")

    if (unnat_res == 1):
        [mae_num, parent, rank, tors, use_rings, group, tors_ring_num] = \
            pl.FindCoreAA(mae_min_file, user_fixed_bonds, use_rings, resname, use_mult_lib, user_core_atom, user_tors)
        tors_ring_num = []
        for t in tors: tors_ring_num.append(0);
    else:
        print('FINDING CORE')
        if (grow == 1 and user_core_atom == -1): user_core_atom = -2
        #######Assign_rank--> Extremely slow!
        [mae_num, parent, rank, tors, use_rings, group, back_tors, tors_ring_num] = \
            pl.FindCore(mae_min_file, user_fixed_bonds, use_rings, resname, \
                     use_mult_lib, user_core_atom, user_tors, back_tors, max_tors, R_group_root_atom_name)
    
    if (use_rings == 1):
        print("Found flexible rings")

    newtors = []
    if (unnat_res == 1 or grow == 1 ):
        newtors = pl.ReorderTorsionsAA(tors, mae_num)


    #Change from mae files atom numbering to the template file ones
    #Convert Torsions to match the template file atom numbering
    #Ring numbers don't have to be changed
    [mae2temp, temp2mae] = pl.MatchTempMaeAtoms(mae_min_file, template_file)

    old_atom_num = [];
    new_tors = [];
    new_back_tors = [];
    for i in mae_num:
        old_atom_num.append(-100)
    for i in range(len(mae2temp)):
        old_atom_num[i] = mae2temp[mae_num[i]]
    for i in range(len(tors)):
        temp = [mae2temp[tors[i][0]], mae2temp[tors[i][1]]]
        new_tors.append(temp)
    for i in range(len(back_tors)):
        temp = [mae2temp[back_tors[i][0]], mae2temp[back_tors[i][1]]]
        new_back_tors.append(temp)
    tors = []
    for i in range(len(new_tors)):
        temp = [old_atom_num.index(new_tors[i][0]), old_atom_num.index(new_tors[i][1])]
        tors.append(temp)
    back_tors = []
    for i in range(len(new_back_tors)):
        temp = [old_atom_num.index(new_back_tors[i][0]), old_atom_num.index(new_back_tors[i][1])]
        back_tors.append(temp)


    #Make (or read) original tempalte file
    print('\n')
    print('CREATE ROTAMER TEMPLATE FILE: {}'.format(output_template_file))
    names = pl.ReorderTemplate(old_atom_num, parent, rank, template_file, output_template_file, mae_file, 
        R_group_root_atom_name=R_group_root_atom_name)

    [tors, tors_ring_num, zmat_atoms] = pl.FindTorsAtom(tors, tors_ring_num, parent)
    #Eliminate Torsions in the backbone (included when entire rings are included in the torsions)
    [tors, tors_ring_num, zmat_atoms] = pl.EliminateBackboneTors(tors, tors_ring_num, zmat_atoms, rank)

    pl.replace_vdwr_from_library(output_template_file)

    if (unnat_res == 1 or grow == 1):
        mynonstandard = pl.TetherRotBonds(mae_file, chain, resno, log_file, newtors)
        mynonstandard.output_rotbonds(R_group_root_atom_name=R_group_root_atom_name)


    else:

        # Reorder the torsions
        for i in range(len(tors)):
            tors[i].sort()
        for i in range(len(tors)):
            for j in range(i + 1, len(tors)):
                if (tors[i] > tors[j]):
                    temp = tors[i];
                    tors[i] = tors[j];
                    tors[j] = temp;
                    temp = tors_ring_num[i];
                    tors_ring_num[i] = tors_ring_num[j];
                    tors_ring_num[j] = temp

        [tors, tors_ring_num, zmat_atoms] = pl.FindTorsAtom(tors, tors_ring_num, parent)


    ################################CHANGE MACROMODEL CONFORMATIONAL SEARCH######################
    #Run the conformational Search
    conf_root = root + "_conf"
    if (conf_file == conf_root + '-out.mae'):
        raise Exception('Must use different name for conformational file')

    if (conf_file == ''):
        run_conf = 1
    else:
        run_conf = 0

    back_lib = "";
    if (unnat_res != 1):  
        if (conf_file != 'none'):
            rotamers_file = pl.make_libraries(resname, conf_file, root, names, zmat_atoms, group, use_rings, use_mult_lib,
                           output_template_file, gridres, debug, out_rot)
            print("\n")
            print("CREATE ROTAMER LIBRARY")
            print(rotamers_file)
            print("\n")


        
        else:
            if (len(zmat_atoms) > 0):
                ring_libs = pl.build_ring_libs(mae_min_file, root, resname, tors, \
                                            tors_ring_num, names, rank, parent, old_atom_num, mae2temp, gridres,
                                            files2clean, debug)
     
            else:
                ring_libs = []
                print("No rotatable sidechains found")
            
            rotamers_file = pl.find_build_lib(resname, mae_min_file, root, tors, names, group, gridres, gridres_oh, use_rings, back_lib,
                                  tors_ring_num, ring_libs, debug, out_rot)
            print("\n")
            print("CREATE ROTAMER LIBRARY")
            print(rotamers_file)
            print("\n")




    if (clean):
        for file in files2clean:
            print('Removing Intermediate File: {}'.format(file))
            os.remove(file)

    return output_template_file, rotamers_file



def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("mae_file", type=str, help="ligand maestro mae file")
    parser.add_argument("--core", type=int, help="Give one atom of the core section", default=-1)
    parser.add_argument("--mtor", type=int, help="Gives the maximum number of torsions allowed in each group.  Will freeze bonds to extend the core if necessary.", default=-1)
    parser.add_argument("--n", type=int, help="Maximum Number of Entries in Rotamer File", default=1000)
    parser.add_argument("--mae_charges", help="Use charges in mae", action='store_true')
    parser.add_argument("--clean", help="Whether to clean up all the intermediate files", action='store_true')
    parser.add_argument("--gridres", help="Rotamer resolution", type=str, default="10.0")
    parser.add_argument("--out_temp", type=str, help="Template output path")
    parser.add_argument("--out_rot", type=str, help="Rotamer ouput path")
    args = parser.parse_args()

    
    return args.mae_file, args.mtor, args.n, args.core, args.mae_charges, args.clean, args.out_temp, args.out_rot, args.gridres

if __name__ == "__main__":
    mae_file, mtor, n, core, mae_charge, clean, out_temp, out_rot, gridres = parse_args() 
    template, rotamers_file = main(mae_file, mtor, n, core, mae_charge, clean, out_temp, out_rot, gridres=gridres)
    

    print("########################################################################")
    print("\n{} template and {} rotamer library has been successfully created in {}\n").format(
        template,rotamers_file, os.getcwd())
    print("########################################################################")


    

