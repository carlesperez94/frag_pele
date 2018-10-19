import sys
import os
import ntpath
import prody


def superimpose_backbones(pdb_target, pdb_reference):
    # Loading PDB files
    target = prody.parsePDB(pdb_target)
    reference = prody.parsePDB(pdb_reference)
    # Selection of the backbone
    target_backbone = target.select("backbone")
    reference_backbone = reference.select("backbone")
    # Superimpose the target to the reference backbone
    transformation = prody.calcTransformation(target_backbone, reference_backbone)
    prody.applyTransformation(transformation, target)
    return target, reference


def ligands_rmsd_calculator(pdb_target, pdb_reference, write2report=False, write2pdb=False, ligand_chain="L"):
    """

    :param pdb_target: problem pdb file
    :param pdb_reference: reference pdb file
    :param write2report: if True export results in a file
    :param write2pdb: pdb file with the result of the superposition between "pdb_target" and "pdb_reference"
    :param ligand_chain: name of the chain of the ligand
    :return: superpose the backbone of the pdb_target to the pdb_reference and computes the RMSD of the ligand
    """
    target, reference = superimpose_backbones(pdb_target, pdb_reference)
    target_ligand = target.select("chain {}".format(ligand_chain))
    reference_ligand = reference.select("chain {}".format(ligand_chain))
    # Compute the RMSD of the ligand
    RMSD = prody.calcRMSD(reference_ligand, target_ligand)
    if write2report:
        filename = "{}_to_{}_report.rmsd".format(pdb_target[0:3], pdb_reference[0:4])
        if os.path.exists(filename):
            with open(filename, "a") as report:
                report.write("\nRMSD of {} to  {}:  {}".format(pdb_target, pdb_reference, round(float(RMSD), 3)))
        else:
            with open(filename, "w") as report:
                report.write("RMSD of {} to  {}:  {}".format(pdb_target, pdb_reference, round(float(RMSD), 3)))

    if write2pdb:
        prody.writePDB("{}_to_{}.pdb".format(pdb_target[0:3], pdb_reference[0:4]), target+reference)


def sidechains_rmsd_calculator(pdb_target, pdb_reference, radii, path, write2report=True, ligand_chain="L"):
    """
    :param pdb_target: problem pdb file
    :param pdb_reference: reference pdb file
    :param radii: area that we want to select around the ligand
    :param path: output path
    :param write2report: if true extract a report file
    :param ligand_chain: name of the chain of the ligand
    :return: superpose the backbone of the pdb_target to the pdb_reference and computes the RMSD for each side
    chain in the selection area
    """
    target, reference = superimpose_backbones(pdb_target, pdb_reference)
    selected_area_target = reference.select("protein and (within {} of chain {})".format(radii, ligand_chain))
    unique_residues_target = set(selected_area_target.getResnums())
    list_of_results = []
    for residue_target in unique_residues_target:
        res_selected_target = target.select("protein and resnum {}".format(residue_target))
        res_selected_reference = reference.select("protein and resnum {}".format(residue_target))
        target_CA = target.select("protein and resnum {} and name CA".format(residue_target))
        reference_CA = reference.select("protein and resnum {} and name CA".format(residue_target))
        RMSD = prody.calcRMSD(res_selected_reference, res_selected_target)
        distance_bet_CA = prody.calcRMSD(reference_CA, target_CA)
        residue_information = (residue_target, res_selected_target.getResnames()[0], RMSD, distance_bet_CA)
        list_of_results.append(residue_information)
    list_of_results = sorted(list_of_results)
    pdb_reference_name = os.path.splitext(ntpath.basename(pdb_reference))[0]
    pdb_target_name = os.path.splitext(ntpath.basename(pdb_target))[0]
    if write2report:
        filename = "{}_to_{}_sidech.csv".format(pdb_target_name, pdb_reference_name)
        filename = os.path.join(path, filename)
        with open(filename, "w") as report:
            for result in list_of_results:
                report.write("{:4d}\t{}\t{:5.3f}\t{:5.3f}\t{:5.3f}\n".format(result[0], result[1], float(result[2]),
                             float(result[3]), (float(result[2]) - float(result[3]))))


sidechains_rmsd_calculator("/home/carlespl/project/growing/grow/4DJU_4DJV/selected_result_fragment_4djvC15C5/sel_0_best_structure.pdb",
                           "/home/carlespl/project/growing/grow/4DJU_4DJV/4dju_prepared_w.pdb", 6,
                           "/home/carlespl/project/growing/grow/4DJU_4DJV/selected_result_fragment_4djvC15C5")

