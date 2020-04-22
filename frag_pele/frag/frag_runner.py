# Python Imports
import os
import time
import glob
import shutil
import subprocess
import logging

# Third-Party Imports

# Project Imports
from frag_pele.Growing import template_fragmenter, simulations_linker
from frag_pele.Growing import add_fragment_from_pdbs, bestStructs
from frag_pele.Analysis import analyser
from frag_pele.Helpers import helpers, center_of_mass
from frag_pele.Helpers import clusterizer, folder_handler, constraints, check_constants
import frag_pele.constants as c
from frag_pele.Banner import Detector
from frag_pele.frag.ClusterParameters.ClusterParameters import ClusterParameters
from frag_pele.frag.PeleParameters.pele_parameters import PeleParameters

FilePath = os.path.abspath(__file__)
PackagePath = os.path.dirname(FilePath)

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

# Get current path
curr_dir = os.path.abspath(os.path.curdir)


class Frag:

    def __init__(self):
        self._initial_time = 0
        # todo put here peleParameters

    def run_frag(self, pele_parameters: PeleParameters, cluster_parameters: ClusterParameters, complex_pdb,
                 fragment_pdb, core_atom, fragment_atom, growing_steps, criteria,
                 plop_path, sch_python, restart, ID,
                 h_core=None, h_frag=None, c_chain="L", f_chain="L", rotamers="30.0",
                 mae=False, rename=False, threshold_clash=1.7, explorative=False,
                 sampling_control=None, only_prepare=False, only_grow=False,
                 no_check=False):
        """
        Description: FrAG is a Fragment-based ligand growing software which performs automatically the addition of several
        fragments to a core structure of the ligand in a protein-ligand complex.
        # todo parameters:

        :param growing_steps: 
        :param cluster_parameters:
        :param pele_parameters:
        :param complex_pdb: Path to the PDB file which must contain a protein-ligand complex. Its ligand will be
        used as the core structure. Remember to rename the ligand chain with a different character in
        order to detect it.
        :type complex_pdb: str
        :param fragment_pdb: Path to the PDB file containing the fragment.
        :type fragment_pdb: str
        :param core_atom: PDB-atom-name of the atom of the core that will be used as starting point to grow the fragment.
        :type core_atom: str (max length: 4)
        :param fragment_atom: PDB-atom-name of the atom of the fragment that will bond the core.
        :type fragment_atom: str (max length: 4)
        :param criteria: Name of the column of the report file to select the structures that will spawn in the
        next GS. Additionally, this parameter will be the selection criteria to extract the best
        structure after completing the growing.
        :type criteria: str
        :param plop_path: Absolute path to PlopRotTemp.py.
        :type plop_path: str
        :param sch_python: Absolute path to Schrodinger's python.
        :type sch_python: str
        :param restart: If set FrAG will find your last iteration completed and will restart from this point.
        :type restart: bool
        :param ID: Name to identify the new ligand in certain folder.
        :type ID: str
        :param h_core: PDB-atom-name of the hydrogen bonded to the selected heavy atom of the core that will be replaced for
         the fragment, being the initial point of the growing.
        :type h_core: str (max length: 4)
        :param h_frag: PDB-atom-name of the hydrogen bonded to the selected heavy atom of the fragment that will removed to
        bond the fragment with the core.
        :type h_frag: str (max length: 4)
        :param c_chain: Chain name for the ligand in the complex_pdb.
        :type c_chain: str (max length: 1)
        :param f_chain: Chain name for the ligand in the fragment_pdb.
        :type f_chain: str (max length: 1)
        :param mae: If set, the output structures will be saved in "mae." format.
        :type mae: bool
        :param rename: If set, the pdb atom names of the fragment will be renamed with "G+atom_number".
        :type rename: bool
        :param threshold_clash: Distance that will be used to detect contacts between atoms of the fragment and atoms
        of the core.
        :type threshold_clash: float
        :param sampling_control: templatized control file to be used in the sampling simulation.
        :type sampling_control: str
        :return:

        """
        # Extract pele parameters
        pele_params_path, pele_params_archives, pele_params_sim_values = pele_parameters.extract_parameters()

        # Check harcoded path in constants.py
        self._check_hardcoded_path_in_constants(no_check)

        # Time computations
        self._start_timer()

        # Global variable to keep info
        simulation_info = []  # todo: not used, delete it

        # Path definition
        plop_relative_path = os.path.join(PackagePath, plop_path)  # todo: extract
        current_path = os.path.abspath(".")

        pdb_basename = complex_pdb.split(".pdb")[0]  # Get the name of the pdb without extension
        if "/" in pdb_basename:
            pdb_basename = pdb_basename.split("/")[-1]  # And if it is a path, get only the name

        working_dir = os.path.join(current_path, "{}_{}".format(pdb_basename, ID))
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)  # Creating a working directory for each PDB-fragment combination
        pdbout_folder = os.path.join(working_dir, cluster_parameters.pdbout)
        path_to_templates_generated = os.path.join(working_dir,
                                                   "DataLocal/Templates/OPLS2005/HeteroAtoms/templates_generated")
        path_to_templates = os.path.join(working_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms")
        path_to_lib = os.path.join(working_dir, "DataLocal/LigandRotamerLibs")

        # Creation of output folder
        folder_handler.check_and_create_DataLocal(working_dir=working_dir)

        # Creating constraints
        const = "\n".join(constraints.retrieve_constraints(complex_pdb, {}, {}, 5, 5, 10))

        # Creating symbolic links
        helpers.create_symlinks(pele_params_path.data, os.path.join(working_dir, 'Data'))
        helpers.create_symlinks(pele_params_path.documents, os.path.join(working_dir, 'Documents'))

        #  ---------------------------------------Pre-growing part - PREPARATION -------------------------------------------
        fragment_names_dict, hydrogen_atoms, pdb_to_initial_template, pdb_to_final_template, pdb_initialize, \
        core_original_atom, fragment_original_atom = add_fragment_from_pdbs.main(complex_pdb, fragment_pdb, core_atom,
                                                                                 fragment_atom, growing_steps, h_core=h_core,
                                                                                 h_frag=h_frag, core_chain=c_chain,
                                                                                 fragment_chain=f_chain, rename=rename,
                                                                                 threshold_clash=threshold_clash,
                                                                                 output_path=working_dir,
                                                                                 only_grow=only_grow)

        # Create the templates for the initial and final structures
        template_resnames = []
        for pdb_to_template in [pdb_to_initial_template, pdb_to_final_template]:
            if not only_grow:
                cmd = "{} {} {} {} {} {}".format(sch_python, plop_relative_path, os.path.join(working_dir,
                                                                                              add_fragment_from_pdbs.c.PRE_WORKING_DIR,
                                                                                              pdb_to_template), rotamers,
                                                 path_to_templates_generated, path_to_lib)

                try:
                    subprocess.call(cmd.split())
                except OSError:
                    raise OSError(
                        "Path {} not foud. Change schrodinger path under frag_pele/constants.py".format(sch_python))
            template_resname = add_fragment_from_pdbs.extract_heteroatoms_pdbs(
                os.path.join(working_dir, add_fragment_from_pdbs.
                             c.PRE_WORKING_DIR, pdb_to_template),
                False, c_chain, f_chain)
            template_resnames.append(template_resname)

        # Set box center from ligand COM
        resname_core = template_resnames[0]
        center = center_of_mass.center_of_mass(os.path.join(working_dir, c.PRE_WORKING_DIR, "{}.pdb".format(resname_core)))

        # Get template filenames
        template_initial, template_final = ["{}z".format(resname.lower()) for resname in template_resnames]

        if only_prepare:
            print("Files of {} prepared".format(ID))
            return

        # --------------------------------------------GROWING SECTION-------------------------------------------------------
        # Lists definitions

        templates = ["{}_{}".format(os.path.join(path_to_templates_generated, template_final), n) for n in
                     range(0, growing_steps + 1)]

        results = [os.path.join(working_dir, c.OUTPUT_FOLDER, str(n)) for n in range(0, growing_steps + 1)]

        pdbs = [pdb_initialize if n == 0 else "{}_{}".format(n, pdb_initialize) for n in range(0, growing_steps + 1)]

        pdb_selected_names = ["initial_0_{}.pdb".format(n) for n in range(0, pele_params_sim_values.cpus - 1)]

        # Generate starting templates
        template_fragmenter.main(template_initial_path=os.path.join(path_to_templates_generated, template_initial),
                                 template_grown_path=os.path.join(path_to_templates_generated, template_final),
                                 step=1, total_steps=growing_steps, hydrogen_to_replace=core_original_atom,
                                 core_atom_linker=core_atom,
                                 tmpl_out_path=os.path.join(path_to_templates_generated, "{}_0".format(template_final)))

        # Make a copy in the main folder of Templates in order to use it as template for the simulation
        shutil.copy(os.path.join(path_to_templates_generated, "{}_0".format(template_final)),
                    os.path.join(path_to_templates, template_final))  # Replace the original template in the folder

        # Clear PDBs folder
        if not restart:
            list_of_subfolders = glob.glob("{}*".format(pdbout_folder))
            for subfolder in list_of_subfolders:
                shutil.rmtree(subfolder)

        # Simulation loop - LOOP CORE
        for i, (template, pdb_file, result) in enumerate(zip(templates, pdbs, results)):

            # Only if reset
            if restart:
                if os.path.exists(os.path.join(pdbout_folder, "{}".format(i))) and os.path.exists(
                        os.path.join(pdbout_folder,
                                     "{}".format(i),
                                     "initial_0_0.pdb")):
                    print("STEP {} ALREADY DONE, JUMPING TO THE NEXT STEP...".format(i))
                    continue
            # Otherwise start from the beggining
            pdb_input_paths = ["{}".format(os.path.join(pdbout_folder, str(i - 1), pdb_file)) for pdb_file in
                               pdb_selected_names]
            # Banned dihedrals will be checked here
            if cluster_parameters.banned:
                pdbs_with_banned_dihedrals = Detector.check_folder(folder=os.path.join(pdbout_folder, str(i - 1)),
                                                                   threshold=cluster_parameters.limit,
                                                                   dihedrals=cluster_parameters.banned,
                                                                   lig_chain=c_chain,
                                                                   processors=pele_params_sim_values.cpus)
                pdb_input_paths = [pdb_file for pdb_file, flag in pdbs_with_banned_dihedrals.items() if flag]

            # Control file modification
            overlapping_factor = float(pele_params_sim_values.min_overlap) + (((float(pele_params_sim_values.max_overlap)
                                                                                - float(pele_params_sim_values.min_overlap)) * i) / growing_steps)
            overlapping_factor = "{0:.2f}".format(overlapping_factor)

            if i != 0:
                # Check atom overlapping
                pdbs_with_overlapping = clusterizer.check_atom_overlapping(pdb_input_paths)
                pdb_input_paths_checked = []
                for pdb in pdb_input_paths:
                    if pdb not in pdbs_with_overlapping:
                        pdb_input_paths_checked.append(pdb)
                simulation_file = simulations_linker.control_file_modifier(pele_params_archives.control_file, pdb=pdb_input_paths_checked, step=i,
                                                                           license=pele_params_path.pele_license,
                                                                           working_dir=working_dir,
                                                                           overlap=overlapping_factor, results_path=result,
                                                                           steps=pele_params_sim_values.steps,
                                                                           chain=c_chain, constraints=const, center=center,
                                                                           temperature=pele_params_sim_values.temperature, seed=pele_params_sim_values.seed,
                                                                           steering=pele_params_sim_values.steering,
                                                                           translation_high=pele_params_sim_values.translation_high,
                                                                           translation_low=pele_params_sim_values.translation_low,
                                                                           rotation_high=pele_params_sim_values.rotation_high,
                                                                           rotation_low=pele_params_sim_values.rotation_low,
                                                                           radius=pele_params_sim_values.radius_box)
            else:
                logger.info(c.SELECTED_MESSAGE.format(pele_params_archives.control_file, pdb_initialize, result, i))
                simulation_file = simulations_linker.control_file_modifier(pele_params_archives.control_file, pdb=[pdb_initialize], step=i,
                                                                           license=pele_params_path.pele_license,
                                                                           working_dir=working_dir,
                                                                           overlap=overlapping_factor, results_path=result,
                                                                           steps=pele_params_sim_values.steps,
                                                                           chain=c_chain, constraints=const, center=center,
                                                                           temperature=pele_params_sim_values.temperature, seed=pele_params_sim_values.seed,
                                                                           steering=pele_params_sim_values.steering,
                                                                           translation_high=pele_params_sim_values.translation_high,
                                                                           translation_low=pele_params_sim_values.translation_low,
                                                                           rotation_high=pele_params_sim_values.rotation_high,
                                                                           rotation_low=pele_params_sim_values.rotation_low,
                                                                           radius=pele_params_sim_values.radius_box)

            logger.info(c.LINES_MESSAGE)
            if i != 0:
                template_fragmenter.main(template_initial_path=os.path.join(path_to_templates_generated, template_initial),
                                         template_grown_path=os.path.join(path_to_templates_generated, template_final),
                                         step=i + 1, total_steps=growing_steps, hydrogen_to_replace=core_original_atom,
                                         core_atom_linker=core_atom,
                                         tmpl_out_path=os.path.join(path_to_templates, template_final))

            # Make a copy of the template file in growing_templates folder
            shutil.copy(os.path.join(path_to_templates, template_final), template)

            # Creating results folder
            folder_handler.check_and_create_results_folder(result, working_dir)
            # ------SIMULATION PART------
            # Change directory to the working one
            os.chdir(working_dir)
            simulations_linker.simulation_runner(pele_params_path.pele_dir, simulation_file, pele_params_sim_values.cpus)
            logger.info(c.LINES_MESSAGE)
            logger.info(c.FINISH_SIM_MESSAGE.format(result))
            # Before selecting a step from a trajectory we will save the input PDB file in a folder
            folder_handler.check_and_create_pdb_clusters_folder(pdbout_folder, i)

            # ---------------------------------------------------CLUSTERING-------------------------------------------------
            # Transform column name of the criteria to column number
            result_abs = os.path.abspath(result)
            logger.info("Looking structures to cluster in '{}'".format(result_abs))
            column_number = clusterizer.get_column_num(result_abs, criteria, pele_params_archives.report)
            # Selection of the trajectory used as new input
            clusterizer.cluster_traject(str(template_resnames[1]), pele_params_sim_values.cpus - 1, column_number, cluster_parameters.distance_contact,
                                        cluster_parameters.clusterThreshold, "{}*".format(os.path.join(result_abs, pele_params_archives.traject)),
                                        os.path.join(pdbout_folder, str(i)), os.path.join(result_abs),
                                        cluster_parameters.epsilon, pele_params_archives.report, cluster_parameters.condition, cluster_parameters.metricweights, cluster_parameters.nclusters)
        # ----------------------------------------------------EQUILIBRATION-------------------------------------------------
        # Set input PDBs
        pdb_inputs = ["{}".format(os.path.join(pdbout_folder, str(growing_steps), pdb_file)) for pdb_file in
                      pdb_selected_names]
        if cluster_parameters.banned:
            pdbs_with_banned_dihedrals = Detector.check_folder(folder=os.path.join(pdbout_folder, str(growing_steps)),
                                                               threshold=cluster_parameters.limit,
                                                               dihedrals=cluster_parameters.banned,
                                                               lig_chain=c_chain,
                                                               processors=pele_params_sim_values.cpus)
            pdb_inputs = [pdb_file for pdb_file, flag in pdbs_with_banned_dihedrals.items() if flag]
        if not os.path.exists(os.path.join(working_dir, "sampling_result")):  # Create the folder if it does not exist
            os.mkdir(os.path.join(working_dir, "sampling_result"))
        # Modify the control file to increase the steps TO THE SAMPLING SIMULATION
        if sampling_control:
            simulation_file = simulations_linker.control_file_modifier(sampling_control, pdb=pdb_inputs, step=growing_steps,
                                                                       license=pele_params_path.pele_license, working_dir=working_dir,
                                                                       overlap=pele_params_sim_values.max_overlap,
                                                                       results_path=os.path.join(working_dir,
                                                                                                 "sampling_result"),
                                                                       steps=pele_params_sim_values.pele_eq_steps, chain=c_chain,
                                                                       constraints=const, center=center,
                                                                       temperature=pele_params_sim_values.temperature, seed=pele_params_sim_values.seed,
                                                                       steering=pele_params_sim_values.steering,
                                                                       translation_high=pele_params_sim_values.translation_high,
                                                                       translation_low=pele_params_sim_values.translation_low,
                                                                       rotation_high=pele_params_sim_values.rotation_high,
                                                                       rotation_low=pele_params_sim_values.rotation_low,
                                                                       radius=pele_params_sim_values.radius_box)
        elif explorative and not sampling_control:
            simulation_file = simulations_linker.control_file_modifier(pele_params_archives.control_file, pdb=pdb_inputs, license=pele_params_path.pele_license,
                                                                       working_dir=working_dir, step=growing_steps,
                                                                       overlap=pele_params_sim_values.max_overlap,
                                                                       results_path=os.path.join(working_dir,
                                                                                                 "sampling_result"),
                                                                       steps=pele_params_sim_values.pele_eq_steps,
                                                                       chain=c_chain, constraints=const, center=center,
                                                                       temperature=pele_params_sim_values.temperature, seed=pele_params_sim_values.seed,
                                                                       steering=2,
                                                                       translation_high=0.5,
                                                                       translation_low=0.3,
                                                                       rotation_high=0.4,
                                                                       rotation_low=0.15,
                                                                       radius=25)
        else:
            simulation_file = simulations_linker.control_file_modifier(pele_params_archives.control_file, pdb=pdb_inputs, step=growing_steps,
                                                                       license=pele_params_path.pele_license, overlap=pele_params_sim_values.max_overlap,
                                                                       working_dir=working_dir,
                                                                       results_path=os.path.join(working_dir,
                                                                                                 "sampling_result"),
                                                                       steps=pele_params_sim_values.pele_eq_steps, chain=c_chain,
                                                                       constraints=const, center=center,
                                                                       temperature=pele_params_sim_values.temperature, seed=pele_params_sim_values.seed,
                                                                       steering=pele_params_sim_values.steering,
                                                                       translation_high=pele_params_sim_values.translation_high,
                                                                       translation_low=pele_params_sim_values.translation_low,
                                                                       rotation_high=pele_params_sim_values.rotation_high,
                                                                       rotation_low=pele_params_sim_values.rotation_low,
                                                                       radius=pele_params_sim_values.radius_box)

        # EQUILIBRATION SIMULATION
        if not (restart and os.path.exists("selected_result")):
            shutil.copy(os.path.join(path_to_templates_generated, template_final), path_to_templates)
            logger.info(".....STARTING EQUILIBRATION.....")
            simulations_linker.simulation_runner(pele_params_path.pele_dir, simulation_file, pele_params_sim_values.cpus)
        os.chdir(curr_dir)
        equilibration_path = os.path.join(working_dir, "sampling_result")
        # SELECTION OF BEST STRUCTURES
        selected_results_path = os.path.join(working_dir, "top_result")
        if not os.path.exists(selected_results_path):  # Create the folder if it does not exist
            os.mkdir(selected_results_path)
        best_structure_file, all_output_files = bestStructs.main(criteria, selected_results_path, path=equilibration_path,
                                                                 n_structs=50)

        shutil.copy(os.path.join(selected_results_path, best_structure_file), os.path.join(working_dir, c.PRE_WORKING_DIR,
                                                                                           selected_results_path + ".pdb"))
        # COMPUTE AND SAVE THE SCORE
        analyser.analyse_at_epoch(report_prefix=pele_params_archives.report, path_to_equilibration=equilibration_path, execution_dir=curr_dir,
                                  column=criteria, quantile_value=0.25)

        # MOVE FROM PDB TO MAE
        if mae:
            if sch_python.endswith("python"):
                schrodinger_path = os.path.dirname(os.path.dirname(sch_python))
            elif sch_python.endswith("run"):
                schrodinger_path = os.path.dirname(sch_python)
            python_file = os.path.join(os.path.dirname(FilePath), "Analysis/output_files.py")
            for outputfile in all_output_files:
                filename = os.path.join(selected_results_path, outputfile)
                command = "{} {} {} --schr {} {}".format(sch_python, python_file, filename, schrodinger_path, "--remove")
                subprocess.call(command.split())

        # COMPUTE TIME
        end_time = time.time()
        total_time = (end_time - self._initial_time) / 60
        logging.info("Growing of {} in {} min".format(fragment_pdb, total_time))

        return fragment_names_dict

    def _check_hardcoded_path_in_constants(self, no_check):
        if not no_check:
            check_constants.check()

    def _start_timer(self):
        self._initial_time = time.time()

    # todo: WIP
    def _define_paths(self, plop_path, complex_pdb, pdbout, ID):
        # Path definition
        plop_relative_path = os.path.join(PackagePath, plop_path)  # todo: extract
        current_path = os.path.abspath(".")

        pdb_basename = complex_pdb.split(".pdb")[0]  # Get the name of the pdb without extension
        if "/" in pdb_basename:
            pdb_basename = pdb_basename.split("/")[-1]  # And if it is a path, get only the name

        working_dir = os.path.join(current_path, "{}_{}".format(pdb_basename, ID))
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)  # Creating a working directory for each PDB-fragment combination
        pdbout_folder = os.path.join(working_dir, pdbout)
        path_to_templates_generated = os.path.join(working_dir,
                                                   "DataLocal/Templates/OPLS2005/HeteroAtoms/templates_generated")
        path_to_templates = os.path.join(working_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms")
        path_to_lib = os.path.join(working_dir, "DataLocal/LigandRotamerLibs")
