import os
import time
import shutil
import subprocess
import frag_pele.constants as c

DATA_PATH = "Data/Templates/OPLS2005/Protein"

def create_template(pdb_file, sch_path=c.SCHRODINGER,
                    plop_script_path="../PlopRotTemp_S_2017/ligand_prep.py", rotamers="10.0",
                    out_templates_path=".", path_to_lib=".", cov_res=None, work_dir="."):
    if "/" in pdb_file:
        filename = pdb_file.split("/")[-1]
        filename = filename.split(".pdb")[0].lower()
    else:
        filename = pdb_file.split(".pdb")[0].lower()

    if cov_res:
       if os.path.exists(os.path.join(work_dir, DATA_PATH, filename)):
           print("Template found in {}".format(os.path.join(DATA_PATH, filename)))
           shutil.copy(os.path.join(work_dir, DATA_PATH, filename), 
                       os.path.join(work_dir, "DataLocal/Templates/OPLS2005/Protein/templates_generated"))
           return filename
       else:
           schrodinger_path = sch_path
           sch_python = os.path.join(sch_path, "utilities/python")
           prepare_pdb(pdb_file, "GRW.pdb", schrodinger_path)
           cmd = "{} {} {} {} {} {}".format(sch_python, plop_script_path, "GRW.pdb", rotamers,
                                            os.path.join(work_dir, 
                                            "DataLocal/Templates/OPLS2005/Protein/templates_generated"), 
                                            path_to_lib)

           print("Running PlopRotTemp...")
           print(cmd)
           try:
              subprocess.call(cmd.split())
           except OSError:
               raise OSError("Path {} not foud. Change schrodinger path under frag_pele/constants.py".format(sch_python))
    else:
        cmd = "{} {} {} {} {} {}".format(sch_python, plop_script_path, pdb_file, rotamers,
                                            out_templates_path, path_to_lib)
        print("Running PlopRotTemp...")
        print(cmd)
        try:
            subprocess.call(cmd.split())
        except OSError:
            raise OSError("Path {} not foud. Change schrodinger path under frag_pele/constants.py".format(sch_python))
    return filename    

def prepare_pdb(pdb_in, pdb_out, sch_path):
    command = [os.path.join(sch_path, "utilities/prepwizard"), pdb_in, pdb_out, "-noepik", "-noprotassign",
               "-noimpref", "-noccd", "-NOJOBID"]
    print(command)
    subprocess.call(command)
