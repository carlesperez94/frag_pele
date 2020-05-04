import pytest
import shutil
import glob
import subprocess
import os
import frag_pele.main as mn

FLAGS = [
'"radius" : 3,',
'"steeringUpdateFrequency": 3,',
'translationRange": 33',
'translationRange": 31',
'"overlapFactor": 3',
'"temperature": 100001',
'"seed" : 3',
'"overlapFactor": 3.00'
]


def test_double():
    mn.main("1w7h_preparation_structure_2w.pdb", "serie_file_carbonil.conf", cpus=4, iterations=1, steps=1, temperature=1500, pele_eq_steps=1, contrl="test.conf", test=True)
    assert glob.glob("1w7h_preparation_structure_2w_carbonylC6-H7O1-C1/top_result/*BindingEnergy*.pdb")

def test_HT():
    mn.main("1w7h_preparation_structure_2w.pdb", "serie_file_carbonil.conf", cpus=4, iterations=1, steps=1, pele_eq_steps=1, contrl="test.conf", protocol="HT", mae=True, test=True)
    assert glob.glob("1w7h_preparation_structure_2w_aminoC1N1/top_result/*traject*.mae")

def test_sequential():
    mn.main("1w7h_preparation_structure_2w.pdb", "sequential.conf", cpus=4, iterations=1, steps=1, pele_eq_steps=1, contrl="test.conf", test=True)
    assert glob.glob("1w7h_preparation_structure_2w_aminoC1N1/top_result/*BindingEnergy*.pdb") and glob.glob("top_result_aminoC1N1phenylC6C1/top_result/*BindingEnergy*.pdb")

def test_flags():
    mn.main("1w7h_preparation_structure_2w2.pdb", "serie_file2.conf", cpus=4, iterations=1, steps=1, pele_eq_steps=1, temperature=100001, criteria="currentEnergy", 
    rotamers=60, resfold="bests", seed=3, steering=33, contrl="test.conf", test=True, debug=True,
    translation_high=33, translation_low=31, rotation_high=33, rotation_low=31,
    radius_box=3, pdbout="best.pdb", min_overlap=3, max_overlap=3, c_chain="Z", f_chain="Z")
    # Check flags
    with open("1w7h_preparation_structure_2w2_phenyl2C6C1/control_template.conf", "r") as f:
        lines = f.readlines()
        errors = [flag for flag in FLAGS if flag not in lines]
        assert errors
    with open("1w7h_preparation_structure_2w2_phenyl2C6C1/DataLocal/LigandRotamerLibs/3IP.rot.assign", "r") as f:
        lines = f.readlines()
        if "60" in lines[1]:
            assert True
        else:
            assert False

def test_restart(directory="1w7h_preparation_structure_2w_carbonylC6-H7O1-C1"):
    os.chdir("data")
    #make clean
    if os.path.exists(directory):
        shutil.rmtree(directory)
    #create working folder
    shutil.copytree("original", directory)

    # Test restart if simulation is done
    mn.main("1w7h_preparation_structure_2w.pdb", "serie_file_carbonil.conf", cpus=4, iterations=1, steps=1, temperature=1500, pele_eq_steps=1, contrl="test.conf", test=True, restart=True)
    assert glob.glob(os.path.join(directory, "top_result/*BindingEnergy*.pdb"))

    # Test restart if equilibration is not done
    shutil.rmtree(os.path.join(directory, "top_result/"))
    mn.main("1w7h_preparation_structure_2w.pdb", "serie_file_carbonil.conf", cpus=4, iterations=1, steps=1, temperature=1500, pele_eq_steps=1, contrl="test.conf", test=True, restart=True)
    assert glob.glob(os.path.join(directory, "top_result/*BindingEnergy*.pdb"))

    # Test restart if growing steps did not finish
    mn.main("1w7h_preparation_structure_2w.pdb", "serie_file_carbonil.conf", cpus=4, iterations=1, steps=1, temperature=1500, pele_eq_steps=1, contrl="test.conf", test=True, restart=True)
    subprocess.call("bash test_restart.sh".split())
    shutil.rmtree(os.path.join(directory, "top_result/"))
    shutil.rmtree(os.path.join(directory, "sampling_result/"))
    shutil.rmtree(os.path.join(directory, "clustering_PDBs/1/"))
    shutil.rmtree(os.path.join(directory, "growing_steps/1/"))
    mn.main("1w7h_preparation_structure_2w.pdb", "serie_file_carbonil.conf", cpus=4, iterations=1, steps=1, temperature=1500, pele_eq_steps=1, contrl="test.conf", test=True, restart=True)
    subprocess.call("bash test_restart.sh".split())
    assert glob.glob(os.path.join(directory, "top_result/*BindingEnergy*.pdb"))

    #Remove working folder
    shutil.rmtree("1w7h_preparation_structure_2w_carbonylC6-H7O1-C1")
    os.chdir("..")
