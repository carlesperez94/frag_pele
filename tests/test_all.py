import pytest
import shutil
import glob
import subprocess
import os


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
    subprocess.call("bash test_double_bond.sh".split())
    assert glob.glob("1w7h_preparation_structure_2w_carbonylC6-H7O1-C1/top_result/*BindingEnergy*.pdb")

def test_HT():
    subprocess.call("bash test_HT.sh".split())
    assert glob.glob("1w7h_preparation_structure_2w_aminoC1N1/top_result/*traject*.mae")

def test_sequential():
    subprocess.call("bash test_sequential.sh".split())
    assert glob.glob("1w7h_preparation_structure_2w_aminoC1N1/top_result/*BindingEnergy*.pdb") and glob.glob("top_result_aminoC1N1phenylC6C1/top_result/*BindingEnergy*.pdb")

def test_flags():
    subprocess.call("bash test_flags.sh".split())
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
    # Test restart if simulation is don
    subprocess.call("bash test_restart.sh".split())
    assert glob.glob(os.path.join(directory, "top_result/*BindingEnergy*.pdb"))
    # Test restart if equilibration is not done
    shutil.rmtree(os.path.join(directory, "top_result/"))
    subprocess.call("bash test_restart.sh".split())
    assert glob.glob(os.path.join(directory, "top_result/*BindingEnergy*.pdb"))
    # Test restart if growing steps did not finish
    subprocess.call("bash test_restart.sh".split())
    shutil.rmtree(os.path.join(directory, "top_result/"))
    shutil.rmtree(os.path.join(directory, "sampling_result/"))
    shutil.rmtree(os.path.join(directory, "clustering_PDBs/1/"))
    shutil.rmtree(os.path.join(directory, "growing_steps/1/"))
    subprocess.call("bash test_restart.sh".split())
    assert glob.glob(os.path.join(directory, "top_result/*BindingEnergy*.pdb"))
    #Remove working folder
    shutil.rmtree("1w7h_preparation_structure_2w_carbonylC6-H7O1-C1")
    os.chdir("..")
