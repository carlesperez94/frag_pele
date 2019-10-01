import pytest
import subprocess
import os



def test_HT():
    subprocess.call("bash test_HT.sh".split())
    assert os.path.exists("selected_result_aminoC1N1/epochesampling_result_aminoC1N1_trajectory_1.10_BindingEnergy-52.54600000000001.pdb") \
    or os.path.exists("selected_result_aminoC1N1/epochsampling_result_aminoC1N1_trajectory_1.10_BindingEnergy-49.3207.pdb")  \
    or os.path.exists("selected_result_aminoC1N1/epochsampling_result_aminoC1N1_trajectory_1.11_BindingEnergy-48.7485.pdb")

def test_mae():
    subprocess.call("bash test_mae.sh".split())
    assert os.path.exists("selected_result_aminoC1N1/trajectory_1.10.mae")

def test_rotamers():
    subprocess.call("bash test_rotamers.sh".split())
    
    with open("DataLocal/LigandRotamerLibs/3IP.rot.assign", "r") as f:
        lines = f.readlines()
        if "60" in lines[1]:
            assert True
        else:
            assert False

def test_criteria():
    subprocess.call("bash test_criteria.sh".split())
    assert os.path.exists("selected_result_aminoC1N1/epochsampling_result_aminoC1N1_trajectory_1.12_currentEnergy-13514.3.pdb") \
    or os.path.exists("selected_result_aminoC1N1/epochsampling_result_aminoC1N1_trajectory_1.10_currentEnergy-13516.9.pdb")

def test_chainame():
    subprocess.call("bash test_chain_core_frag.sh".split())
    assert os.path.exists("selected_result_phenyl2C6C1/epochesampling_result_phenyl2C6C1_trajectory_1.1_BindingEnergy-72.5229.pdb") \
    or os.path.exists("selected_result_phenyl2C6C1/epochesampling_result_phenyl2C6C1_trajectory_1.1_BindingEnergy-46.979.pdb") \
    or os.path.exists("selected_result_phenyl2C6C1/epochsampling_result_phenyl2C6C1_trajectory_1.1_BindingEnergy-45.3478.pdb")

def test_temperature():
    subprocess.call("bash test_temperature.sh".split())

    with open("control_template.conf", "r") as f:
        for line in f.readlines():
            if "temperature" in line:
                if "100000" in line:
                    assert True
                else:
                    assert False

def test_EX():
    subprocess.call("bash test_EX.sh".split())
    assert os.path.exists("selected_result_aminoC1N1/epochsampling_result_aminoC1N1_trajectory_1.11_BindingEnergy-48.7485.pdb")
