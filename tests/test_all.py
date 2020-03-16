import pytest
import glob
import subprocess
import os


def test_double():
    subprocess.call("bash test_double_bond.sh".split())
    assert glob.glob("1w7h_preparation_structure_2w_carbonylC6-H7O1-C1/top_result/*BindingEnergy*.pdb")

def test_HT():
    subprocess.call("bash test_HT.sh".split())
    assert glob.glob("1w7h_preparation_structure_2w_aminoC1N1/top_result/*BindingEnergy*.pdb")

def test_mae():
    subprocess.call("bash test_mae.sh".split())
    assert glob.glob("1w7h_preparation_structure_2w_aminoC1N1/top_result/*.mae")

def test_rotamers():
    subprocess.call("bash test_rotamers.sh".split())
    
    with open("1w7h_preparation_structure_2w_aminoC1N1/DataLocal/LigandRotamerLibs/3IP.rot.assign", "r") as f:
        lines = f.readlines()
        if "60" in lines[1]:
            assert True
        else:
            assert False

def test_criteria():
    subprocess.call("bash test_criteria.sh".split())
    assert glob.glob("1w7h_preparation_structure_2w_aminoC1N1/top_result/*currentEnergy*.pdb")

def test_chainame():
    subprocess.call("bash test_chain_core_frag.sh".split())
    assert glob.glob("1w7h_preparation_structure_2w2_phenyl2C6C1/top_result/*BindingEnergy*.pdb")

def test_temperature():
    subprocess.call("bash test_temperature.sh".split())

    with open("1w7h_preparation_structure_2w2_phenyl2C6C1/control_template.conf", "r") as f:
        for line in f.readlines():
            if "temperature" in line:
                if "100000" in line:
                    assert True
                else:
                    assert False
