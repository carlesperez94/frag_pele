import pytest
import subprocess
import os



def test_sequential():
    subprocess.call("bash test_sequential.sh".split())
    assert os.path.exists("selected_result_aminoC1N1phenylN1-H2C1-H1aminoC6N1/epochequilibration_result_aminoC1N1phenylN1-H2C1-H1aminoC6N1_trajectory_1.1_BindingEnergy-44.1931.pdb")

def test_mae():
    subprocess.call("bash test_sequential_mae.sh".split())
    assert os.path.exists("selected_result_aminoC1N1phenylN1-H2C1-H1aminoC6N1/trajectory_1.1.mae")

def test_criteria():
    subprocess.call("bash test_sequential_criteria.sh".split())
    assert os.path.exists("selected_result_aminoC1N1phenylN1-H2C1-H1aminoC6N1/epochequilibration_result_aminoC1N1phenylN1-H2C1-H1aminoC6N1_trajectory_1.1_currentEnergy-13516.2.pdb")

#def test_chainame():
#    subprocess.call("bash test_chain_core_frag.sh".split())
#    assert os.path.exists("selected_result_phenyl2C6C1/epochequilibration_result_phenyl2C6C1_trajectory_1.1_BindingEnergy-72.5229.pdb") \
#    or os.path.exists("selected_result_phenyl2C6C1/epochequilibration_result_phenyl2C6C1_trajectory_1.1_BindingEnergy-46.979.pdb")
