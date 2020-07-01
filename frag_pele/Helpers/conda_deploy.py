import os

PYTHONS = ["3.8", "3.7", "3.6"]

for python in PYTHONS:
    print("conda build -c conda-forge -c rdkit -c nostrumbiodiscovery conda_recipe/ --python={}".format(python))
    os.system("conda build -c conda-forge -c rdkit -c nostrumbiodiscovery conda_recipe/ --python={}".format(python))
