#!/bin/bash
cd conda_recipe/
conda-build -c conda-forge -c nostrumbiodiscovery -c rdkit -c carlesperez94 --python=3.8 .
conda-build -c conda-forge -c nostrumbiodiscovery -c rdkit -c carlesperez94 --python=3.7 .
conda-build -c conda-forge -c nostrumbiodiscovery -c rdkit -c carlesperez94 --python=3.6 .
conda-build -c conda-forge -c nostrumbiodiscovery -c rdkit -c carlesperez94 --python=3.5 .

