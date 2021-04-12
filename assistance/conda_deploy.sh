#!/bin/bash
cd conda_recipe/
conda-build -c carlesperez94 -c conda-forge -c nostrumbiodiscovery -c rdkit -c omnia -c martimunicoy --python=3.8 .
conda-build -c carlesperez94 -c conda-forge -c nostrumbiodiscovery -c rdkit -c omnia -c martimunicoy --python=3.7 .
conda-build -c carlesperez94 -c conda-forge -c nostrumbiodiscovery -c rdkit -c omnia -c martimunicoy --python=3.6 .
conda-build -c carlesperez94 -c conda-forge -c nostrumbiodiscovery -c rdkit -c omnia -c martimunicoy --python=3.5 .

