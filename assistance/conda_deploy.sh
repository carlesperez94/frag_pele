#!/bin/bash
cd conda_recipe/
conda-build -c conda-forge -c nostrumbiodiscovery -c rdkit -c carlesperez94 .

