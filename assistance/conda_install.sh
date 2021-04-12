#!/bin/bash
conda create -n frag python=3.7
conda activate frag
conda install -c carlesperez94 -c conda-forge -c nostrumbiodiscovery -c rdkit -c omnia -c martimunicoy frag_pele
