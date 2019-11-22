#!/bin/bash
conda create -n frag python=3.7
source activate frag
conda install -c NostrumBioDiscovery -c conda-forge -c rdkit frag_pele
cd tests/
python -m pytest test_all.py
