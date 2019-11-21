#!/bin/bash
conda create -n frag_V2.0 -c carlesperez94 -c NostrumBioDiscovery -c conda-forge -c carlesperez94 -c rdkit frag_pele --yes
source activate frag_V2.0
cd tests/
python -m pytest test_all.py
