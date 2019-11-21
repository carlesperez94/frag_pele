#!/bin/bash
conda install -c NostrumBioDiscovery -c conda-forge -c carlesperez94 -c rdkit frag_pele
conda install -c conda-forge pytest
cd tests/
python -m pytest test_all.py
