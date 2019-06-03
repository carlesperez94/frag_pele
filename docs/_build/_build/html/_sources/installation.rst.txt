============
Installation
============

.. toctree::
   :maxdepth: 2

Conda (Recommended)
====================

If you want to use a separate environment for frag:
    
    conda create -n frag python=3.7 (If you want to use a separate environment for frag)

    source activate frag

Else start here:

    conda install -c NostrumBioDiscovery -c conda-forge frag_pele

    change PELE schrodinger & mpirun under /conda/env/lib/pythonX/site-packages/frag_pele/constants.py

Pypi
========

Ongoing

Source code
=============

git clone https://github.com/carlesperez94/frag_pele.git

cd frag_pele

pip install numpy cython (in case not to have them)

python setup.py install

change PELE schrodinger & mpirun under /site-packages/frag_pele/constants.py
