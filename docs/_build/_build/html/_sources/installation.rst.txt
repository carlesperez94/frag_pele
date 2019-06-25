============
Installation
============

.. toctree::
   :maxdepth: 2

Conda (Recommended)
====================

If you want to use a separate environment for frag:
    
    conda create -n frag python=3.7 --yes

    source activate frag

Else start here:

    conda install -c NostrumBioDiscovery -c conda-forge frag_pele --yes

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

Test it works (fast test)
===========================

git clone https://github.com/carlesperez94/frag_pele.git

cd frag_pele/frag_pele/example/1W7H_example/

python -m frag_pele.main -cp 1w7h_preparation_structure_2w.pdb -x 1 --steps 1 -sef --pele_eq_steps 1 sequential_frag.conf --cpus 2 --steps 1 --temp 1000000
