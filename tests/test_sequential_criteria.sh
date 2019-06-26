#!/bin/bash

python -m frag_pele.main -cp 1w7h_preparation_structure_2w.pdb -sef sequential.conf -x 1 --steps 1 -es 1 --temp 100000 --cpus 2 -cr currentEnergy
