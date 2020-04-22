#!/bin/bash

python -m frag_pele.main -cp 1w7h_preparation_structure_2w2.pdb -sef serie_file2.conf -cc Z -fc Z --cpus 4 -x 1 --steps 1 -es 1  --temp 100001 -cr currentEnergy -rot 60 -r "bests" -sd 3 -st 3 -trh 33 -roth 33 -trl 31 -rotl 31 -rad 3 -pdbf "best.pdb" --debug -miov 3 -maov 3
