In "run_example" you will find an example of command to run the program. Things that you need to run it properly:

command: python3 ../../grow_for_pele.py -i mbez -f pyjz -oa "_H8_" -fa "_C8_" -c control_template.conf -p 181l_prepared_met.pdb

- Initial ligand template ("DataLocal/Templates/OPLS2005/Heteroatoms/growing_templates/mbez"). 
- Final ligand template ("DataLocal/Templates/OPLS2005/Heteroatoms/growing_templates/pyjz").
- Control file template ("control_template.conf"). Control file templatized. The results path will change along different steps of the simulation.
- PDB file with the dummy-like atoms("181l_prepared_met.pdb"). Ensure that the ligand name is the same that you have in the final ligand template.
- Original atom ("_H8_"). Atom PDB name of the initial ligand template whose parameters we want to use to transform another atom of the final template.
- Final atom ("_C8_"). Atom PDB name of the final ligand template that we want to transform into another of the initial template. 


- 
