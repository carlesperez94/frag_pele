import os
import frag_pele.Growing.template_fragmenter as tf

BACKBONE_ATOMS = ['_N__', '_C__', '_CA_']

def replace_atom_template(original, to_replace, atomname):
    original_indx = original.find_index_of_atom_name(atomname)
    atom_to_replace = original.list_of_atoms[original_indx]
    to_replace_indx = to_replace.find_index_of_atom_name(atomname) 
    to_replace.replace_atom(to_replace_indx, atom_to_replace)   

def replace_bond_template(original, to_replace, atomname):
    origin_atom_indx = original.find_index_of_atom_name(atomname)
    origin_bond_indx = original.find_bond_from_atom(origin_atom_indx)
    repl_atom_indx = to_replace.find_index_of_atom_name(atomname)
    repl_bond_indx = to_replace.find_bond_from_atom(repl_atom_indx)
    for o_indx in origin_bond_indx:
        bond_to_replace = original.list_of_bonds[o_indx]
        or_atom_1 = original.list_of_atoms[o_indx[0]].pdb_atom_name
        or_atom_2 = original.list_of_atoms[o_indx[1]].pdb_atom_name
        for r_indx in repl_bond_indx:
            rep_atom_1 = to_replace.list_of_atoms[r_indx[0]].pdb_atom_name
            rep_atom_2 = to_replace.list_of_atoms[r_indx[1]].pdb_atom_name
            if (or_atom_1 == rep_atom_1 and or_atom_2 == rep_atom_2) or (or_atom_1 == rep_atom_2 and or_atom_2 == rep_atom_1):
                print("Replacing bond between {} and {}. Indexes {}".format(or_atom_1, or_atom_2, r_indx))
                bond_to_replace.atom1 = r_indx[0] 
                bond_to_replace.atom2 = r_indx[1]
                to_replace.replace_bond(r_indx, bond_to_replace)


def replace_theta_template(original, to_replace, atomname):
    origin_atom_indx = original.find_index_of_atom_name(atomname)
    origin_theta_indx = original.find_theta_from_atom(origin_atom_indx)
    repl_atom_indx = to_replace.find_index_of_atom_name(atomname)
    repl_theta_indx = to_replace.find_theta_from_atom(repl_atom_indx)
    for o_indx in origin_theta_indx:
        theta_to_replace = original.list_of_thetas[o_indx]
        or_atom_1 = original.list_of_atoms[o_indx[0]].pdb_atom_name
        or_atom_2 = original.list_of_atoms[o_indx[1]].pdb_atom_name
        or_atom_3 = original.list_of_atoms[o_indx[2]].pdb_atom_name
        for r_indx in repl_theta_indx:
            rep_atom_1 = to_replace.list_of_atoms[r_indx[0]].pdb_atom_name
            rep_atom_2 = to_replace.list_of_atoms[r_indx[1]].pdb_atom_name
            rep_atom_3 = to_replace.list_of_atoms[r_indx[2]].pdb_atom_name
            if (or_atom_1 == rep_atom_1 and or_atom_2 == rep_atom_2 and or_atom_3 == rep_atom_3) or \
               (or_atom_1 == rep_atom_3 and or_atom_2 == rep_atom_2 and or_atom_3 == rep_atom_1):
                print("Replacing theta between {}, {} and {}. Indexes {}".format(or_atom_1, or_atom_2, or_atom_3, r_indx))
                theta_to_replace.atom1 = r_indx[0]
                theta_to_replace.atom2 = r_indx[1]
                theta_to_replace.atom3 = r_indx[2]
                to_replace.replace_theta(r_indx, theta_to_replace)


def correct_template(template, aminoacid_path="Data/Templates/OPLS2005/Protein/gly", work_dir="."):
    templ = tf.TemplateImpact(template)
    templ.erease_atom_from_template('_HN_')
    templ.erease_atom_from_template('_HXT')
    templ_to_copy = tf.TemplateImpact(os.path.join(work_dir, aminoacid_path))
    for name in BACKBONE_ATOMS:
        replace_atom_template(templ_to_copy, templ, name)
        replace_bond_template(templ_to_copy, templ, name)
        replace_theta_template(templ_to_copy, templ, name)
        templ.write_template_to_file(template)

def delete_atom_from_template(template, atomname):
    templ = tf.TemplateOPLS2005(template)
    if " " in atomname:
        atomname = atomname.replace(" ", "_")
    templ.erease_atom_from_template(atomname)
    templ.write_template_to_file(template)
