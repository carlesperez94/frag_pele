import frag_pele.Growing.template_fragmenter as tf

BACKBONE_ATOMS = ['_N__', '_C__', '_CA_']

def correct_template(template):
    templ = tf.TemplateOPLS2005(template)
    templ.erease_atom_from_template('_HN_')
    templ.erease_atom_from_template('_HXT')
    templ.write_template_to_file("tst")


correct_template("/gpfs/projects/bsc72/FragPELE/FragPELE2.3.0_testing/frag_pele/tests/receptor_145_cys_frag_resSGC4/DataLocal/Templates/OPLS2005/Protein/templates_generated/grwz")
