import sys
import os
import argparse

AMINOACIDS = ["VAL", "ASN", "GLY", "LEU", "ILE",
              "SER", "ASP", "LYS", "MET", "GLN",
              "TRP", "ARG", "ALA", "THR", "PRO",
              "PHE", "GLU", "HIS", "HIP", "TYR",
              "CYS", "HID"]

TER_CONSTR = 5

BACK_CONSTR = 0.5

CONSTR_ATOM = '''{{ "type": "constrainAtomToPosition", "springConstant": {0}, "equilibriumDistance": 0.0, "constrainThisAtom": "{1}:{2}:{3}" }},'''

CONSTR_DIST = '''{{ "type": "constrainAtomsDistance", "springConstant": {}, "equilibriumDistance": {}, "constrainThisAtom": "{}:{}:{}", "toThisOtherAtom": "{}:{}:{}" }},'''

CONSTR_CALPHA = '''{{ "type": "constrainAtomToPosition", "springConstant": {2}, "equilibriumDistance": 0.0, "constrainThisAtom": "{0}:{1}:_CA_" }},'''

class ConstraintBuilder(object):

    def __init__(self, pdb, gaps, metals):
        self.pdb = pdb
        self.gaps = gaps
        self.metals = metals

    def parse_atoms(self, interval=10):
        residues = {}
        initial_res = None
        with open(self.pdb, "r") as pdb:
            for line in pdb:
                resname = line[16:21].strip()
                atomtype = line[11:16].strip()
                resnum = line[22:26].strip()
                chain = line[20:23].strip()
                if line.startswith("ATOM") and resname in AMINOACIDS and atomtype == "CA":
                    try:
                        if not initial_res:
                            residues["initial"] = [chain, line[22:26].strip()]
                            initial_res = True
                            continue
                        # Apply constraint every 10 residues
                        elif int(resnum) % interval != 1:
                            residues["terminal"] = [chain, line[22:26].strip()]
                        elif int(resnum) % interval == 1 and line.startswith("ATOM") and resname in AMINOACIDS and atomtype == "CA":
                            residues[resnum] = chain
                    except ValueError:
                        continue
        return residues

    def build_constraint(self, residues, BACK_CONSTR=BACK_CONSTR, TER_CONSTR=TER_CONSTR):

        init_constr = ['''"constraints":[''', ]

        back_constr = [CONSTR_CALPHA.format(chain, resnum, BACK_CONSTR) for resnum, chain in residues.items() if resnum.isdigit()]

        gaps_constr = self.gaps_constraints()

        metal_constr = self.metal_constraints()

        terminal_constr = [CONSTR_CALPHA.format(residues["initial"][0], residues["initial"][1], TER_CONSTR), CONSTR_CALPHA.format(residues["terminal"][0], residues["terminal"][1], TER_CONSTR).strip(",")]

        final_constr = ["],"]

        constraints = init_constr + back_constr + gaps_constr + metal_constr + terminal_constr + final_constr

        return constraints

    def gaps_constraints(self):
        #self.gaps = {}
        gaps_constr = []
        for chain, residues in self.gaps.items():
            gaps_constr = [CONSTR_ATOM.format(TER_CONSTR, chain, terminal, "_CA_") for terminals in residues for terminal in terminals]
        return gaps_constr

    def metal_constraints(self):

        metal_constr = []
        for metal, ligands in self.metals.items():
            metal_name, chain, metnum = metal.split(" ")
            for ligand in ligands:
                ligand_info, bond_lenght = ligand
                resname, resnum, chain, ligname = ligand_info.split(" ")
                metal_constr.append(CONSTR_DIST.format(TER_CONSTR, bond_lenght, chain, resnum, ligname, chain, metnum, metal_name))
        return metal_constr


def retrieve_constraints(pdb_file, gaps, metal, back_constr=BACK_CONSTR, ter_constr=TER_CONSTR, interval=10):
    constr = ConstraintBuilder(pdb_file, gaps, metal)
    residues = constr.parse_atoms(interval=interval)
    constraints = constr.build_constraint(residues, back_constr, ter_constr)
    return constraints

def parseargs():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('pdb', type=str, help='pdb to create the contraints on')
    parser.add_argument('conf', help='Control file to fill in. It need to templetazide with $CONSTRAINTS')
    parser.add_argument('--interval', type=int, help="Every how many CA to constraint")
    parser.add_argument('--ca', type=float, help="Constraint value to use on backbone CA", default=BACK_CONSTR)
    parser.add_argument('--terminal', type=float, help="Constraint value to use on terminal CA", default=TER_CONSTR)
    args = parser.parse_args()
    return os.path.abspath(args.pdb), os.path.abspath(args.conf), args.interval, args.conf, args.ca, args.terminal

if __name__ == "__main__":
    pdb, conf, interval, conf, back_constr, ter_constr = parseargs()
    constraints = retrieve_constraints(pdb, {}, {}, back_constr, ter_constr, interval)
