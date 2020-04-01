import logging

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

# Definition of reggex patterns
HEADER_OPLS2005 = "* LIGAND DATABASE FILE (OPLS2005)\n*\n"
PATTERN_OPLS2005_RESX_HEADER = "{:5} {:6d} {:6d} {:7d} {:7d} {:8d} \n"
PATTERN_OPLS2005_RESX_LINE = "{:5d} {:5d} {:1}   {:4} {:4} {:5d} {: >11.5f} {: >11.5f} {: >11.5f}\n"
PATTERN_OPLS2005_NBON = "{:5d} {: >8.4f} {: >8.4f} {: >10.6f} {: >8.4f} {: >8.4f} {: >13.9f} {: >13.9f}\n"
PATTERN_OPLS2005_BOND = "{:5d} {:5d} {:>9.3f} {:>6.3f}\n"
PATTERN_OPLS2005_THETA = "{:5d} {:5d} {:5d} {:>11.5f} {: >11.5f}\n"
PATTERN_OPLS2005_PHI = "{:5d} {:5d} {: 5d} {:5d} {:>9.5f} {: >4.1f} {: >3.1f}\n"


class Atom:
    """A class which contains all the information and properties of an Atom, and several methods to build templates from
    this data (currently only in OPLS2005)."""
    def __init__(self, atom_id, parent_id, location, atom_type, pdb_atom_name, unknown, x_zmatrix=0, y_zmatrix=0,
                 z_zmatrix=0, sigma=0, epsilon=0, charge=0, radnpSGB=0, radnpType=0, sgbnpGamma=0, sgbnpType=0,
                 is_linker=False, is_fragment=False):
        """
        :param atom_id: ID of the atom in the template.
        :type atom_id: int
        :param parent_id: ID ot the parent atom.
        :type parent_id: int
        :param location: Location of the atom. M for backbone, S for side chain; ligand links and non-ligand links use
        these classification for some checks. However, protein links define its backbone atoms independently of this
        flag; protein residue templates have M only for N, C and CA, under whatever name is actually used for those
        atoms.
        :type location: str
        :param atom_type: Atom type. For example: CA (carbon aromatic), CT (carbon terminal), etc.
        :type atom_type: str
        :param pdb_atom_name: PDB atom name.
        :type pdb_atom_name: str
        :param unknown: Nobody knows what is this...
        :type unknown: int
        :param x_zmatrix: Coord X of the Z-matrix.
        :type x_zmatrix: float
        :param y_zmatrix: Coord Y of the Z-matrix.
        :type y_zmatrix: float
        :param z_zmatrix: Coord Z of the Z-matrix.
        :type z_zmatrix: float
        :param sigma: sigma value, used to compute Van Der Waals terms. Units in Armstrong.
        :type sigma: float
        :param epsilon: epsilon value, used to compute Van Der Waals terms. Units in kcal/mol.
        :type epsilon: float
        :param charge: charge value, used to compute electrostatic potentials. Units are elementary charge.
        :type charge: float
        :param radnpSGB: radii of non polar SGB. Atomic radii used to calculate the surface of the molecules when
        obtaining SGB Born radii.
        :type radnpSGB: float
        :param radnpType: radii of non polar Type. Atomic radii used to calculate SASA in the non-polar term of the SGB
        and VDGBNP models
        :type radnpType: float
        :param sgbnpGamma: SGB non polar Gamma. Gamma parameter of the nonpolar model.
        :type sgbnpGamma: float
        :param sgbnpType: SGB non polar Type. Alpha parameter for the nonpolar model
        :param is_linker: Flag set when the atom is linking the fragment and the core.
        :type is_linker: bool
        :param is_fragment: Flag set when the atom is of the fragment.
        :type is_fragment: bool
        """
        self.atom_id = int(atom_id)
        self.parent_id = int(parent_id)
        self.location = str(location)
        self.atom_type = str(atom_type)
        self.pdb_atom_name = str(pdb_atom_name)
        self.unknown = int(unknown)
        self.x_zmatrix = float(x_zmatrix)
        self.y_zmatrix = float(y_zmatrix)
        self.z_zmatrix = float(z_zmatrix)
        self.sigma = float(sigma)
        self.epsilon = float(epsilon)
        self.charge = float(charge)
        self.radnpSGB = float(radnpSGB)
        self.radnpType = float(radnpType)
        self.sgbnpGamma = float(sgbnpGamma)
        self.sgbnpType = float(sgbnpType)
        self.bonds = []
        self.thetas = []
        self.phis = []
        self.iphis = []
        self.is_fragment = bool(is_fragment)
        self.is_linker = bool(is_linker)

    def write_resx(self):
        return PATTERN_OPLS2005_RESX_LINE.format(self.atom_id, self.parent_id, self.location, self.atom_type.strip(),
                                                 self.pdb_atom_name, self.unknown, self.x_zmatrix,
                                                 self.y_zmatrix, self.z_zmatrix)

    def write_nbon(self):
        return PATTERN_OPLS2005_NBON.format(self.atom_id, self.sigma, self.epsilon, self.charge, self.radnpSGB,
                                            self.radnpType, self.sgbnpGamma, self.sgbnpType)


class Bond:
    def __init__(self, atom1, atom2, spring, eq_dist, is_fragment=False, is_linker=False):
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.spring = float(spring)
        self.eq_dist = float(eq_dist)
        self.is_fragment = bool(is_fragment)
        self.is_linker = bool(is_linker)

    def write_bond(self):
        return PATTERN_OPLS2005_BOND.format(self.atom1, self.atom2, self.spring, self.eq_dist)


class Theta:
    def __init__(self, atom1, atom2, atom3, spring, eq_angle, is_fragment=False):
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.atom3 = int(atom3)
        self.spring = float(spring)
        self.eq_angle = float(eq_angle)
        self.is_fragment = bool(is_fragment)

    def write_theta(self):
        return PATTERN_OPLS2005_THETA.format(self.atom1, self.atom2, self.atom3, self.spring, self.eq_angle)


class Phi:
    def __init__(self, atom1, atom2, atom3, atom4, constant, prefactor, nterm, improper, is_fragment=False):
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.atom3 = int(atom3)
        self.atom4 = int(atom4)
        self.constant = float(constant)
        self.prefactor = float(prefactor)
        self.nterm = float(nterm)
        self.improper = bool(improper)
        self.is_fragment = bool(is_fragment)

    def write_phi(self):
        if not self.improper:
            return PATTERN_OPLS2005_PHI.format(self.atom1, self.atom2, self.atom3, self.atom4, self.constant,
                                               self.prefactor, self.nterm)

    def write_iphi(self):
        if self.improper:
            return PATTERN_OPLS2005_PHI.format(self.atom1, self.atom2, self.atom3, self.atom4, self.constant,
                                               self.prefactor, self.nterm)


class TemplateOPLS2005:
    def __init__(self, path_to_template):
        self.path_to_template = path_to_template
        self.template_name = ""
        self.num_nbon_params = 0
        self.num_bond_params = 0
        self.num_angle_params = 0
        self.num_dihedr_params = 0
        self.num_nonnull = 0
        self.list_of_atoms = {}
        self.list_of_bonds = {}
        self.list_of_thetas = {}
        self.list_of_phis = []
        self.list_of_iphis = []
        self.unique_atoms = []
        self.read_template()

    def read_template(self):
        template = file_to_list_of_lines(self.path_to_template)
        for line in template[2:3]:
            self.template_name = get_string_from_line(line=line, index_initial=0, index_final=5)
            self.num_nbon_params = int(get_string_from_line(line=line, index_initial=6, index_final=11))
            self.num_bond_params = int(get_string_from_line(line=line, index_initial=13, index_final=17))
            self.num_angle_params = int(get_string_from_line(line=line, index_initial=18, index_final=24))
            self.num_dihedr_params = int(get_string_from_line(line=line, index_initial=25, index_final=31))
            self.num_nonnull = int(get_string_from_line(line=line, index_initial=32, index_final=39))
        for line in template[3:]:
            if line.startswith("NBON"):
                index = template.index(line)
                break
            try:
                atom_id = get_string_from_line(line=line, index_initial=0, index_final=6)
                parent_id = get_string_from_line(line=line, index_initial=6, index_final=11)
                location = get_string_from_line(line=line, index_initial=12, index_final=13)
                atom_type = get_string_from_line(line=line, index_initial=15, index_final=20)
                pdb_atom_name = get_string_from_line(line=line, index_initial=21, index_final=25)
                unknown = get_string_from_line(line=line, index_initial=26, index_final=31)
                x_zmatrix = get_string_from_line(line=line, index_initial=32, index_final=43)
                y_zmatrix = get_string_from_line(line=line, index_initial=44, index_final=55)
                z_zmatrix = get_string_from_line(line=line, index_initial=56, index_final=67)
                atom = Atom(atom_id=atom_id, parent_id=parent_id, location=location, atom_type=atom_type,
                            pdb_atom_name=pdb_atom_name, unknown=unknown, x_zmatrix=x_zmatrix, y_zmatrix=y_zmatrix,
                            z_zmatrix=z_zmatrix)
                self.list_of_atoms.setdefault(atom.atom_id, atom)
                if pdb_atom_name not in self.unique_atoms:
                    self.unique_atoms.append(pdb_atom_name)
                else:
                    raise ValueError("ERROR: PDB ATOM NAME {} ALREADY EXISTS in the template {}!".format(pdb_atom_name,
                                                                                              self.path_to_template))
            except ValueError:
                raise ValueError(
                     "Unexpected type in line {} of {}\n{}".format(template.index(line), self.path_to_template, line))

        for line in template[index + 1:]:
            if line.startswith("BOND"):
                index = template.index(line)
                break
            try:
                id = int(get_string_from_line(line=line, index_initial=0, index_final=6))
                self.list_of_atoms[id].sigma = float(get_string_from_line(line=line, index_initial=7, index_final=14))
                self.list_of_atoms[id].epsilon = float(get_string_from_line(line=line, index_initial=15, index_final=23))
                self.list_of_atoms[id].charge = float(get_string_from_line(line=line, index_initial=24, index_final=34))
                self.list_of_atoms[id].radnpSGB = float(get_string_from_line(line=line, index_initial=35, index_final=43))
                self.list_of_atoms[id].radnpType = float(get_string_from_line(line=line, index_initial=44, index_final=52))
                self.list_of_atoms[id].sgbnpGamma = float(get_string_from_line(line=line, index_initial=53,
                                                                               index_final=66))
                self.list_of_atoms[id].sgbnpType = float(get_string_from_line(line=line, index_initial=67, index_final=80))
            except ValueError:
                raise ValueError(
                    "Unexpected type in line {} of {}\n{}".format(template.index(line), self.path_to_template, line))

        for line in template[index + 1:]:
            if line.startswith("THET"):
                index = template.index(line)
                break
            try:
                id_atom1 = int(get_string_from_line(line=line, index_initial=0, index_final=6))
                id_atom2 = int(get_string_from_line(line=line, index_initial=6, index_final=12))
                spring = get_string_from_line(line=line, index_initial=13, index_final=21)
                eq_dist = get_string_from_line(line=line, index_initial=23, index_final=29)
                # Create bond instance
                bond = Bond(atom1=id_atom1, atom2=id_atom2, spring=spring, eq_dist=eq_dist)
                self.list_of_bonds.setdefault((id_atom1, id_atom2), bond)
                # Set which atom is bonded with
                self.list_of_atoms[id_atom1].bonds.append(bond)

            except ValueError:
                raise ValueError(
                    "Unexpected type in line {} of {}\n{}".format(template.index(line), self.path_to_template, line))
        for line in template[index + 1:]:
            if line.startswith("PHI"):
                index = template.index(line)
                break
            try:
                id_atom1 = int(get_string_from_line(line=line, index_initial=0, index_final=6))
                id_atom2 = int(get_string_from_line(line=line, index_initial=6, index_final=12))
                id_atom3 = int(get_string_from_line(line=line, index_initial=13, index_final=18))
                spring = get_string_from_line(line=line, index_initial=19, index_final=29)
                eq_angle = get_string_from_line(line=line, index_initial=31, index_final=40)
                # Create bond instance
                theta = Theta(atom1=id_atom1, atom2=id_atom2, atom3=id_atom3, spring=spring, eq_angle=eq_angle)
                self.list_of_thetas.setdefault((id_atom1, id_atom2, id_atom3), theta)
                self.list_of_atoms[id_atom1].thetas.append(theta)

            except ValueError:
                raise ValueError(
                    "Unexpected type in line {} of {}\n{}".format(template.index(line), self.path_to_template,
                                                                  line))
        for line in template[index + 1:]:
            if line.startswith("IPHI"):
                index = template.index(line)
                break
            try:
                id_atom1 = int(get_string_from_line(line=line, index_initial=0, index_final=5))
                id_atom2 = int(get_string_from_line(line=line, index_initial=6, index_final=11))
                id_atom3 = int(get_string_from_line(line=line, index_initial=12, index_final=17))
                id_atom4 = int(get_string_from_line(line=line, index_initial=18, index_final=23))
                constant = get_string_from_line(line=line, index_initial=25, index_final=32)
                preafactor = get_string_from_line(line=line, index_initial=33, index_final=38)
                nterm = get_string_from_line(line=line, index_initial=39, index_final=42)
                # Create bond instance
                phi = Phi(atom1=id_atom1, atom2=id_atom2, atom3=id_atom3, atom4=id_atom4, constant=constant,
                          prefactor=preafactor, nterm=nterm, improper=False)
                self.list_of_phis.append(phi)
                self.list_of_atoms[id_atom1].phis.append(phi)

            except ValueError:
                raise ValueError(
                    "Unexpected type in line {} of {}\n{}".format(template.index(line), self.path_to_template,
                                                                  line))
        for line in template[index + 1:]:
            if line.startswith("END"):
                break
            try:
                id_atom1 = int(get_string_from_line(line=line, index_initial=0, index_final=6))
                id_atom2 = int(get_string_from_line(line=line, index_initial=7, index_final=12))
                id_atom3 = int(get_string_from_line(line=line, index_initial=13, index_final=18))
                id_atom4 = int(get_string_from_line(line=line, index_initial=19, index_final=24))
                constant = get_string_from_line(line=line, index_initial=26, index_final=34)
                preafactor = get_string_from_line(line=line, index_initial=34, index_final=39)
                nterm = get_string_from_line(line=line, index_initial=40, index_final=43)
                # Create bond instance
                phi = Phi(atom1=id_atom1, atom2=id_atom2, atom3=id_atom3, atom4=id_atom4, constant=constant,
                          prefactor=preafactor, nterm=nterm, improper=True)
                self.list_of_iphis.append(phi)

            except ValueError:
                raise ValueError(
                    "Unexpected type in line {} of {}\n{}".format(template.index(line), self.path_to_template,
                                                                  line))

    def write_header(self):
        return HEADER_OPLS2005+PATTERN_OPLS2005_RESX_HEADER.format(self.template_name, self.num_nbon_params,
                                                                   self.num_bond_params, self.num_angle_params,
                                                                   self.num_dihedr_params, self.num_nonnull)

    def write_xres(self):
        content = []
        for n in range(1, len(self.list_of_atoms)+1):
            line = self.list_of_atoms[n].write_resx()
            content.append(line)
        return "".join(content)

    def write_nbon(self):
        content = []
        for n in range(1, len(self.list_of_atoms)+1):
            line = self.list_of_atoms[n].write_nbon()
            content.append(line)
        return "".join(content)

    def write_bond(self):
        content = []
        for key in self.list_of_bonds.keys():
            line = self.list_of_bonds[key].write_bond()
            content.append(line)
        return "".join(content)

    def write_theta(self):
        content = []
        for key in self.list_of_thetas.keys():
            line = self.list_of_thetas[key].write_theta()
            content.append(line)
        return "".join(content)

    def write_phis(self):
        content = []
        for phi in self.list_of_phis:
            line = phi.write_phi()
            content.append(line)
        return "".join(content)

    def write_iphis(self):
        content = []
        for phi in self.list_of_iphis:
            line = phi.write_iphi()
            content.append(line)
        return "".join(content)

    def write_template(self):
        header = self.write_header()
        xres = self.write_xres()
        nbon_header = "NBON\n"
        nbon = self.write_nbon()
        bond_header = "BOND\n"
        bond = self.write_bond()
        theta_header = "THET\n"
        theta = self.write_theta()
        phis_header = "PHI\n"
        phis = self.write_phis()
        iphis_header = "IPHI\n"
        iphis = self.write_iphis()
        ending = "END"

        return header+xres+nbon_header+nbon+bond_header+bond+theta_header+theta+phis_header+phis+iphis_header+iphis+ending

    def write_template_to_file(self, template_new_name=None):
        if not template_new_name:
            name = self.template_name.lower()+"z"
        else:
            name = template_new_name
        with open(name, "w") as template:
            template.write(self.write_template())

    def get_list_of_fragment_atoms(self):
        atoms = []
        for key, atom in self.list_of_atoms.items():
            if atom.is_fragment:
                atoms.append((key, atom))
        return atoms
    
    def get_list_of_core_atoms(self):
        atoms = []
        for key, atom in self.list_of_atoms.items():
            if not atom.is_fragment:
                atoms.append((key, atom))
        return atoms
 
    def get_list_of_fragment_bonds(self):
        bonds = []
        for key, bond in self.list_of_bonds.items():
            if bond.is_fragment:
                bonds.append((key, bond))
        return bonds

    def get_list_of_fragment_thetas(self):
        thetas = []
        for key, theta in self.list_of_thetas.items():
            if theta.is_fragment:
                thetas.append((key, theta))
        return thetas

    def get_list_of_fragment_phis(self):
        phis = []
        for phi in self.list_of_phis:
            if phi.is_fragment:
                phis.append(phi)
        return phis

    def get_list_of_fragment_iphis(self):
        iphis = []
        for iphi in self.list_of_iphis:
            if iphi.is_fragment:
                iphis.append(iphi)
        return iphis


class ReduceProperty:
    def __init__(self, template, lambda_to_reduce, template_core=None):
        self.template = template
        self.lambda_to_reduce = lambda_to_reduce
        self.template_core = template_core

    def reduce_epsilons(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            epsilon = atom.epsilon
            result = function(epsilon)
            self.template.list_of_atoms[key].epsilon = result

    def reduce_sigmas(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            sigma = atom.sigma
            result = function(sigma)
            self.template.list_of_atoms[key].sigma = result

    def reduce_charges(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            charge = atom.charge
            result = function(charge)
            self.template.list_of_atoms[key].charge = result
    
    def reduce_sgbnpGamma(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            sgbnpGamma = atom.sgbnpGamma
            result = function(sgbnpGamma)
            self.template.list_of_atoms[key].sgbnpGamma = result

    def reduce_sgbnpType(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            sgbnpType = atom.sgbnpType
            result = function(sgbnpType)
            self.template.list_of_atoms[key].sgbnpType = result

    def reduce_radnpSGB(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            radnpSGB = atom.radnpSGB
            result = function(radnpSGB)
            self.template.list_of_atoms[key].radnpSGB = result

    def reduce_radnpType(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            radnpType = atom.radnpType
            result = function(radnpType)
            self.template.list_of_atoms[key].radnpType = result

    def reduce_bond_eq_dist(self, function):
        bonds = self.template.get_list_of_fragment_bonds()
        for key, bond in bonds:
            eq_dist = bond.eq_dist
            result = function(eq_dist)
            self.template.list_of_bonds[key].eq_dist = result

    def modify_core_epsilons(self, function):
        atoms_grw = self.template.get_list_of_core_atoms()
        atoms_ini = self.template_core.get_list_of_core_atoms()
        for atom_grow in atoms_grw:
            index_g, atom_g = atom_grow
            for atom_core in atoms_ini:
               index_c, atom_c = atom_core
               if atom_g.pdb_atom_name == atom_c.pdb_atom_name:
                   result = function(atom_c.epsilon, atom_g.epsilon)
                   self.template.list_of_atoms[index_g].epsilon = result

    def modify_core_sigmas(self, function):
        atoms_grw = self.template.get_list_of_core_atoms()
        atoms_ini = self.template_core.get_list_of_core_atoms()
        for atom_grow in atoms_grw:
            index_g, atom_g = atom_grow
            for atom_core in atoms_ini:
               index_c, atom_c = atom_core
               if atom_g.pdb_atom_name == atom_c.pdb_atom_name:
                   result = function(atom_c.sigma, atom_g.sigma)
                   self.template.list_of_atoms[index_g].sigma = result
    
    def modify_core_charges(self, function):
        atoms_grw = self.template.get_list_of_core_atoms()
        atoms_ini = self.template_core.get_list_of_core_atoms()
        for atom_grow in atoms_grw:
            index_g, atom_g = atom_grow
            for atom_core in atoms_ini:
               index_c, atom_c = atom_core
               if atom_g.pdb_atom_name == atom_c.pdb_atom_name:
                   result = function(atom_c.charge, atom_g.charge)
                   self.template.list_of_atoms[index_g].charge = result

    def modify_core_bond_eq_dist(self, function):
        bonds_grw = self.template
        bonds_ini = self.template_core
        for bond_grw in bonds_grw.list_of_bonds.items():
            key_grw, bond_g = bond_grw
            # Get PDB atom names of the bonding atoms
            atoms_g = [bonds_grw.list_of_atoms[bond_g.atom1].pdb_atom_name, 
                       bonds_grw.list_of_atoms[bond_g.atom2].pdb_atom_name]
            for bond_ini in bonds_ini.list_of_bonds.items():
                key_ini, bond_i = bond_ini
                atoms_i = [bonds_ini.list_of_atoms[bond_i.atom1].pdb_atom_name, 
                       bonds_ini.list_of_atoms[bond_i.atom2].pdb_atom_name]
                if sorted(atoms_i) == sorted(atoms_g) and not bond_g.is_linker:
                    result = function(bond_i.eq_dist, bond_g.eq_dist)
                    bond_g.eq_dist = result
                
                
                


    def modify_core_sgbnpGamma(self, function):
        atoms_grw = self.template.get_list_of_core_atoms()
        atoms_ini = self.template_core.get_list_of_core_atoms()
        for atom_grow in atoms_grw:
            index_g, atom_g = atom_grow
            for atom_core in atoms_ini:
               index_c, atom_c = atom_core
               if atom_g.pdb_atom_name == atom_c.pdb_atom_name:
                   result = function(atom_c.sgbnpGamma, atom_g.sgbnpGamma)
                   self.template.list_of_atoms[index_g].sgbnpGamma = result   

    def modify_core_sgbnpType(self, function):
        atoms_grw = self.template.get_list_of_core_atoms()
        atoms_ini = self.template_core.get_list_of_core_atoms()
        for atom_grow in atoms_grw:
            index_g, atom_g = atom_grow
            for atom_core in atoms_ini:
               index_c, atom_c = atom_core
               if atom_g.pdb_atom_name == atom_c.pdb_atom_name:
                   result = function(atom_c.sgbnpType, atom_g.sgbnpType)
                   self.template.list_of_atoms[index_g].sgbnpType = result

    def modify_core_radnpSGB(self, function):
        atoms_grw = self.template.get_list_of_core_atoms()
        atoms_ini = self.template_core.get_list_of_core_atoms()
        for atom_grow in atoms_grw:
            index_g, atom_g = atom_grow
            for atom_core in atoms_ini:
               index_c, atom_c = atom_core
               if atom_g.pdb_atom_name == atom_c.pdb_atom_name:
                   result = function(atom_c.radnpSGB, atom_g.radnpSGB)
                   self.template.list_of_atoms[index_g].radnpSGB = result

    def modify_core_radnpType(self, function):
        atoms_grw = self.template.get_list_of_core_atoms()
        atoms_ini = self.template_core.get_list_of_core_atoms()
        for atom_grow in atoms_grw:
            index_g, atom_g = atom_grow
            for atom_core in atoms_ini:
               index_c, atom_c = atom_core
               if atom_g.pdb_atom_name == atom_c.pdb_atom_name:
                   result = function(atom_c.radnpType, atom_g.radnpType)
                   self.template.list_of_atoms[index_g].radnpType = result


class ReduceLinearly(ReduceProperty):
    def __init__(self, template, lambda_to_reduce, template_core=None):
        ReduceProperty.__init__(self, template, lambda_to_reduce, template_core)

    def reduce_value(self, value):
        result = value * self.lambda_to_reduce
        return result
    
    def reduce_value_from_diference(self, value_init, value_grown):
        diff = value_grown - value_init
        result = (diff * self.lambda_to_reduce) + value_init
        return result
     

class ReduceExponentially(ReduceProperty):
    def __init__(self, template, lambda_to_reduce):
        ReduceProperty.__init__(self, template, lambda_to_reduce)

    def reduce_value(self, value):
        result = value * self.lambda_to_reduce ** 2
        return result


def file_to_list_of_lines(file_path):
    with open(file_path, "r") as template:
        content = template.readlines()
    return content


def get_string_from_line(line, index_initial, index_final):
    string = line[index_initial:index_final]
    return string.strip()


def find_equal_pdb_atom_names(template1, template2):
    pdb_atom_names_tmpl_1 = [template1.list_of_atoms[n].pdb_atom_name for n in range(1, len(template1.list_of_atoms)+1)]
    pdb_atom_names_tmpl_2 = [template2.list_of_atoms[n].pdb_atom_name for n in range(1, len(template2.list_of_atoms)+1)]
    return list(set(pdb_atom_names_tmpl_1).intersection(pdb_atom_names_tmpl_2))


def detect_atoms(template_initial, template_grown, hydrogen_to_replace):
    fragment_atoms = []
    core_atoms_in = []
    core_atoms_grown = []
    core_atoms_names = find_equal_pdb_atom_names(template_initial, template_grown)
    # First fill the list of the template grown
    for key, atom in template_grown.list_of_atoms.items():
        pdb_atom_name = atom.pdb_atom_name
        if pdb_atom_name not in core_atoms_names:
            fragment_atoms.append(atom)
        elif hydrogen_to_replace in pdb_atom_name:
            fragment_atoms.append(atom)
        elif pdb_atom_name in core_atoms_names:
            core_atoms_grown.append(atom)
    # Then the ones of the initial template
    for key, atom in template_initial.list_of_atoms.items():
        pdb_atom_name = atom.pdb_atom_name
        if pdb_atom_name in core_atoms_names:
            core_atoms_in.append(atom)
    
    return fragment_atoms, core_atoms_in, core_atoms_grown


def set_fragment_atoms(list_of_fragment_atoms):
    for atom in list_of_fragment_atoms:
        atom.is_fragment = True


def detect_fragment_bonds(list_of_fragment_atoms, template_grown):
    fragment_bonds = []
    fragment_indexes = []
    for atom in list_of_fragment_atoms:
        fragment_indexes.append(atom.atom_id)
    fragment_indexes = list(set(fragment_indexes))
    for key, bond in template_grown.list_of_bonds.items():
        if key[0] in fragment_indexes and key[1] in fragment_indexes:
            fragment_bonds.append(bond)
    return fragment_bonds


def set_linker_bond(template):
    for key, bond in template.list_of_bonds.items():
        if template.list_of_atoms[bond.atom1].is_linker or template.list_of_atoms[bond.atom2].is_linker: 
            #Only the core one will be detected, thus we must check if the bond has fragment atoms to select the linker
            if template.list_of_atoms[bond.atom1].is_fragment or template.list_of_atoms[bond.atom2].is_fragment: 
                bond.is_linker=True


def set_fragment_bonds(list_of_fragment_bonds):
    for bond in list_of_fragment_bonds:
        bond.is_fragment = True


def set_connecting_atom(template_grown, pdb_atom_name):
    for key, atom in template_grown.list_of_atoms.items():
        if pdb_atom_name in atom.pdb_atom_name:
            atom.is_linker = True


def reduce_fragment_parameters_linearly(template_object, lambda_to_reduce):
    reductor = ReduceLinearly(template_object, lambda_to_reduce)
    reductor.reduce_sigmas(reductor.reduce_value)
    reductor.reduce_epsilons(reductor.reduce_value)
    reductor.reduce_charges(reductor.reduce_value)
    reductor.reduce_bond_eq_dist(reductor.reduce_value)
    reductor.reduce_radnpSGB(reductor.reduce_value)
    reductor.reduce_radnpType(reductor.reduce_value)
    reductor.reduce_sgbnpGamma(reductor.reduce_value)
    reductor.reduce_sgbnpType(reductor.reduce_value)


def modify_core_parameters_linearly(template_grow, lambda_to_reduce, template_core):
    reductor = ReduceLinearly(template_grow, lambda_to_reduce, template_core)
    reductor.modify_core_epsilons(reductor.reduce_value_from_diference)
    reductor.modify_core_sigmas(reductor.reduce_value_from_diference)
    reductor.modify_core_charges(reductor.reduce_value_from_diference)
    reductor.modify_core_bond_eq_dist(reductor.reduce_value_from_diference)
    reductor.modify_core_radnpSGB(reductor.reduce_value_from_diference)
    reductor.modify_core_radnpType(reductor.reduce_value_from_diference)
    reductor.modify_core_sgbnpGamma(reductor.reduce_value_from_diference)
    reductor.modify_core_sgbnpType(reductor.reduce_value_from_diference)


def main(template_initial_path, template_grown_path, step, total_steps, hydrogen_to_replace, core_atom_linker,
         tmpl_out_path):
    """
    Module to modify templates, currently working in OPLS2005. This main function basically compares two templates;
    an initial and a grown one, extracting the atoms of the fragment (that have been grown). Then, it uses this data
    to modify Linearly different attributes of the template, particularly, sigmas, charges, bond equilibrium distance,
    and the radius non polar SGB from atoms and bonds of the fragment. This modification is performed according to a
    lambda parameter that is computed dividing the current step by the total number of steps. Finally, the template is
    modified and written again to an output file.
    :param template_initial_path: Path to an OPLS2005 template of the core ligand.
    :type template_initial_path: str
    :param template_grown_path: Path to an OPLS2005 template of the ligand with the fragment added to the core.
    :type template_grown_path: str
    :param step: Current step of the total steps.
    :type step: int
    :param total_steps: Total number of steps.
    :type total_steps: int
    :param hydrogen_to_replace: PDB atom name of the hydrogen that will be replaced for the linking atom of the fragment.
    :type hydrogen_to_replace: str
    :param core_atom_linker: PDB atom name of the core that is linking the fragment.
    :type core_atom_linker: str
    :param tmpl_out_path: Output path for the template modified.
    :type tmpl_out_path: str
    :return: None
    """
    lambda_to_reduce = float(step/(total_steps+1))
    templ_ini = TemplateOPLS2005(template_initial_path)
    for bond in templ_ini.list_of_bonds:
        key, bond_cont = bond
    templ_grw = TemplateOPLS2005(template_grown_path)
    fragment_atoms, core_atoms_in, core_atoms_grown = detect_atoms(template_initial=templ_ini, 
                                                                   template_grown=templ_grw,
                                                                   hydrogen_to_replace=hydrogen_to_replace)
    set_fragment_atoms(list_of_fragment_atoms=fragment_atoms)
    set_connecting_atom(template_grown=templ_grw, pdb_atom_name=core_atom_linker)
    fragment_bonds = detect_fragment_bonds(list_of_fragment_atoms=fragment_atoms, template_grown=templ_grw)
    set_fragment_bonds(list_of_fragment_bonds=fragment_bonds)
    set_linker_bond(templ_grw)
    modify_core_parameters_linearly(templ_grw, lambda_to_reduce, templ_ini)
    reduce_fragment_parameters_linearly(templ_grw, lambda_to_reduce)
    templ_grw.write_template_to_file(template_new_name=tmpl_out_path)



