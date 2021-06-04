import os, argparse
from lib_prep.pdb_modifier import PDB


def parseArguments():
    """
        Parse the command-line options
    """
    desc = "It process PDB files putting the ligand part of a covalent ligand into the residue that it is attached with."
    parser = argparse.ArgumentParser(description=desc)
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("-i", "--input_pdb", type=str,  
                                help="Path to the input PDB complex.")
    required_named.add_argument("-rc", "--residue_chain", type=str,   
                                help="Chain of the residue that has the ligand covalently attached")
    required_named.add_argument("-rn", "--residue_number", type=str,
                                help="Residue number of the residue that has the ligand covalently attached")
    required_named.add_argument("-o", "--output_pdb", type=str,
                                help="Output PDB path.")
    parser.add_argument("-lr", "--ligand_resname", default="UNK",
                        help="Residue name of the ligand that is covalently attached to the previous residue")
    parser.add_argument("-lc", "--ligand_chain", default=" ",
                        help="Chain of the ligand that is covalently attached to the previous residue")
    parser.add_argument("-v", "--verbose", default=False, action='store_true',
                        help="Set true to activate all unneded prints.")
    parser.add_argument("-e", "--extract_ligand", default=True, action='store_false',
                        help="Set true to inactivate the extraction of a ligand in a PDB file.")
    args = parser.parse_args()

    return args.input_pdb, args.residue_chain, args.residue_number, args.output_pdb, \
           args.ligand_resname, args.ligand_chain, args.verbose, args.extract_ligand


class CovCorrector:
    def __init__(self, input_pdb, residue_chain, residue_number, ligand_resname='UNK', ligand_chain=' ',
                 verbose=False):
        """
        Class to correct PDB files with covalent ligands from Maestro's outputs.

        Parameters
        ----------
        input_pdb : str
            path to input PDB file
        residue_chain : str
            chain of the residue that has the ligand covalently attached
        residue_number : int
            number of the residue that has the ligand covalently attached
        ligand_resname : str
            residue name of the ligand that is covalently attached to the previous residue
        ligand_chain : str
            chain of the ligand that is covalently attached to the previous residue
        verbose : bool
            if true it will activate all the prints
        ----------

        """
        self._input_pdb = input_pdb
        self._residue_chain = residue_chain
        self._residue_number = residue_number
        if isinstance(self._residue_number, str):
            self._residue_number = int(self._residue_number)
        self._ligand_resname = ligand_resname
        self._ligand_chain = ligand_chain
        self._verbose = verbose
        self._pdb = PDB(in_pdb=self._input_pdb, chain=self._ligand_chain, resname=self._ligand_resname)
        self._ligand_lines = None
        self._ligand = None
        self._residue = None
        self._residue_type = None
        self._get_ligand_section()
        self._get_residue_info()

    def _get_ligand_section(self):
        """
        It gets and prints the lines of the PDB file which contains the assigned ligand_resname.
        It will fill the attribute self._ligand_lines and self._ligand
        """
        self._ligand_lines = self._pdb.get_atoms_of_resname(self._ligand_resname)
        self._ligand = ''.join(self._ligand_lines)
        if self._verbose:
            print('Ligand lines: \n' + self._ligand)

    def _get_residue_info(self):
        """
        It obtains the residue type of the residue that has the ligand covalently attached and fills
        the variable self.residue_type
        """
        residue_lines = self._pdb.get_residue(self._residue_chain, str(self._residue_number))
        self.residue = ''.join(residue_lines)
        self.residue_type = ''.join(residue_lines[0][17:20])
        if self._verbose:
            print('Residue lines: \n' + self.residue)
            print('Residue type: {}'.format(self.residue_type))

    def _correct_ligand_to_be_residue(self):
        """
        It corrects the ligand lines to transform it in to residue lines
        """
        new_ligand_lines = []
        for lig_line in self._ligand_lines:
            lig_line = list(lig_line)
            lig_line[0:6] = "ATOM  "  # Atom section
            lig_line[17:20] = self.residue_type  # Residue name
            lig_line[21:22] = self._residue_chain  # Chain
            lig_line[22:26] = "{:>4}".format(self._residue_number)  # Residue number, must be right-aligned
            new_ligand_lines.append(''.join(lig_line))
        self._ligand = ''.join(new_ligand_lines)
        self._ligand_lines = new_ligand_lines
        if self._verbose:
            print('New ligand lines:\n' + self._ligand)

    def _join_ligand_and_residue(self, reindexing=True):
        """
        It appends ligand lines onto the residue, and it can reindex atom ids.

        Parameters
        ----------
        reindexing : bool
            if it is true, it will reindex the atom IDs of the joining result
            
        Returns
        -------
        joining_result : list
            lines of the amino-acid with the ligand attached, all in the same residue
        """
        counter = 1
        joining_result = []
        for res_line in self.residue.split('\n')[0:-1]:  # The last line is empty because of the splitting
            res_line = list(res_line)
            if reindexing:
                res_line = list(res_line)
                res_line[6:11] = "{:>5}".format(counter)
                counter = counter + 1
            res_line[17:20] = self._ligand_resname
            res_line = ''.join(res_line)
            joining_result.append(res_line+'\n')
        for lig_line in self._ligand_lines:
            lig_line = list(lig_line)
            if reindexing:
                lig_line = list(lig_line)
                lig_line[6:11] = "{:>5}".format(counter)
                counter = counter + 1
            lig_line[17:20] = self._ligand_resname
            lig_line = ''.join(lig_line)
            joining_result.append(lig_line)
        return joining_result

    def correct(self, reindexing=True):
        """
        It joins the ligand and the amino-acid of a covalent ligand into a single residue, and both are moved to the
        residue position of the PDB file. This makes it compatible with PELE.

        Parameters
        ----------
        reindexing : bool
            if it is true, it will reindex the atom IDs of the joining result
        """
        global residue_idx
        new_pdb_lines = []
        self._correct_ligand_to_be_residue()
        residue_corrected = self._join_ligand_and_residue(reindexing)
        index_filled = False
        for n, line in enumerate(self._pdb.lines):
            if line.startswith('ATOM') and \
                    line[21:22] == self._residue_chain and \
                    int(line[22:26].strip()) == self._residue_number:
                if not index_filled:
                    residue_idx = n
                    index_filled = True
                continue
            if line.startswith('HETATM') and \
                    line[21:22] == self._ligand_chain and \
                    line[17:20].strip() == self._ligand_resname:
                continue
            new_pdb_lines.append(line)
        new_pdb_lines[residue_idx:residue_idx] = residue_corrected  # Insert the residue in the correct index
        self._pdb.lines = new_pdb_lines
        self._pdb.content = ''.join(new_pdb_lines)

    def write_file(self, output_file, extract_ligand=True):
        """
        It writes the content into an output PDB file.

        Parameters
        ----------
        output_file : str
            path of the output PDB file
        extract_ligand : bool
            if true it extracts the ligand in a separated PDB file
        """
        with open(output_file, "w") as out_pdb:
            out_pdb.write(self._pdb.content)
            print("PDB saved in {}.".format(output_file))
        if extract_ligand:
            lig_lines = self._pdb.get_atoms_of_resname(f"{self._ligand_resname}")
            ligand_block = "".join(lig_lines)
            out_path = os.path.dirname(output_file)
            with open(os.path.join(out_path, f"{self._ligand_resname}.pdb"), "w") as out_lig:
                out_lig.write(ligand_block)
                print("Ligand PDB saved in {}.".format(os.path.join(out_path, 
                                                       f"{self._ligand_resname}.pdb")))

def run(input_pdb, residue_chain, residue_number, out_pdb, ligand_resname='UNK', ligand_chain=' ', 
        verbose=False, extract_ligand=True):
    """
    It process PDB files putting the ligand part of a covalent ligand into the residue that it is attached with.

    Parameters
        ----------
        input_pdb : str
            path to input PDB file
        residue_chain : str
            chain of the residue that has the ligand covalently attached
        residue_number : int
            number of the residue that has the ligand covalently attached
        out_pdb : str
            path of the output PDB file
        ligand_resname : str
            residue name of the ligand that is covalently attached to the previous residue
        ligand_chain : str
            chain of the ligand that is covalently attached to the previous residue
        verbose : bool
            if true it will activate all the prints
        extract_ligand : bool
            if true it extracts the ligand in a separated PDB file
    """
    corrector = CovCorrector(input_pdb=input_pdb, residue_chain=residue_chain, residue_number=residue_number,
                             ligand_resname=ligand_resname, ligand_chain=ligand_chain, verbose=verbose)
    corrector.correct()
    corrector.write_file(output_file=out_pdb, extract_ligand=extract_ligand)


if __name__ == '__main__':
    input_pdb, residue_chain, residue_number, out, ligand_resname, ligand_chain, \
    verbose, extract_ligand = parseArguments()
    run(input_pdb, residue_chain, residue_number, out, ligand_resname, ligand_chain,
        verbose, extract_ligand)



