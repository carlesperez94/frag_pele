import os
import argparse
import schrodinger.structure as st


class OutputFile(object):
    
    def __init__(self, inputfile):
        self.inputfile = inputfile
        self._format = self.inputfile.split(".")[-1]

    
    def get_property(self, Property):
        return self._get_property(Property)


    def set_property(self, Property, value, typedata="r", author="fragpele"):
        self._typedata = typedata
        self._author = author
        return self._set_property(Property, value)



class MaeFile(OutputFile):
    
    
    def __init__(self, inputfile):
        OutputFile.__init__(self, inputfile)
        self._structure = next(st.StructureReader(self.inputfile))
        self._properties = self._structure.property

    def _get_property(self, Property):
        try:
            return self._properties[Property]
        except KeyError:
            raise KeyError("Property {} not found in {}".format(Property, self.inputfile))

    def _set_property(self, Property, value):
        if Property in self._structure.property:
            self._structure.property[Property] = value
        else:
            propertyname = "{}_{}_{}_{}".format(self._typedata, self._author, Property, value)
            self._structure.property[propertyname] =  value
            self._properties = self._structure.property
        return self._properties


def pdb_to_mae(pdb_inputfile, schr_path, mae_output_file=None, remove=False):
    # Get output name from file
    file_info = pdb_inputfile.split("_")
    if not mae_output_file:
        dirpath = os.path.join(os.path.dirname(pdb_inputfile))
        filename = "{}_{}.mae".format(file_info[-3], file_info[-2])
        mae_output_file = os.path.join(dirpath, filename)

    # Get properties from name
    properties = {"BindingEnergy": float(file_info[-1].replace("BindingEnergy", "").replace(".pdb", "")),
                  "epoch": int(file_info[-2].split(".")[0]), 
                  "trajectory": file_info[-2].split(".")[1]
                  }

    # Convert to mae file
    pdb_convert = os.path.join(schr_path, "utilities/pdbconvert")
    os.system("{} -ipdb {} -omae {}".format(pdb_convert, pdb_inputfile, mae_output_file))

    # Set properties to mae_file
    mae_file = MaeFile(mae_output_file)
    for Property, value in properties.items():
        if type(value) == str:
            proper = mae_file.set_property(Property, value, typedata="s")
        if type(value) == int:
            proper = mae_file.set_property(Property, value, typedata="i")
        if type(value) == float:
            proper = mae_file.set_property(Property, value, typedata="r")

    # output maestro file
    with st.StructureWriter(mae_output_file) as writer:
        writer.append(mae_file._structure)

    if remove:
        os.remove(pdb_inputfile)


def add_args(parser):
    parser.add_argument('inputfile', type=str, help="Pdb input file")
    parser.add_argument('--schr', type=str, help="schrodinger root path")
    parser.add_argument('--remove', action="store_true", help="Remove inputfile at exit")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Write mae and sdf files with certain properties')
    add_args(parser)
    args = parser.parse_args()
    pdb_to_mae(args.inputfile, args.schr, remove=args.remove)
