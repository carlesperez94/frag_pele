from __future__ import absolute_import, division, print_function, unicode_literals
from AdaptivePELE.validator import validatorBlockNames
from AdaptivePELE.constants import blockNames
import json
import warnings
import numbers
import sys
import ast
try:
    # Check if the basestring type if available, this will fail in python3
    basestring
except NameError:
    basestring = str


def validate(control_file):
    """
        Validate an AdaptivePELE control file to ensure that there are no errors

        :param controlFile: Adaptive sampling control file
        :type controlFile: str

        :raise ValueError: If a error is found in the control file
    """
    isCorrect = True
    with open(control_file, 'r') as f:
        jsonFile = f.read()
    try:
        parsedJSON = json.loads(jsonFile)
    except ValueError:
        print(jsonFile)
        raise ValueError("Invalid JSON file!")

    for block in dir(validatorBlockNames.ControlFileParams):
        if block.startswith('__'):
            continue
        block_obj = getattr(validatorBlockNames,
                            eval("validatorBlockNames.ControlFileParams.%s" %
                                                         block))
        controlfile_obj = eval("blockNames.ControlFileParams.%s" % block)

        if block == "generalParams":
            try:
                blockCorrect = validateGeneralBlock(block_obj,
                                                    parsedJSON[controlfile_obj])
                isCorrect = isCorrect and blockCorrect
            except KeyError:
                warnings.warn("Block %s not found in control file!" %
                              controlfile_obj)
                isCorrect = False
        else:
            try:
                blockCorrect = validateBlock(block_obj,
                                             parsedJSON[controlfile_obj])
                isCorrect = isCorrect and blockCorrect
            except KeyError:
                warnings.warn("Block %s not found in control file!" %
                              controlfile_obj)
                isCorrect = False

    if isCorrect:
        print("Congratulations! No errors found in your control file!")
    else:
        print(jsonFile)
        raise ValueError("There are errors in your control file!!!")
    return True


def validateBlock(blockName, controlFileBlock):
    """
        Validate the a block of the control file to ensure that
        there are no errors. Raise a warning if an error is found.

        :param blockName: Dictionary containing the parameters and possible
            types for the block
        :type blockName: dict
        :param controlFileBlock: Dictionary containing the parameters specified
            in the control file for the block
        :type controlFileBlock: dict

        :returns: bool -- Wheter and error has been found in a block
            block
    """
    isCorrect = True
    try:
        blockType = controlFileBlock["type"]
    except KeyError:
        isCorrect = False
        warnings.warn("Missing mandatory parameter type in block")

    # Check if type selected is valid
    if not isinstance(blockType, basestring):
        warnings.warn("Type for %s should be %s and instead is %s" %
                      (blockType, 'str', type(blockType).__name__))
        isCorrect = False

    # check for mandatory parameters
    try:
        for mandatory, value in blockName.types[blockType].items():
            try:
                if not isinstance(controlFileBlock['params'][mandatory], eval(value)):
                    warnings.warn("Type for %s should be %s and instead is %s" %
                                  (mandatory, value, type(controlFileBlock['params'][mandatory]).__name__))
                    isCorrect = False
            except KeyError as err:
                warnings.warn("%s missing: Mandatory parameter %s in %s not found." %
                              (str(err), mandatory, blockName.__name__))
                isCorrect = False
    except KeyError as err:
        warnings.warn("Missing %s: Type %s in %s not found." %
                      (str(err), blockType, blockName.__name__))
        isCorrect = False
    # check rest of parameters specified
    try:
        for param, value in controlFileBlock["params"].items():
            try:
                if not isinstance(value, eval(blockName.params[param])):
                    warnings.warn("Type for %s should be %s and instead is %s" %
                                  (param, blockName.params[param],
                                   type(value).__name__))
                    isCorrect = False
            except KeyError:
                warnings.warn("Parameter %s in block %s not recognized." %
                              (param, blockName.__name__))
                isCorrect = False
    except KeyError as err:
        warnings.warn("Missing %s in %s" % (str(err), blockName.__name__))
        isCorrect = False

    for block in dir(blockName):
        if not block.startswith('__') and block not in ["params", "types"]:
            # The parameters blocks for density and threshold calculator are
            # not mandatory
            if block not in controlFileBlock:
                continue
            try:
                types_dict = eval("blockName.%s" % block)["types"]
                params_dict = eval("blockName.%s" % block)["params"]
            except KeyError as err:
                warnings.warn("Type %s in %s not found." %
                              (str(err), block))
                isCorrect = False
            try:
                blockType = controlFileBlock[block]["type"]
            except KeyError:
                isCorrect = False
                warnings.warn("Missing mandatory parameter type in block %s" %
                              block)
            if blockType not in types_dict:
                warnings.warn("Type %s in %s not found." %
                              (blockType, blockName.__name__))
                isCorrect = False
            if not isinstance(blockType, basestring):
                warnings.warn("Type for %s should be %s and instead is %s" %
                              (blockType, 'str', type(blockType).__name__))
                isCorrect = False
            # check rest of parameters specified
            # Do a get on the "params" block and return an empty list if not found
            paramsControlFile = controlFileBlock[block].get("params", {})
            for param, value in paramsControlFile.items():
                try:
                    if not isinstance(value, eval(params_dict[param])):
                        warnings.warn("Type for %s should be %s and instead is %s" %
                                      (param, params_dict[param], type(value).__name__))
                        isCorrect = False
                except KeyError:
                    warnings.warn("Parameter %s not recognized." % param)
                    isCorrect = False

    return isCorrect


def validateGeneralBlock(blockName, controlFileBlock):
    """
        Validate the generalParams block of the control file to ensure that
        there are no errors. Raise a warning if an error is found.

        :param blockName: Dictionary containing the parameters and possible
            types for the block
        :type blockName: dict
        :param controlFileBlock: Dictionary containing the parameters specified
         in the control file for the block
        :type controlFileBlock: dict

        :returns: bool -- Wheter and error has been found in the generalParams
            block
    """
    isCorrect = True
    for key, value in blockName.mandatory.items():
        try:
            if not isinstance(controlFileBlock[key], eval(value)):
                warnings.warn("Type for %s should be %s and instead is %s" %
                              (key, value, type(controlFileBlock[key]).__name__))
                isCorrect = False
        except KeyError:
            warnings.warn("Mandatory parameter %s in GeneralParams not found." %
                          key)
            isCorrect = False

    for key, value in controlFileBlock.items():
        try:
            if not isinstance(value, eval(blockName.params[key])):
                warnings.warn("Type for %s should be %s and instead is %s" %
                              (key, value, type(blockName.params[key]).__name__))
                isCorrect = False
        except KeyError:
            warnings.warn("Parameter %s in GeneralParams not recognized." %
                          key)
            isCorrect = False
    return isCorrect

if __name__ == "__main__":

    if len(sys.argv) > 1:
        validate(sys.argv[1])
    else:
        controlFiles = ["tests/data/3ptb_data/integrationTest%i.conf" % i for i in range(1, 4)]
        for contfile in controlFiles:
            print("Validating control file %s" % contfile)
            validate(contfile)
