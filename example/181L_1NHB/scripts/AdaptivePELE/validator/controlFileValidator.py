import validatorBlockNames
from AdaptivePELE.constants import blockNames
import json
import warnings
import numbers
import sys


def validate(control_file):
    isCorrect = True
    with open(control_file, 'r') as f:
        jsonFile = f.read()
    try:
        print jsonFile
        parsedJSON = json.loads(jsonFile)
    except ValueError:
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
        print "Congratulations! No errors found in your control file!"
    else:
        raise ValueError("There are errors in your control file!!!")
    return True


def validateBlock(blockName, controlFileBlock):
    isCorrect = True
    blockType = controlFileBlock["type"]
    # Check if type selected is valid
    if not isinstance(blockType, basestring):
        warnings.warn("Type for %s should be %s and instead is %s" %
                      (blockType, 'basestring', type(blockType).__name__))
        isCorrect = False

    # check for mandatory parameters
    try:
        for mandatory, value in blockName.types[blockType].iteritems():
            try:
                if not isinstance(controlFileBlock['params'][mandatory], eval(value)):
                    warnings.warn("Type for %s should be %s and instead is %s" %
                                  (mandatory, value, type(controlFileBlock['params'][mandatory]).__name__))
                    isCorrect = False
            except KeyError as err:
                warnings.warn("%s missing: Mandatory parameter %s in %s not found." %
                              (err.message, mandatory, blockName.__name__))
                isCorrect = False
    except KeyError as err:
        warnings.warn("Missing %s: Type %s in %s not found." %
                      (err.message, blockType, blockName.__name__))
        isCorrect = False
    # check rest of parameters specified
    try:
        for param, value in controlFileBlock["params"].iteritems():
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
        warnings.warn("Missing %s in %s" % (err.message, blockName.__name__))
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
                              (err.message, block))
                isCorrect = False
            blockType = controlFileBlock[block]["type"]
            if blockType not in types_dict:
                warnings.warn("Type %s in %s not found." %
                              (blockType, blockName.__name__))
                isCorrect = False
            if not isinstance(blockType, basestring):
                warnings.warn("Type for %s should be %s and instead is %s" %
                              (blockType, 'basestring', type(blockType).__name__))
                isCorrect = False
            # check rest of parameters specified
            # Do a get on the "params" block and return an empty list if not found
            paramsControlFile = controlFileBlock[block].get("params", {})
            for param, value in paramsControlFile.iteritems():
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
    isCorrect = True
    for key,value in blockName.mandatory.iteritems():
        try:
            if not isinstance(controlFileBlock[key], eval(value)):
                warnings.warn("Type for %s should be %s and instead is %s" %
                              (key, value, type(controlFileBlock[key]).__name__))
                isCorrect = False
        except KeyError:
            warnings.warn("Mandatory parameter %s in GeneralParams not found." %
                          key)
            isCorrect = False

    for key, value in controlFileBlock.iteritems():
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
            print "Validating control file %s" % contfile
            validate(contfile)
