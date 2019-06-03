import os
import frag_pele.constants as cs


def check():
    if not os.path.exists(cs.SCHRODINGER_PY_PATH):
        raise OSError("Schrodinger path {} not found. Change the harcoded path under frag_pele/constants.py".format(cs.SCHRODINGER_PY_PATH))
    if not os.path.exists(cs.PATH_TO_PELE):
        raise OSError("Pele path {} not found. Change the harcoded path under frag_pele/constants.py".format(cs.PATH_TO_PELE))
    elif not os.path.exists(cs.PATH_TO_PELE_DATA):
        raise OSError("Pele data path {} not found. Change the harcoded path under frag_pele/constants.py".format(cs.PATH_TO_PELE_DATA))
    elif not os.path.exists(cs.PATH_TO_PELE_DOCUMENTS):
        raise OSError("Pele documents path {} not found. Change the harcoded path under frag_pele/constants.py".format(cs.PATH_TO_PELE_DOCUMENTS))
    elif not os.path.exists(cs.PATH_TO_LICENSE):
        raise OSError("Pele license path {} not found. Change the harcoded path under frag_pele/constants.py".format(cs.PATH_TO_LICENSE))
