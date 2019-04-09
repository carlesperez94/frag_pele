#!/bin/bash

echo "Schrodinger Path (/opt/schrodinger2018_2/):"
read SCHR

echo "PELE Path (/opt/PELErev1245/):"
read PELE

echo "PELE Binary Path (/opt/PELErev1245/bin):"
read PELE_BIN

echo "PELE Licenses Path (/opt/PELErev1245/licenses):"
read PELE_LICENSES

python FrAG_PELE/FrAG/Helpers/helpers.py --schr $SCHR --pele $PELE --pele_exec $PELE_BIN --pele_license $PELE_LICENSES

pip install FrAG_PELE/
