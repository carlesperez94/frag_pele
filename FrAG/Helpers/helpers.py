import os
from string import Template
import shutil
import argparse


def create_symlinks(src, dst):
    if not os.path.islink(dst):
        os.symlink(src, dst)


def installer(schr, pele, pele_exec, pele_license):
    file_input = 'FrAG_PELE/FrAG/constants.py'
    shutil.copy('FrAG_PELE/FrAG/Templates/constants.py', file_input)
    d = {"SCHRODINGER":schr, "PELE":pele, "PELE_BIN":pele_exec, "LICENSE":pele_license }
    filein = open(file_input)
    src = Template( filein.read() )
    installation_content = src.safe_substitute(d)
    filein.close()
    with open(file_input, "w") as f:
        f.write(installation_content)


def parse_args():
    parser = argparse.ArgumentParser(description='Installation Script')
    parser.add_argument('--schr', type=str, help='SCRHODINGER main path. i.e /opt/apps/schrodinger-2017/'),
    parser.add_argument('--pele', type=str, help='PELE main path. i.e /opt/apps/PELErev1234/'),
    parser.add_argument('--pele_exec', type=str, help='PELE bin path. i.e /opt/apps/PELErev1234/bin/Pele_mpi'),
    parser.add_argument('--pele_license', type=str, help='PELE licenses PATH. i.e /opt/apps/PELErev12345/licenses/')
    args = parser.parse_args()
    return args.schr, args.pele, args.pele_exec, args.pele_license

if __name__ == "__main__":
    schr, pele, pele_exec, pele_license = parse_args()
    installer(schr, pele, pele_exec, pele_license)
