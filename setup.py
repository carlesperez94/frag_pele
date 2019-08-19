import subprocess
import numpy
import sys
import shutil
# subprocess.call("pip install numpy cython".split())
from setuptools import setup, find_packages, Command
from string import Template
# To use a consistent encoding
from codecs import open
from os import path
from distutils.extension import Extension
from setuptools.command.install import install
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True
from distutils.command.sdist import sdist as _sdist

# Run the following line to compile atomset package
# python setup.py build_ext --inplace
import frag_pele.constants as constants


class PreInstallCommand(install):
    description = "Installer"
    user_options = install.user_options + [
        ('schr=', None, 'SCRHODINGER main path. i.e /opt/apps/schrodinger-2017/'),
        ('pele=', None, 'PELE main path. i.e /opt/apps/PELErev1234/'),
        ('pele-exec=', None, 'PELE bin path. i.e /opt/apps/PELErev1234/bin/Pele_mpi'),
        ('pele-license=', None, 'PELE licenses PATH. i.e /opt/apps/PELErev12345/licenses/')
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.schr = None
        self.pele = None
        self.pele_exec = None
        self.pele_license = None

    def finalize_options(self):
        install.finalize_options(self)
        #if not self.schr:
        #    raise ValueError("Define --schr path. Check --help-commands for more help")
        #if not self.pele:
        #    raise ValueError("Define --pele path. Check --help-commands for more help")
        #if not self.pele_exec:
        #    raise ValueError("Define --pele-exec path. Check --help-commands for more help")
        #if not self.pele_license:
        #    raise ValueError("Define --pele-license path. Check --help-commands for more help")
        #if not self.mpirun:
        #    raise ValueError("Define --mpirun path. Check --help-commands for more help")

    def run(self):
        #print("Cythonazing")
        #subprocess.call("python frag_pele/setup.py build_ext --inplace".split())
        #print("Installing packages")
        #subprocess.call("pip install {}".format(" ".join(packages)).split())
        print("Setting environmental variables")
        installer(self.schr, self.pele, self.pele_exec, self.pele_license)
        print("Install")
        install.run(self)

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        install.run(self)


def installer(schr, pele, pele_exec, pele_license):
    file_input = 'frag_pele_pele/constants.py'
    shutil.copy('frag_pele_pele/Templates/constants.py', file_input)
    d = {"SCHRODINGER":schr, "PELE":pele, "PELE_BIN":pele_exec, "LICENSE":pele_license }
    filein = open(file_input)
    src = Template( filein.read() )
    installation_content = src.safe_substitute(d)
    filein.close()
    with open(file_input, "w") as f:
        f.write(installation_content)

        

here = path.abspath(path.dirname(__file__))
ext_modules = []
cmdclass = {}
#cmdclass.update({'install': PreInstallCommand})


class sdist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are
        # up-to-date
        from Cython.Build import cythonize
        cythonize(['cython/mycythonmodule.pyx'])
        _sdist.run(self)
        cmdclass['sdist'] = sdist

        

if use_cython:
    ext_modules += [
        Extension("frag_pele.AdaptivePELE_repo.AdaptivePELE.atomset.atomset", ["frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset/atomset.pyx"], include_dirs=["frag_pele/AdaptivePELE_repo/AdaptivePELE", "frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset"]),
        Extension("frag_pele.AdaptivePELE_repo.AdaptivePELE.atomset.SymmetryContactMapEvaluator", ["frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset/SymmetryContactMapEvaluator.pyx"], include_dirs=["frag_pele/AdaptivePELE_repo/AdaptivePELE", "frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset"]),
        Extension("frag_pele.AdaptivePELE_repo.AdaptivePELE.atomset.RMSDCalculator", ["frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset/RMSDCalculator.pyx"], include_dirs=["frag_pele/AdaptivePELE_repo/AdaptivePELE", "frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset"]),
        Extension("frag_pele.AdaptivePELE_repo.AdaptivePELE.freeEnergies.utils", ["frag_pele/AdaptivePELE_repo/AdaptivePELE/freeEnergies/utils.pyx"], include_dirs=["frag_pele/AdaptivePELE_repo/AdaptivePELE", "frag_pele/AdaptivePELE_repo/AdaptivePELE/freeEnergies"])
    ]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("frag_pele.AdaptivePELE_repo.AdaptivePELE.atomset.atomset", ["frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset/atomset.c"], include_dirs=["frag_pele/AdaptivePELE_repo/AdaptivePELE", "frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset"]),
        Extension("frag_pele.AdaptivePELE_repo.AdaptivePELE.atomset.SymmetryContactMapEvaluator", ["frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset/SymmetryContactMapEvaluator.c"], include_dirs=["frag_pele/AdaptivePELE_repo/AdaptivePELE", "frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset"]),
        Extension("frag_pele.AdaptivePELE_repo.AdaptivePELE.atomset.RMSDCalculator", ["frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset/RMSDCalculator.c"], include_dirs=["frag_pele/AdaptivePELE_repo/AdaptivePELE", "frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset"]),
        Extension("frag_pele.AdaptivePELE_repo.AdaptivePELE.freeEnergies.utils", ["frag_pele/AdaptivePELE_repo/AdaptivePELE/freeEnergies/utils.c"], include_dirs=["frag_pele/AdaptivePELE_repo/AdaptivePELE", "frag_pele/AdaptivePELE_repo/AdaptivePELE/freeEnergies"])
    ]

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
setup(
    name="frag_pele",
    version="1.1.0.4",
    description='FrAG, a new tool for in silico hit-to-lead drug design, capable of growing a frag_pelement into a core while exploring the protein-ligand conformational space',
    long_description=long_description,
    url="https://github.com/carlesperez94/frag_pele/",
    author='Carles Perez Lopez, Daniel Soler Viladrich',
    author_email='daniel.soler@nostrumbiodiscovery.com, carlesperez@gmail.com',
    license='',
    packages=find_packages(exclude=['docs', 'tests']),
    package_data={"frag_pele/AdaptivePELE_repo/AdaptivePELE/atomset": ['*.pxd'], "Templates": ["*.pdb", "*.conf"] },
    include_package_data=True,
    include_dirs=[numpy.get_include()],
    install_requires=['cython', 'numpy',  'scipy', 'matplotlib', 'biopython ', 'pandas',  'prody', 'six'],
    cmdclass=cmdclass,
    ext_modules=ext_modules,  # accepts a glob pattern
    #include_dirs=[numpy.get_include()],
    classifiers=(
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research"
    ),
    project_urls={
    'Documentation': 'https://carlesperez94.github.io/frag_pele/',
    'Source': 'https://carlesperez94.github.io/frag_pele/',
'Tracker': 'https://github.com/carlesperez94/frag_pele/issues',
},
)

