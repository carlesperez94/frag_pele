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


here = path.abspath(path.dirname(__file__))
ext_modules = []
cmdclass = {}
#cmdclass.update({'install': PreInstallCommand})

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
setup(
    name="FrAG",
    version="2.0.0",
    description='FrAG, a new tool for in silico hit-to-lead drug design, capable of growing a frag_pelement into a core while exploring the protein-ligand conformational space',
    long_description=long_description,
    url="https://github.com/carlesperez94/frag_pele/",
    author='Carles Perez Lopez, Daniel Soler Viladrich',
    author_email='daniel.soler@nostrumbiodiscovery.com, carlesperez@gmail.com',
    license='',
    packages=find_packages(exclude=['docs', 'tests']),
    include_package_data=True,
    include_dirs=[numpy.get_include()],
    install_requires=['cython', 'numpy',  'scipy', 'matplotlib', 'biopython ', 'pandas',  'prody==1.10', 'six', 'pytest'],
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

