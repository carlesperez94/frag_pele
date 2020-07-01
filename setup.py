import numpy
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import frag_pele as fp

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

here = path.abspath(path.dirname(__file__))
ext_modules = []
cmdclass = {}

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
setup(
    name="frag_pele",
    version=fp.__version__,
    description='FragPELE, a new tool for in silico hit-to-lead drug design, capable of growing a frag_pelement into a core while exploring the protein-ligand conformational space',
    long_description=long_description,
    url="https://github.com/carlesperez94/frag_pele/",
    author='Carles Perez Lopez, Daniel Soler Viladrich',
    author_email='daniel.soler@nostrumbiodiscovery.com, carlesperez94@gmail.com',
    license='',
    packages=find_packages(exclude=['docs', 'tests']),
    include_package_data=True,
    include_dirs=[numpy.get_include()],
    install_requires=['cython', 'numpy',  'scipy', 'matplotlib', 'biopython ', 'pandas',  'prody==1.10', 'pytest',
    'AdaptivePELE', 'lib_prep', 'mdtraj'],
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

