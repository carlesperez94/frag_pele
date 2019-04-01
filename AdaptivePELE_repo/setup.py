import numpy
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from distutils.extension import Extension
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True
from distutils.command.sdist import sdist as _sdist
import AdaptivePELE as a

# Run the following line to compile atomset package
# python setup.py build_ext --inplace


class sdist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are
        # up-to-date
        from Cython.Build import cythonize
        cythonize(['cython/mycythonmodule.pyx'])
        _sdist.run(self)
        cmdclass['sdist'] = sdist

here = path.abspath(path.dirname(__file__))
ext_modules = []
cmdclass = {}
# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

if use_cython:
    ext_modules += [
        Extension("AdaptivePELE.atomset.atomset", ["AdaptivePELE/atomset/atomset.pyx"], include_dirs=["AdaptivePELE", "AdaptivePELE/atomset"]),
        Extension("AdaptivePELE.atomset.SymmetryContactMapEvaluator", ["AdaptivePELE/atomset/SymmetryContactMapEvaluator.pyx"], include_dirs=["AdaptivePELE", "AdaptivePELE/atomset"]),
        Extension("AdaptivePELE.atomset.RMSDCalculator", ["AdaptivePELE/atomset/RMSDCalculator.pyx"], include_dirs=["AdaptivePELE", "AdaptivePELE/atomset"]),
        Extension("AdaptivePELE.freeEnergies.utils", ["AdaptivePELE/freeEnergies/utils.pyx"], include_dirs=["AdaptivePELE", "AdaptivePELE/freeEnergies"])
    ]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("AdaptivePELE.atomset.atomset", ["AdaptivePELE/atomset/atomset.c"], include_dirs=["AdaptivePELE", "AdaptivePELE/atomset"]),
        Extension("AdaptivePELE.atomset.SymmetryContactMapEvaluator", ["AdaptivePELE/atomset/SymmetryContactMapEvaluator.c"], include_dirs=["AdaptivePELE", "AdaptivePELE/atomset"]),
        Extension("AdaptivePELE.atomset.RMSDCalculator", ["AdaptivePELE/atomset/RMSDCalculator.c"], include_dirs=["AdaptivePELE", "AdaptivePELE/atomset"]),
        Extension("AdaptivePELE.freeEnergies.utils", ["AdaptivePELE/freeEnergies/utils.c"], include_dirs=["AdaptivePELE", "AdaptivePELE/freeEnergies"])
    ]

setup(
    name="AdaptivePELE",
    version="%s" % a.__version__,
    description='Enhanced sampling of molecular simulations',
    long_description=long_description,
    url="https://github.com/AdaptivePELE/AdaptivePELE",
    author='Daniel Lecina, Joan Francesc Gilabert',
    author_email='danilecina@gmail.com, cescgina@gmail.com',
    license='',
    packages=find_packages(exclude=['docs', 'tests']),
    package_data={"AdaptivePELE/atomset": ['*.pxd']},
    install_requires=['numpy', 'mdtraj', 'scipy'],
    cmdclass=cmdclass,
    ext_modules=ext_modules,  # accepts a glob pattern
    include_dirs=[numpy.get_include()],
    classifiers=(
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research"
    ),
    project_urls={
    'Documentation': 'https://adaptivepele.github.io/AdaptivePELE//',
    'Source': 'https://github.com/AdaptivePELE/AdaptivePELE',
    'Tracker': 'https://github.com/AdaptivePELE/AdaptivePELE/issues',
},

)
