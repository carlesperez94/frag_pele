import numpy
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension


ext_modules = [
        Extension("atomset", ["atomset.pyx"], include_dirs = ["."]),
        Extension("SymmetryContactMapEvaluator", ["SymmetryContactMapEvaluator.pyx"], include_dirs = ["."]),
        Extension("RMSDCalculator", ["RMSDCalculator.pyx"], include_dirs = ["."])
            ]
setup(
    ext_modules = cythonize(ext_modules),  # accepts a glob pattern
    include_dirs=[numpy.get_include()]

)
