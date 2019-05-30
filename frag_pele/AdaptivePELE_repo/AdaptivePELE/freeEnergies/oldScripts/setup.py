"""
build with python setup.py build_ext --inplace
"""
from distutils.core import setup
from Cython.Build import cythonize

import numpy as np

setup(
    ext_modules=cythonize("revTransitionMatrix.pyx"),
    include_dirs=[np.get_include()]
)

