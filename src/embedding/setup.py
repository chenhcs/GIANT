from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize("encoder/gensimmod/model/word2vec_inner.pyx"),
    include_dirs=[numpy.get_include()]
)
