from distutils.core import setup, Extension
import numpy

setup(
    ext_modules=[
        Extension("curvefitting_cython", ["curvefitting_cython.c"],
                  include_dirs=[numpy.get_include()]),
    ],
)
