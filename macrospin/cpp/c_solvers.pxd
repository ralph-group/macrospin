#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cython imports for Kernels from C++ macrospin library
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014-2015 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# distutils: language = c++
cimport cython

from macrospin.cpp.c_types cimport *
from macrospin.cpp.c_kernels cimport Kernel

cdef extern from "cpp/solvers.h":
    cdef:
        void euler_step(Kernel*)
        void huen_step(Kernel*)
        void rk23_step(Kernel*)
        void rk4_step(Kernel*)
        void rk45_step(Kernel*)