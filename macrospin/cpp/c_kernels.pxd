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

cdef extern from "cpp/kernels.h":
    cdef cppclass Step:
        float3 m, torque
        float t, eps
    cdef cppclass Kernel:
        Step previous, current
        void (*step_func)(Kernel*);
        float dt;
        void evolve(float*, long, long)
        void set_step_func(char*)
        float3 field(float, float3)
        float3 torque(float, float3)
        float energy(float, float3)
    ctypedef void (*pStepFunc)(Kernel*)
    cdef cppclass BasicKernel:
        Step previous, current
        pStepFunc step_func;
        float dt, alpha;
        float3 hext, Nd;
        void evolve(float*, long, long)
        void set_step_func(char*)
        float3 field(float, float3)
        float3 torque(float, float3)
        float energy(float, float3)
