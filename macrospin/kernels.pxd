# distutils: language = c++
cimport cython
import numpy as np
cimport numpy as np

from macrospin.types cimport *

cdef struct Step:
    float3 m, torque
    float t, eps


cdef class Kernel:
    cdef:
        public object parameters
        Step previous, current
        void (*step_func)(Kernel)
        public float dt
    cdef float3 field(self, float t, float3 m)
    cdef float3 torque(self, float t, float3 m)
    cdef float energy(self, float t, float3 m)


cdef class BasicKernel(Kernel):
    cdef:
        public float alpha
        float3 hext, Nd
    cdef float3 field(self, float t, float3 m)
    cdef float3 torque(self, float t, float3 m)
    cdef float energy(self, float t, float3 m)


cdef class AnisotropyKernel(BasicKernel):
    cdef:
        public float hu1, hu2, hc1, hc2
        float3 u, c1, c2, c3
    cdef float3 field(self, float t, float3 m)
    cdef float3 torque(self, float t, float3 m)
    cdef float energy(self, float t, float3 m)