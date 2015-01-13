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
    cdef void _evolve(self, float* moments_ptr, long external_steps, long internal_steps)
    cdef float3 field(self, float t, float3 m)
    cdef float3 torque(self, float t, float3 m)
    cdef float energy(self, float t, float3 m)


cdef class BasicKernel(Kernel):
    cdef:
        public float alpha
        float3 hext, Nd
        
    cdef void _evolve(self, float* moments_ptr, long external_steps, long internal_steps)
    cdef float3 field(self, float t, float3 m)
    cdef float3 torque(self, float t, float3 m)
    cdef float energy(self, float t, float3 m)
