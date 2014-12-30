#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Kernel Cython classes for efficiently evolving the macrospin
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# distutils: language = c++
cimport cython
import numpy as np
cimport numpy as np
from libcpp cimport queue

cdef extern from "float3.h":
    cdef cppclass float3:
        float3(x, y, z)
        float3()
        float x, y, z
        float3 operator+(float3)
        float3 operator-(float3)
        float3 operator*(float3)
        float3 operator*(float)
        float3 operator/(float)
        float3 cross(float3)
        float dot(float3)
        float mag()
        void normalize()

@cython.boundscheck(False)
@cython.wraparound(False)
cdef float3 make_float3(float[::1] l):
    cdef float3 r
    r.x, r.y, r.z = l[0], l[1], l[2]
    return r


class Kernel(object):

    def __call__(self, time=None, timeout=None):
        """ Run the simulation for a given time, or until the moment
        stabilizes (None)
        """
        self._stopped = False
        pass

    def stop(self):
        """ Stops the kernel by setting the flag
        """
        self._stopped = True


class BasicKernel(Kernel):
    pass