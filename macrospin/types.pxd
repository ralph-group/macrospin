#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cython imports for C++ macrospin types
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014-2015 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# distutils: language = c++
cimport cython
import numpy as np
cimport numpy as np

cdef extern from "_types.h":
    cdef cppclass float3:
        float3(x, y, z)
        float3()
        float x, y, z
        float3 operator+(float3)
        float3 operator-()
        float3 operator-(float3)
        float3 operator*(float3)
        float3 operator*(float)
        float3 operator/(float)
        float3 cross(float3)
        float dot(float3)
        float mag()
        void normalize()
    cdef:
        float3 operator*(float, float3)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline float3 make_float3(float[::1] l):
    cdef float3 r
    r.x, r.y, r.z = l[0], l[1], l[2]
    return r


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.ndarray make_array(float3 f):
    return np.array([f.x, f.y, f.z], dtype=np.float32)