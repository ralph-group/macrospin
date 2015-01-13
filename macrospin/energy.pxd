from macrospin.types cimport *


cdef inline float zeeman(float3 m, float3 h):
    """ Returns the Zeeman energy from a field """
    return -m.dot(h)