#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Kernel Cython classes for efficiently evolving the macrospin
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014-2015 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# distutils: language = c++
from macrospin.float3 cimport *


cdef inline float3 demagnetization(float3 m, float3 Nd):
    return -m*Nd


cdef inline float3 uniaxial_anisotropy(float3 m, float3 eu, float hu1, float hu2):
    cdef float m_eu = m.dot(eu)
    return eu*(m_eu*hu1) + eu*(m_eu*m_eu*m_eu*hu2)
    

cdef inline float3 cubic_anisotropy(float3 m, float3 c1, float3 c2, float3 c3,
     float hc1, float hc2):
    cdef:
        float m_c1 = m.dot(c1)
        float m_c2 = m.dot(c2)
        float m_c3 = m.dot(c3)
        float3 h
    h =     hc1*(m_c2**2 + m_c3**2)*(m_c1*c1)
    h = h + hc1*(m_c1**2 + m_c3**2)*(m_c2*c2)
    h = h + hc1*(m_c1**2 + m_c2**2)*(m_c3*c3)

    h = h + hc2*(m_c2**2 * m_c3**2)*(m_c1*c1)
    h = h + hc2*(m_c1**2 * m_c3**2)*(m_c2*c2)
    h = h + hc2*(m_c1**2 * m_c2**2)*(m_c3*c3)
    return h
