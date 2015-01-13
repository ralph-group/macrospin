from macrospin.types cimport *


cdef inline float3 landau_lifshitz(float3 m, float3 heff, float alpha) {
    """ Returns the Landau-Lifshitz torque
    
    m - moment unit vector
    heff - effective magnetic field
    alpha - Gilbert damping parameter

    """
    cdef float3 hxm = heff.cross(m)
    return hxm + alpha*m.cross(hxm)


cdef inline float3 slonczewski(float3 m, float3 Jc, float stt) {
    """ Returns the Slonczewski spin-transfer torque

    m - moment unit vector
    Jc - current density vector
    stt - torque prefactor (pre-calculated)

    """
    cdef float3 p = -Jc*stt
    return m.cross(p.cross(m))
