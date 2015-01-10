#ifndef __FIELD_H_INCLUDE__
#define __FIELD_H_INCLUDE__

#include "types.h"

namespace field {

inline float3 demagnetization(float3 m, float3 Nd) {
	/* Returns the demagnetization field based on the diagonal
	elements of the demagnetization tensor

	m - moment unit vector
	Nd - diagonal elements of demagnetization tensor

	*/
	return -m*Nd;
}

inline float3 uniaxial_anisotropy(float3 m, float3 u, float hu1, float hu2) {
	/* Returns the uniaxial anisotropy field

	m - moment unit vector
	u - uniaxial anisotropy unit vector
	hu1 - normalized uniaxial anisotropy field (1st order)
	hu2 - normalized uniaxial anisotropy field (2nd order)

	*/
	float m_u = m.dot(u);
	return u*(m_u*hu1) + u*(m_u*m_u*m_u*hu2);
}

inline float3 cubic_anisotropy(float3 m, float3 c1, float3 c2, float3 c3,
	float hc1, float hc2) {
	/* Returns the cubic anisotropy field

	m - moment unit vector
	c1, c2, c3 - orthogonal cubic axis unit vectors
	hc1 - normalized cubic anisotropy field (1st order)
	hc2 - normalized cubic anisotropy field (2nd order)

	*/
    float m_c1 = m.dot(c1);
    float m_c2 = m.dot(c2);
    float m_c3 = m.dot(c3);
    float3 h;
    h =  hc1*(m_c2*m_c2 + m_c3*m_c3)*(m_c1*c1);
    h += hc1*(m_c1*m_c1 + m_c3*m_c3)*(m_c2*c2);
    h += hc1*(m_c1*m_c1 + m_c2*m_c2)*(m_c3*c3);

    h += hc2*(m_c2*m_c2 * m_c3*m_c3)*(m_c1*c1);
    h += hc2*(m_c1*m_c1 * m_c3*m_c3)*(m_c2*c2);
    h += hc2*(m_c1*m_c1 * m_c2*m_c2)*(m_c3*c3);
    return h;
}

}

#endif