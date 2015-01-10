#ifndef __TORQUE_H_INCLUDED__
#define __TORQUE_H_INCLUDED__

#include "types.h"

namespace torque {

inline float3 landau_lifshitz(float3 m, float3 heff, float alpha) {
	/* Returns the Landau-Lifshitz torque
	
	m - moment unit vector
	heff - effective magnetic field
	alpha - Gilbert damping parameter

	*/
	float3 hxm = heff.cross(m);
	return hxm + alpha*m.cross(hxm);
}

inline float3 slonczewski(float3 m, float3 Jc, float stt) {
	/* Returns the Slonczewski spin-transfer torque

	m - moment unit vector
	Jc - current density vector
	stt - torque prefactor (pre-calculated)

	*/
	float3 p = -Jc*stt;
	return m.cross(p.cross(m));
}

}

#endif