#ifndef __ENERGY_H_INCLUDE__
#define __ENERGY_H_INCLUDE__

#include "types.h"

namespace energy {

inline float zeeman(float3 m, float3 h) {
	/* Returns the Zeeman energy of 

	*/
	return -m.dot(h);
}

}

#endif