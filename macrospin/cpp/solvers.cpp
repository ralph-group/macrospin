#include "solvers.h"
#include "kernels.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void euler_step(Kernel* kernel) {
	/* Takes one step using the Euler method on a Kernel

	kernel - the Kernel object (or child object)

	*/
	kernel->current.torque = kernel->torque(kernel->previous.t, kernel->previous.m);
	kernel->current.m = kernel->previous.m + kernel->dt*kernel->current.torque;
	kernel->current.t = kernel->previous.t + kernel->dt;
	kernel->current.m.normalize();
	kernel->current.eps = 0.0; // TODO: Include error
}


void rk23_step(Kernel* kernel) {
	/* Takes one step using the Bogacki-Shampine method (Runga-Kutta RK23) on a Kernel

	kernel - the Kernel object (or child object)

	*/

	float3 k1 = kernel->previous.torque;
	float3 k2 = kernel->torque(kernel->previous.t + kernel->dt/2.0, kernel->previous.m + kernel->dt*k1/2.0);
	float3 k3 = kernel->torque(kernel->previous.t + 3.0*kernel->dt/2.0, kernel->previous.m + 3.0*kernel->dt*k2/2.0);

	kernel->current.m = kernel->previous.m + 2.0*kernel->dt*k1/9.0 + kernel->dt*k2/3.0 + 4*kernel->dt*k3/9.0;
	kernel->current.m.normalize();
	kernel->current.t = kernel->previous.t + kernel->dt;
	kernel->current.torque = kernel->torque(kernel->current.t, kernel->current.m);
	kernel->current.eps = 0.0; // TODO: Include error
}

#ifdef __cplusplus
}
#endif