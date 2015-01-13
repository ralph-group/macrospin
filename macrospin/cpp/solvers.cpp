#include "solvers.h"
#include "kernels.h"

#ifdef __cplusplus
extern "C" {
#endif

void euler_step(Kernel* kernel) {
	/* Takes one step using the Euler method on a Kernel

	kernel - the Kernel object (or child object)

	*/
	kernel->current.t = kernel->previous.t + kernel->dt;

	kernel->current.torque = kernel->torque(kernel->previous.t, kernel->previous.m);
	kernel->current.m = kernel->previous.m + kernel->dt*kernel->current.torque;
	kernel->current.m.normalize();
	kernel->current.eps = 0.0; // TODO: Include error
}


void huen_step(Kernel* kernel) {
	/* Takes one step using Huen's method

	kernel - the Kernel object (or child object)

	*/
	kernel->current.t = kernel->previous.t + kernel->dt;

	float3 k1 = kernel->previous.torque;
	float3 k2 = kernel->torque(kernel->current.t, kernel->previous.m + kernel->dt*k1);

	kernel->current.m = kernel->previous.m + kernel->dt*(k1 + k2)/2.0;
	kernel->current.m.normalize();
	kernel->current.eps = 0.0; // TODO: Include error
}



void rk23_step(Kernel* kernel) {
	/* Takes one step using the Bogacki-Shampine method (Runga-Kutta RK23) on a Kernel

	kernel - the Kernel object (or child object)

	*/
	kernel->current.t = kernel->previous.t + kernel->dt;

	float3 k1 = kernel->previous.torque;
	float3 k2 = kernel->torque(kernel->previous.t + kernel->dt/2.0, kernel->previous.m + kernel->dt*k1/2.0);
	float3 k3 = kernel->torque(kernel->previous.t + 3.0*kernel->dt/2.0, kernel->previous.m + 3.0*kernel->dt*k2/2.0);

	kernel->current.m = kernel->previous.m + 2.0*kernel->dt*k1/9.0 + kernel->dt*k2/3.0 + 4*kernel->dt*k3/9.0;
	kernel->current.m.normalize();
	kernel->current.torque = kernel->torque(kernel->current.t, kernel->current.m);
	kernel->current.eps = 0.0; // TODO: Include error
}


void rk4_step(Kernel* kernel) {
	/* Takes one step using the Classic 4th order Runga-Kutta method

	kernel - the Kernel object (or child object)

	*/
	kernel->current.t = kernel->previous.t + kernel->dt;

	float3 k1 = kernel->previous.torque
	float3 k2 = kernel->torque(kernel->previous.t + kernel->dt/2.0, kernel->previous.m + kernel->dt*k1/2.0);
	float3 k3 = kernel->torque(kernel->previous.t + kernel->dt/2.0, kernel->previous.m + kernel->dt*k2/2.0);
	float3 k4 = kernel->torque(kernel->current.t, kernel->previous.m + kernel->dt*k3);

	kernel->current.m = kernel->previous.m + kernel->dt*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
	kernel->current.m.normalize();
	kernel->current.torque = kernel->torque(kernel->current.t, kernel->current.m);
	kernel->current.eps = 0.0; // TODO: Include error
}

void rk45_step(Kernel* kernel) {
	/* Takes one step using the Dormand-Prince method

	kernel - the Kernel object (or child object)

	*/
	// TODO: Finished RK45
}


#ifdef __cplusplus
}
#endif