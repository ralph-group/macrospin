#include "kernels.h"
#include "field.h"
#include "energy.h"
#include "torque.h"
#include <string.h> 
#include "solvers.h"

#ifdef __cplusplus
extern "C" {
#endif

void Kernel::evolve(float *moments_ptr, const long external_steps, const long internal_steps) {
    /* Takes steps

    **moments - pointer to array of array of floats
    external_steps -
    internal_steps -

    */

    for(long i = 0; i < external_steps; i++) {
        for(long j = 0; j < internal_steps; j++) {
            this->previous = this->current;
            this->step_func(this);
        }
        moments_ptr[3*i] = this->current.m.x;
        moments_ptr[3*i+1] = this->current.m.y;
        moments_ptr[3*i+2] = this->current.m.z;             
    }

}

void Kernel::set_step_func(char *step_func_name) {
    pStepFunc step_func;
    if (strcmp(step_func_name, "Euler") == 0) {
        step_func = &euler_step;
    } else if (strcmp(step_func_name, "RK23") == 0) {
        step_func = &rk23_step;
    }
    this->step_func = step_func;
}

float3 Kernel::field(float t, float3 m) {
    float3 field;
    return field;
}


float3 Kernel::torque(float t, float3 m) {
    float3 torque;
    return torque;
}

float Kernel::energy(float t, float3 m) {
    /* Returns the free energy for a given moment and time

    m - moment unit vector
    t - simulation time

    */
    return energy::zeeman(m, this->field(t, m));
}




float3 BasicKernel::field(float t, float3 m) {
    /* Returns the effective field for a given moment and time

    m - moment unit vector
    t - simulation time

    */
    return this->hext + field::demagnetization(m, this->Nd);
}

float3 BasicKernel::torque(float t, float3 m) {
    /* Returns the torque for a given moment and time

    m - moment unit vector
    t - simulation time

    */
    float3 heff = this->field(t, m);
    return torque::landau_lifshitz(m, heff, this->alpha);
}

#ifdef __cplusplus
}
#endif