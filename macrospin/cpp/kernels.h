#ifndef __KERNELS_H_INCLUDED__
#define __KERNELS_H_INCLUDED__

#include <functional>

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

class Kernel;

struct Step {
    /* Stores the results of a time step operation

    m - moment unit vector
    t - simulation time
    torque - torque at (t, m)
    eps - error

    */
    float3 m, torque;
    float t, eps;
};

typedef void (*pStepFunc)(Kernel*);

class Kernel {
    public:
        Step previous, current;
        pStepFunc step_func;
        float dt;
        void evolve(float *moments_ptr, const long external_steps, const long internal_steps);
        void set_step_func(char *step_func_name);
        virtual float3 field(float t, float3 m);
        virtual float3 torque(float t, float3 m);
        virtual float energy(float t, float3 m);
};


class BasicKernel: public Kernel {
    public:
        float3 hext, Nd;
        float alpha;
        float3 field(float t, float3 m);
        float3 torque(float t, float3 m);
};

#ifdef __cplusplus
}
#endif

#endif