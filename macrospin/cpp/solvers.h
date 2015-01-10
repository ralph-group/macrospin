#ifndef __SOLVERS_H_INCLUDED__
#define __SOLVERS_H_INCLUDED__

#include <functional>

#include "types.h"
#include "kernels.h"

#ifdef __cplusplus
extern "C" {
#endif

void euler_step(Kernel* kernel);
void huen_step(Kernel* kernel);
void rk23_step(Kernel* kernel);
void rk4_step(Kernel* kernel);
void rk45_step(Kernel* kernel);

#ifdef __cplusplus
}
#endif

#endif