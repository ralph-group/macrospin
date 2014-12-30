#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __shared__
#endif

#include <stdio.h>
#include <pycuda-complex.hpp>

// ============================================================
//                      Helper Functions
// ============================================================
inline __device__ float4 cross(float4 left, float4 right) {
  return make_float4(left.y*right.z - left.z*right.y,
		     left.z*right.x - left.x*right.z,
		     left.x*right.y - left.y*right.x,
		     0.0f);
}
inline __device__ float4 operator+(float4 a, float4 b) {
  return make_float4(a.x+b.x, a.y+b.y, a.z+b.z, 0.0f);
}
inline __device__ float4 mult(float4 a, float4 b) {
  return make_float4(a.x*b.x, a.y*b.y, a.z*b.z, 0.0f);
}
inline __device__ float4 operator-(float4 a, float4 b) {
  return make_float4(a.x-b.x, a.y-b.y, a.z-b.z, 0.0f);
}
inline __device__ float4 operator*(float b, float4 a) {
  return make_float4(a.x*b, a.y*b, a.z*b, 0.0f);
}
inline __device__ float magInv(float4 a) {
  return rsqrtf(a.x*a.x + a.y*a.y + a.z*a.z);
}

// ============================================================
//                      Physics Code
// ============================================================

__global__ void evolve(float4 *m) {

  {{ definitions }}
  {{ index_operations }}

  {{ static_field }}
  {{ demagnetization }}
  {{ uniaxial_anisotropy }}
  {{ cubic_anisotropy }}

  // Start Landau-Lifshitz equation
  float4 hxm = cross(heff, mloc);

  {{ spin_transfer_torque }}

  float4 mxhxm =  cross(mloc, hxm);

  // Compute new moment
  m[i] = mloc + {{ dt }}*(hxm + {{ damping }}*mxhxm);
}

__global__ void normalize(float4 *m) {
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  
  float4 mloc = m[i];
  m[i] = magInv(mloc)*mloc;
}