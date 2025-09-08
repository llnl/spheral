#include "DeviceTest.hh"

#include<stdio.h>

namespace Spheral
{

#if defined(__CUDACC__) or defined(__HIPCC__)
__device__ void add(int a, int b, int *c)
{
  *c = a + b;
}

__global__ void launch(int a, int b, int *c)
{
  *c = a + b;
  //add(a,b,c);
}

#if defined(__CUDACC__)
__host__ int launchCaller(int a, int b)
{
  int c;
  int *d_c;
  cudaMalloc((void**) &d_c, sizeof(int));

  launch<<<1,1>>>(a,b,d_c);
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) 
      printf("Error: %s\n", cudaGetErrorString(err));

  cudaMemcpy(&c, d_c, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_c);
  return c;
}
#endif

#if defined(__HIPCC__)
__host__ int launchCaller(int a, int b)
{
  int c;
  int *d_c;
  (void)hipMalloc((void**) &d_c, sizeof(int));

  launch<<<1,1>>>(a,b,d_c);
  hipError_t err = hipGetLastError();
  if (err != hipSuccess) 
      printf("Error: %s\n", hipGetErrorString(err));

  (void)hipMemcpy(&c, d_c, sizeof(int), hipMemcpyDeviceToHost);
  (void)hipFree(d_c);
  return c;
}
#endif

#else
int launchCaller(int a, int b)
{
  return a + b;
}

#endif

} // namespace Spheral
