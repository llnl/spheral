#include <iostream>
#include "DeviceTestLib/DeviceTest.hh"

int main() {
  int a,b,c;
  a = 3; b = 4;

#if 1
  c = Spheral::launchCaller(a,b);
#else
  int *d_c;
  cudaMalloc((void**) &d_c, sizeof(int));
  launch<<<1,1>>>(a,b,d_c);
  cudaMemcpy(&c, d_c, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_c);
#endif

  std::cout << "C : " << c << "\n";
  
  return EXIT_SUCCESS;
}
