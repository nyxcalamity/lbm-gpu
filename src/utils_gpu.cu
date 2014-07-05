#include <stdio.h>

#include "utils_gpu.h"


int HasCudaGpu(){
	int devices = 0;
	cudaError_t err = cudaGetDeviceCount(&devices);
	devices = (devices > 0 && err == cudaSuccess) ? 1 : 0;
	return devices;
}
