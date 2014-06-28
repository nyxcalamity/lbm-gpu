#include <stdio.h>

#include "utils_gpu.h"


int HasCudaGpu(){
	int devices = 0;
	cudaError_t err = cudaGetDeviceCount(&devices);
	devices = (devices > 0 && err == cudaSuccess) ? 1 : 0;
	//TODO:consider removing printouts
//	if (devices)
//		printf("This computer has CUDA enabled GPU\n");
//	else
//		printf("This computer has no CUDA enabled GPU\n");

	return devices;
}
