#ifndef _UTILS_GPU_H_
#define _UTILS_GPU_H_

#include <stdio.h>

//restricts 3D blocks to have 512 threads (limits: 512 CC<2.x; 1024 CC>2.x)
#define BLOCK_SIZE 8

/**
 * Checks the returned cudaError_t and prints corresponding message in case of error.
 */
#define cudaErrorCheck(ans){ cudaAssert((ans), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, char *file, int line, bool abort=true){
	if (code != cudaSuccess){
		fprintf(stderr,"CUDA Assert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}


//TODO:move lattice constants to gpu constant memory
__device__ static const int LATTICE_VELOCITIES_D[19][3] = {
    {0,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},{0,1,-1},{-1,-1,0},{0,-1,0},{1,-1,0},
    {-1,0,0}, {0,0,0},  {1,0,0}, {-1,1,0},{0,1,0}, {1,1,0},  {0,-1,1},{-1,0,1},
    {0,0,1},  {1,0,1},  {0,1,1}
};
__device__ static const float LATTICE_WEIGHTS_D[19] = {
    1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0,
    2.0/36.0, 12.0/36.0,2.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
    2.0/36.0, 1.0/36.0, 1.0/36.0
};


//TODO:comment what these are and why we need them here
#define LEFT_BOUNDARY_IDX 0
#define RIGHT_BOUNDARY_IDX 1
#define BOTTOM_BOUNDARY_IDX 2
#define TOP_BOUNDARY_IDX 3
#define BACK_BOUNDARY_IDX 4
#define FRONT_BOUNDARY_IDX 5


/**
 * This double array used to store number of pdf's which we need to copy
 * on pdf's in treatBoundary step.
 */
__device__ static const float treat_boundary_indeces[6][5] = {
	{3,7,10,13,17},
	{1,5,8,11,15},
	{4,11,12,13,18},
	{0,5,6,7,14},
	{14,15,16,17,18},
	{0,1,2,3,4}
};


/** Determines if the computer has CUDA enabled GPU and returns true or false */
int HasCudaGpu();

#endif
