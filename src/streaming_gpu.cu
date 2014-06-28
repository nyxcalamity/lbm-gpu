#include <math.h>
#include <stdio.h>

#include "lbm_model.h"
#include "streaming_gpu.h"
#include "utils_gpu.h"


__constant__ int xlength_d;


/**
 * Performs the actual streaming computation
 */
__global__ void DoStreaming(float *stream_field_d, float *collide_field_d){
	//	__syncthreads(); to use after reading data into shared memory
	int x = 1+threadIdx.x+blockIdx.x*blockDim.x;
	int y = 1+threadIdx.y+blockIdx.y*blockDim.y;
	int z = 1+threadIdx.z+blockIdx.z*blockDim.z;
	int step = xlength_d+2, idx = x+y*step+z*step*step, nx, ny, nz, i;

	//check that indices are within the bounds since there could be more threads than needed
	if (x<(step-1) && y<(step-1) && z<(step-1)){
		for(i=0;i<Q_LBM;i++){
			nx=x-LATTICE_VELOCITIES_D[i][0];
			ny=y-LATTICE_VELOCITIES_D[i][1];
			nz=z-LATTICE_VELOCITIES_D[i][2];

			stream_field_d[Q_LBM*idx+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+i];
		}
	}
}


void DoStreamingGpu(float *collide_field, float *stream_field, int *flag_field, int xlength){
	float *collide_field_d=NULL, *stream_field_d=NULL;
	int num_cells = pow(xlength+2, D_LBM);
	size_t computational_field_size = Q_LBM*num_cells*sizeof(float);

	//initialize working data
	cudaErrorCheck(cudaMalloc(&collide_field_d, computational_field_size));
	cudaErrorCheck(cudaMemcpy(collide_field_d, collide_field, computational_field_size, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMalloc(&stream_field_d, computational_field_size));
	cudaErrorCheck(cudaMemcpy(stream_field_d, stream_field, computational_field_size, cudaMemcpyHostToDevice));

	//initialize constant data
	cudaErrorCheck(cudaMemcpyToSymbol(xlength_d, &xlength, sizeof(int), 0, cudaMemcpyHostToDevice));

	//define grid structure
	//NOTE:redundant threads for boundary cells are not accounted for
	dim3 block(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
	dim3 grid((xlength+block.x-1)/block.x, (xlength+block.y-1)/block.y, (xlength+block.z-1)/block.z);

	//perform streaming
	DoStreaming<<<grid,block>>>(stream_field_d, collide_field_d);
	cudaErrorCheck(cudaPeekAtLastError());

	//copy data back to host
	cudaErrorCheck(cudaMemcpy(stream_field, stream_field_d, computational_field_size, cudaMemcpyDeviceToHost));

	//free device memory
	cudaErrorCheck(cudaFree(collide_field_d));
	cudaErrorCheck(cudaFree(stream_field_d));
}
