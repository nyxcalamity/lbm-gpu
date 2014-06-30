#include <math.h>
#include <stdio.h>

#include "lbm_solver_gpu.h"
#include "lbm_model.h"
#include "utils_gpu.h"
#include "utils.h"
#include "cell_computation_gpu.cuh"

__constant__ float tau_d, wall_velocity_d[D_LBM];
__constant__ int xlength_d, num_cells_d;


/**
 * Computes the density from the particle distribution functions stored at currentCell.
 * currentCell thus denotes the address of the first particle distribution function of the
 * respective cell. The result is stored in density.
 */
__device__ void ComputeDensityGpu(float *current_cell, float *density){
    int i; *density=0;
    //TODO:get rid of this loop
    for(i=0;i<Q_LBM;i++)
        *density+=current_cell[i];
    /* TODO:Density should be close to a unit (ρ~1) */
}


/**
 * Computes the velocity within currentCell and stores the result in velocity
 */
__device__ void ComputeVelocityGpu(float *current_cell, float *density, float *velocity){
    int i;
    velocity[0]=0;
    velocity[1]=0;
    velocity[2]=0;

    //TODO:get rid of this loop
    for(i=0;i<Q_LBM;i++){
        velocity[0]+=current_cell[i]*LATTICE_VELOCITIES_D[i][0];
        velocity[1]+=current_cell[i]*LATTICE_VELOCITIES_D[i][1];
        velocity[2]+=current_cell[i]*LATTICE_VELOCITIES_D[i][2];
    }

    velocity[0]/=*density;
    velocity[1]/=*density;
    velocity[2]/=*density;
}

/**
 * Computes the equilibrium distributions for all particle distribution functions of one
 * cell from density and velocity and stores the results in feq.
 */
__device__ void ComputeFeqGpu(float *density, float *velocity, float *feq){
    int i; float s1, s2, s3;
    //TODO:get rid of this loop
    for(i=0;i<Q_LBM;i++){
        s1 = LATTICE_VELOCITIES_D[i][0]*velocity[0]+LATTICE_VELOCITIES_D[i][1]*velocity[1]+
        		LATTICE_VELOCITIES_D[i][2]*velocity[2];
        s2 = s1*s1;
        s3 = velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2];

        feq[i]=LATTICE_WEIGHTS_D[i]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);
        /* TODO:Probability distribution function can not be less than 0 */
    }
}


/**
 * Computes the post-collision distribution functions according to the BGK update rule and
 * stores the results again at the same position.
 */
__device__ void ComputePostCollisionDistributionsGpu(float *current_cell, float *feq){
    int i;
    //TODO:get rid of this loop
    for(i=0;i<Q_LBM;i++){
        current_cell[i]=current_cell[i]-(current_cell[i]-feq[i])/tau_d;

        /* TODO:Probability distribution function can not be less than 0 */
    }
}

// TODO: rename in inv
__device__ int inv2(int i){
    return (Q_LBM-1)-i;
}


/**
 * Performs the actual collision computation
 */
__global__ void DoColision(float *collide_field_d){
	float density, velocity[D_LBM], feq[Q_LBM], *current_cell_s;
	__shared__ float collide_field_s[BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE*Q_LBM];
	//TODO:can be optimized using BLOCK_SIZE constant
	int x = 1+threadIdx.x+blockIdx.x*blockDim.x;
	int y = 1+threadIdx.y+blockIdx.y*blockDim.y;
	int z = 1+threadIdx.z+blockIdx.z*blockDim.z;
	int idx_block = threadIdx.x+threadIdx.y*blockDim.x+threadIdx.z*blockDim.x*blockDim.y;
	int step = xlength_d+2, i;
	int idx_domain = x+y*step+z*step*step;

	//check that indices are within the bounds since there could be more threads than needed
	if (x<(step-1) && y<(step-1) && z<(step-1)){
		//copy current cell values into shared memory
		for(i=0;i<Q_LBM;i++)
			collide_field_s[Q_LBM*idx_block+i]=collide_field_d[Q_LBM*idx_domain+i];

		current_cell_s = &collide_field_s[Q_LBM*idx_block];
		//perform computation
		ComputeDensityGpu(current_cell_s,&density);
		ComputeVelocityGpu(current_cell_s,&density,velocity);
		ComputeFeqGpu(&density,velocity,feq);
		ComputePostCollisionDistributionsGpu(current_cell_s,feq);

		//copy data back
		for(i=0;i<Q_LBM;i++)
			collide_field_d[Q_LBM*idx_domain+i]=collide_field_s[Q_LBM*idx_block+i];
	}
}


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

			stream_field_d[Q_LBM*idx+i]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+i];
		}
	}
}

__global__ void TreatBoundary(float *collide_field_d, int* flag_field_d){
	int x = threadIdx.x+blockIdx.x*blockDim.x;
	int y = threadIdx.y+blockIdx.y*blockDim.y;
	int z = threadIdx.z+blockIdx.z*blockDim.z;
    int nx,ny,nz,i,step=xlength_d+2;
    float density,dot_prod;

    if (x<step && y<step && z<step){
		if(flag_field_d[x+y*step+z*step*step]!=FLUID){
			for(i=0;i<Q_LBM;i++){
				nx=x+LATTICE_VELOCITIES_D[i][0];
				ny=y+LATTICE_VELOCITIES_D[i][1];
				nz=z+LATTICE_VELOCITIES_D[i][2];

				/* We don't need the values outside of our extended domain */
				if(0<nx && nx<step-1 && 0<ny && ny<step-1 && 0<nz && nz<step-1){
					if (flag_field_d[x+y*step+z*step*step]==MOVING_WALL){
						/* Compute density in the neighbour cell */
						ComputeDensityGpu(&collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)],&density);
						/* Compute dot product */
						dot_prod=LATTICE_VELOCITIES_D[i][0]*wall_velocity_d[0]+
								LATTICE_VELOCITIES_D[i][1]*wall_velocity_d[1]+
								LATTICE_VELOCITIES_D[i][2]*wall_velocity_d[2];
						/* Assign the boudary cell value */
						collide_field_d[Q_LBM*(x+y*step+z*step*step)+i]=
								collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)]+
								2*LATTICE_WEIGHTS_D[i]*density*C_S_POW2_INV*dot_prod;
					}else if(flag_field_d[x+y*step+z*step*step]==NO_SLIP){
						collide_field_d[Q_LBM*(x+y*step+z*step*step)+i]=
								collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)];
					}
				}
			}
		}
    }
}


void DoCollisionGpu(float *collide_field, int *flag_field, float tau, int xlength){
	float *collide_field_d=NULL;
	int num_cells = pow(xlength+2, D_LBM);
	size_t collide_field_size = Q_LBM*num_cells*sizeof(float);

	//initialize working data
	cudaErrorCheck(cudaMalloc(&collide_field_d, collide_field_size));
	cudaErrorCheck(cudaMemcpy(collide_field_d, collide_field, collide_field_size, cudaMemcpyHostToDevice));

	//initialize constant data
	cudaErrorCheck(cudaMemcpyToSymbol(tau_d, &tau, sizeof(float), 0, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemcpyToSymbol(xlength_d, &xlength, sizeof(int), 0, cudaMemcpyHostToDevice));

	//define grid structure
	//NOTE:redundant threads for boundary cells are not accounted for
	dim3 block(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
	dim3 grid((xlength+block.x-1)/block.x, (xlength+block.y-1)/block.y, (xlength+block.z-1)/block.z);

	//perform collision
	DoColision<<<grid,block>>>(collide_field_d);
	cudaErrorCheck(cudaPeekAtLastError());

	//copy data back to host
	cudaErrorCheck(cudaMemcpy(collide_field, collide_field_d, collide_field_size, cudaMemcpyDeviceToHost));

	//free device memory
	cudaErrorCheck(cudaFree(collide_field_d));
}


void TreatBoundaryGpu(float *collide_field, int *flag_field, float *wall_velocity, int xlength){
	float *collide_field_d=NULL, data[3];
	int *flag_field_d=NULL, num_cells;
	size_t size;

	for(int i=0;i<D_LBM;i++)
		data[i]=wall_velocity[i];

	cudaMemcpyToSymbol(wall_velocity_d, data, sizeof(data), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(xlength_d, &xlength, sizeof(int), 0, cudaMemcpyHostToDevice);

	num_cells = (xlength+2)*(xlength+2)*(xlength+2);
	size = Q_LBM*num_cells*sizeof(float);

	cudaMalloc(&collide_field_d, size);
	cudaMalloc(&flag_field_d, num_cells*sizeof(int));
	cudaMemcpy(collide_field_d, collide_field, size, cudaMemcpyHostToDevice);
	cudaMemcpy(flag_field_d, flag_field, num_cells*sizeof(int), cudaMemcpyHostToDevice);

	//define grid structure
	//NOTE:redundant threads for boundary cells are not accounted for
	dim3 block(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
	dim3 grid((xlength+2+block.x-1)/block.x, (xlength+2+block.y-1)/block.y, (xlength+2+block.z-1)/block.z);

	TreatBoundary<<<grid,block>>>(collide_field_d, flag_field_d);
	cudaErrorCheck(cudaPeekAtLastError());

	cudaMemcpy(collide_field, collide_field_d, size, cudaMemcpyDeviceToHost);
	cudaFree(collide_field_d);
	cudaFree(flag_field_d);
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
