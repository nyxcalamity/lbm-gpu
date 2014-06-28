#include <math.h>
#include <stdio.h>

#include "lbm_model.h"
#include "collision_gpu.h"
#include "utils_gpu.h"


__constant__ float tau_d;
__constant__ int xlength_d;

/**
 * Computes the density from the particle distribution functions stored at currentCell.
 * currentCell thus denotes the address of the first particle distribution function of the
 * respective cell. The result is stored in density.
 */
__device__ void ComputeDensityGpu(float *current_cell, float *density){
    int i; *density=0;
    for(i=0;i<Q_LBM;i++)
        *density+=current_cell[i];

    /* Density should be close to a unit (ρ~1) */
//    if((*density-1.0)>EPS)
//        ERROR("Density dropped below error tolerance.");
}


/**
 * Computes the velocity within currentCell and stores the result in velocity
 */
__device__ void ComputeVelocityGpu(float *current_cell, float *density, float *velocity){
    int i;
    velocity[0]=0;
    velocity[1]=0;
    velocity[2]=0;

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
    int i;
    float s1, s2, s3;
    for(i=0;i<Q_LBM;i++){
        s1 = LATTICE_VELOCITIES_D[i][0]*velocity[0]+LATTICE_VELOCITIES_D[i][1]*velocity[1]+
        		LATTICE_VELOCITIES_D[i][2]*velocity[2];
        s2 = s1*s1;
        s3 = velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2];

        feq[i]=LATTICE_WEIGHTS_D[i]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

        /* Probability distribution function can not be less than 0 */
//        if (feq[i] < 0)
//            ERROR("Probability distribution function can not be negative.");
    }
}


/**
 * Computes the post-collision distribution functions according to the BGK update rule and
 * stores the results again at the same position.
 */
__device__ void ComputePostCollisionDistributionsGpu(float *current_cell, float *feq){
    int i;
    for(i=0;i<Q_LBM;i++){
        current_cell[i]=current_cell[i]-(current_cell[i]-feq[i])/tau_d;

        /* Probability distribution function can not be less than 0 */
//        if (current_cell[i] < 0)
//            ERROR("Probability distribution function can not be negative.");
    }
}


/**
 * Performs the actual collision computation
 */
__global__ void DoColision(float *collide_field_d){
	//	__syncthreads(); to use after reading data into shared memory
	float density, velocity[D_LBM], feq[Q_LBM], *currentCell;
	int x = 1+threadIdx.x+blockIdx.x*blockDim.x;
	int y = 1+threadIdx.y+blockIdx.y*blockDim.y;
	int z = 1+threadIdx.z+blockIdx.z*blockDim.z;
	int step = xlength_d+2;
	int idx = x+y*step+z*step*step;

	//check that indices are within the bounds since there could be more threads than needed
	if (x<(step-1) && y<(step-1) && z<(step-1)){
		currentCell=&collide_field_d[Q_LBM*idx];
		ComputeDensityGpu(currentCell,&density);
		ComputeVelocityGpu(currentCell,&density,velocity);
		ComputeFeqGpu(&density,velocity,feq);
		ComputePostCollisionDistributionsGpu(currentCell,feq);
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
