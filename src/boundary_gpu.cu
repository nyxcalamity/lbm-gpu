#include "boundary_gpu.h"
#include "utils_gpu.h"
#include "utils.h"

// TODO: put it in one file
//computation enhancing values
#define EPS 0.05
#define C_S_POW2_INV 3.0
#define C_S_POW4_INV 9.0


__constant__ int xlength_d;
__constant__ float wall_velocity_d[3];


// TODO: rename in ComputeDensityGpu
__device__ void ComputeDensityGpu2(float *current_cell, float *density){
    int i; *density=0;
    for(i=0;i<Q_LBM;i++)
        *density+=current_cell[i];

    /* Density should be close to a unit (ρ~1) */
//    if((*density-1.0)>EPS)
//        ERROR("Density dropped below error tolerance.");
}


// TODO: rename in inv
__device__ int inv2(int i){
    return (Q_LBM-1)-i;
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
						ComputeDensityGpu2(&collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)],&density);
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


void TreatBoundaryGpu(float *collide_field, int *flag_field, float *wall_velocity, int xlength){
	float *collide_field_d=NULL;
	int *flag_field_d=NULL;
	int num_cells;
	size_t size;
	float data[3];

	for(int i=0;i<3;i++)
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
	dim3 grid((xlength+block.x-1)/block.x, (xlength+block.y-1)/block.y, (xlength+block.z-1)/block.z);

	TreatBoundary<<<grid,block>>>(collide_field_d, flag_field_d);
	cudaErrorCheck(cudaPeekAtLastError());

	cudaMemcpy(collide_field, collide_field_d, size, cudaMemcpyDeviceToHost);
	cudaFree(collide_field_d);
	cudaFree(flag_field_d);
}
