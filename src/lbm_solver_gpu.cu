#include <math.h>
#include <stdio.h>

#include "lbm_solver_gpu.h"
#include "lbm_model.h"
#include "utils_gpu.h"
#include "utils.h"
#include "cell_computation_gpu.cuh"
#include "boundary.h"


__constant__ float tau_d, wall_velocity_d[D_LBM];
__constant__ int xlength_d, num_cells_d;
__device__ float *stream_field_d, *collide_field_d;

/**
 * Computes the post-collision distribution functions according to the BGK update rule and
 * stores the results again at the same position.
 */
__device__ void ComputePostCollisionDistributionsGpu(float *current_cell, float *feq){
	current_cell[0]=current_cell[0]-(current_cell[0]-feq[0])/tau_d;
	current_cell[1]=current_cell[1]-(current_cell[1]-feq[1])/tau_d;
	current_cell[2]=current_cell[2]-(current_cell[2]-feq[2])/tau_d;
	current_cell[3]=current_cell[3]-(current_cell[3]-feq[3])/tau_d;
	current_cell[4]=current_cell[4]-(current_cell[4]-feq[4])/tau_d;
	current_cell[5]=current_cell[5]-(current_cell[5]-feq[5])/tau_d;
	current_cell[6]=current_cell[6]-(current_cell[6]-feq[6])/tau_d;
	current_cell[7]=current_cell[7]-(current_cell[7]-feq[7])/tau_d;
	current_cell[8]=current_cell[8]-(current_cell[8]-feq[8])/tau_d;
	current_cell[9]=current_cell[9]-(current_cell[9]-feq[9])/tau_d;
	current_cell[10]=current_cell[10]-(current_cell[10]-feq[10])/tau_d;
	current_cell[11]=current_cell[11]-(current_cell[11]-feq[11])/tau_d;
	current_cell[12]=current_cell[12]-(current_cell[12]-feq[12])/tau_d;
	current_cell[13]=current_cell[13]-(current_cell[13]-feq[13])/tau_d;
	current_cell[14]=current_cell[14]-(current_cell[14]-feq[14])/tau_d;
	current_cell[15]=current_cell[15]-(current_cell[15]-feq[15])/tau_d;
	current_cell[16]=current_cell[16]-(current_cell[16]-feq[16])/tau_d;
	current_cell[17]=current_cell[17]-(current_cell[17]-feq[17])/tau_d;
	current_cell[18]=current_cell[18]-(current_cell[18]-feq[18])/tau_d;
	/* TODO:Probability distribution function can not be less than 0 */
}


//__device__ void DoStreaming(int *flag_field_d, int x, int y, int z){
__device__ void DoStreaming(int *flag_field_d, int x, int y, int z){
//	int x = threadIdx.x+blockIdx.x*blockDim.x;
//	int y = threadIdx.y+blockIdx.y*blockDim.y;
//	int z = threadIdx.z+blockIdx.z*blockDim.z;
	int step=xlength_d+2, idx=x+y*step+z*step*step, nx, ny, nz;

	//check that indices are within the bounds since there could be more threads than needed
	if (idx<num_cells_d && flag_field_d[idx]==FLUID){
		nx=x-LATTICE_VELOCITIES_D[0][0];
		ny=y-LATTICE_VELOCITIES_D[0][1];
		nz=z-LATTICE_VELOCITIES_D[0][2];
		stream_field_d[Q_LBM*idx]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)];
		nx=x-LATTICE_VELOCITIES_D[1][0];
		ny=y-LATTICE_VELOCITIES_D[1][1];
		nz=z-LATTICE_VELOCITIES_D[1][2];
		stream_field_d[Q_LBM*idx+1]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+1];
		nx=x-LATTICE_VELOCITIES_D[2][0];
		ny=y-LATTICE_VELOCITIES_D[2][1];
		nz=z-LATTICE_VELOCITIES_D[2][2];
		stream_field_d[Q_LBM*idx+2]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+2];
		nx=x-LATTICE_VELOCITIES_D[3][0];
		ny=y-LATTICE_VELOCITIES_D[3][1];
		nz=z-LATTICE_VELOCITIES_D[3][2];
		stream_field_d[Q_LBM*idx+3]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+3];
		nx=x-LATTICE_VELOCITIES_D[4][0];
		ny=y-LATTICE_VELOCITIES_D[4][1];
		nz=z-LATTICE_VELOCITIES_D[4][2];
		stream_field_d[Q_LBM*idx+4]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+4];
		nx=x-LATTICE_VELOCITIES_D[5][0];
		ny=y-LATTICE_VELOCITIES_D[5][1];
		nz=z-LATTICE_VELOCITIES_D[5][2];
		stream_field_d[Q_LBM*idx+5]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+5];
		nx=x-LATTICE_VELOCITIES_D[6][0];
		ny=y-LATTICE_VELOCITIES_D[6][1];
		nz=z-LATTICE_VELOCITIES_D[6][2];
		stream_field_d[Q_LBM*idx+6]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+6];
		nx=x-LATTICE_VELOCITIES_D[7][0];
		ny=y-LATTICE_VELOCITIES_D[7][1];
		nz=z-LATTICE_VELOCITIES_D[7][2];
		stream_field_d[Q_LBM*idx+7]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+7];
		nx=x-LATTICE_VELOCITIES_D[8][0];
		ny=y-LATTICE_VELOCITIES_D[8][1];
		nz=z-LATTICE_VELOCITIES_D[8][2];
		stream_field_d[Q_LBM*idx+8]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+8];
		nx=x-LATTICE_VELOCITIES_D[9][0];
		ny=y-LATTICE_VELOCITIES_D[9][1];
		nz=z-LATTICE_VELOCITIES_D[9][2];
		stream_field_d[Q_LBM*idx+9]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+9];
		nx=x-LATTICE_VELOCITIES_D[10][0];
		ny=y-LATTICE_VELOCITIES_D[10][1];
		nz=z-LATTICE_VELOCITIES_D[10][2];
		stream_field_d[Q_LBM*idx+10]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+10];
		nx=x-LATTICE_VELOCITIES_D[11][0];
		ny=y-LATTICE_VELOCITIES_D[11][1];
		nz=z-LATTICE_VELOCITIES_D[11][2];
		stream_field_d[Q_LBM*idx+11]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+11];
		nx=x-LATTICE_VELOCITIES_D[12][0];
		ny=y-LATTICE_VELOCITIES_D[12][1];
		nz=z-LATTICE_VELOCITIES_D[12][2];
		stream_field_d[Q_LBM*idx+12]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+12];
		nx=x-LATTICE_VELOCITIES_D[13][0];
		ny=y-LATTICE_VELOCITIES_D[13][1];
		nz=z-LATTICE_VELOCITIES_D[13][2];
		stream_field_d[Q_LBM*idx+13]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+13];
		nx=x-LATTICE_VELOCITIES_D[14][0];
		ny=y-LATTICE_VELOCITIES_D[14][1];
		nz=z-LATTICE_VELOCITIES_D[14][2];
		stream_field_d[Q_LBM*idx+14]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+14];
		nx=x-LATTICE_VELOCITIES_D[15][0];
		ny=y-LATTICE_VELOCITIES_D[15][1];
		nz=z-LATTICE_VELOCITIES_D[15][2];
		stream_field_d[Q_LBM*idx+15]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+15];
		nx=x-LATTICE_VELOCITIES_D[16][0];
		ny=y-LATTICE_VELOCITIES_D[16][1];
		nz=z-LATTICE_VELOCITIES_D[16][2];
		stream_field_d[Q_LBM*idx+16]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+16];
		nx=x-LATTICE_VELOCITIES_D[17][0];
		ny=y-LATTICE_VELOCITIES_D[17][1];
		nz=z-LATTICE_VELOCITIES_D[17][2];
		stream_field_d[Q_LBM*idx+17]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+17];
		nx=x-LATTICE_VELOCITIES_D[18][0];
		ny=y-LATTICE_VELOCITIES_D[18][1];
		nz=z-LATTICE_VELOCITIES_D[18][2];
		stream_field_d[Q_LBM*idx+18]=collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+18];
	}
}


//__device__ void DoCollision(int *flag_field_d, int x, int y, int z){
__global__ void DoCollision(int *flag_field_d){
	int x = threadIdx.x+blockIdx.x*blockDim.x;
	int y = threadIdx.y+blockIdx.y*blockDim.y;
	int z = threadIdx.z+blockIdx.z*blockDim.z;
	//	__shared__ float collide_field_s[BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE*Q_LBM];
	//	int idx_block = threadIdx.x+threadIdx.y*blockDim.x+threadIdx.z*blockDim.x*blockDim.y;
	int step=xlength_d+2, idx=x+y*step+z*step*step;
	float density, velocity[D_LBM], feq[Q_LBM], *current_cell_s;

	//check that indices are within the bounds since there could be more threads than needed
	if (0<x && x<(step-1) && 0<y && y<(step-1) && 0<z && z<(step-1)){
		//copy current cell values into shared memory
//		collide_field_s[Q_LBM*idx_block]=collide_field_d[Q_LBM*idx];
//		collide_field_s[Q_LBM*idx_block+1]=collide_field_d[Q_LBM*idx+1];
//		collide_field_s[Q_LBM*idx_block+2]=collide_field_d[Q_LBM*idx+2];
//		collide_field_s[Q_LBM*idx_block+3]=collide_field_d[Q_LBM*idx+3];
//		collide_field_s[Q_LBM*idx_block+4]=collide_field_d[Q_LBM*idx+4];
//		collide_field_s[Q_LBM*idx_block+5]=collide_field_d[Q_LBM*idx+5];
//		collide_field_s[Q_LBM*idx_block+6]=collide_field_d[Q_LBM*idx+6];
//		collide_field_s[Q_LBM*idx_block+7]=collide_field_d[Q_LBM*idx+7];
//		collide_field_s[Q_LBM*idx_block+8]=collide_field_d[Q_LBM*idx+8];
//		collide_field_s[Q_LBM*idx_block+9]=collide_field_d[Q_LBM*idx+9];
//		collide_field_s[Q_LBM*idx_block+10]=collide_field_d[Q_LBM*idx+10];
//		collide_field_s[Q_LBM*idx_block+11]=collide_field_d[Q_LBM*idx+11];
//		collide_field_s[Q_LBM*idx_block+12]=collide_field_d[Q_LBM*idx+12];
//		collide_field_s[Q_LBM*idx_block+13]=collide_field_d[Q_LBM*idx+13];
//		collide_field_s[Q_LBM*idx_block+14]=collide_field_d[Q_LBM*idx+14];
//		collide_field_s[Q_LBM*idx_block+15]=collide_field_d[Q_LBM*idx+15];
//		collide_field_s[Q_LBM*idx_block+16]=collide_field_d[Q_LBM*idx+16];
//		collide_field_s[Q_LBM*idx_block+17]=collide_field_d[Q_LBM*idx+17];
//		collide_field_s[Q_LBM*idx_block+18]=collide_field_d[Q_LBM*idx+18];

//		current_cell_s = &collide_field_s[Q_LBM*idx_block];
		current_cell_s=&collide_field_d[Q_LBM*idx];
		ComputeDensityGpu(current_cell_s,&density);
		ComputeVelocityGpu(current_cell_s,&density,velocity);
		ComputeFeqGpu(&density,velocity,feq);
		ComputePostCollisionDistributionsGpu(current_cell_s,feq);

		//copy data back
//		collide_field_d[Q_LBM*idx]=collide_field_s[Q_LBM*idx_block];
//		collide_field_d[Q_LBM*idx+1]=collide_field_s[Q_LBM*idx_block+1];
//		collide_field_d[Q_LBM*idx+2]=collide_field_s[Q_LBM*idx_block+2];
//		collide_field_d[Q_LBM*idx+3]=collide_field_s[Q_LBM*idx_block+3];
//		collide_field_d[Q_LBM*idx+4]=collide_field_s[Q_LBM*idx_block+4];
//		collide_field_d[Q_LBM*idx+5]=collide_field_s[Q_LBM*idx_block+5];
//		collide_field_d[Q_LBM*idx+6]=collide_field_s[Q_LBM*idx_block+6];
//		collide_field_d[Q_LBM*idx+7]=collide_field_s[Q_LBM*idx_block+7];
//		collide_field_d[Q_LBM*idx+8]=collide_field_s[Q_LBM*idx_block+8];
//		collide_field_d[Q_LBM*idx+9]=collide_field_s[Q_LBM*idx_block+9];
//		collide_field_d[Q_LBM*idx+10]=collide_field_s[Q_LBM*idx_block+10];
//		collide_field_d[Q_LBM*idx+11]=collide_field_s[Q_LBM*idx_block+11];
//		collide_field_d[Q_LBM*idx+12]=collide_field_s[Q_LBM*idx_block+12];
//		collide_field_d[Q_LBM*idx+13]=collide_field_s[Q_LBM*idx_block+13];
//		collide_field_d[Q_LBM*idx+14]=collide_field_s[Q_LBM*idx_block+14];
//		collide_field_d[Q_LBM*idx+15]=collide_field_s[Q_LBM*idx_block+15];
//		collide_field_d[Q_LBM*idx+16]=collide_field_s[Q_LBM*idx_block+16];
//		collide_field_d[Q_LBM*idx+17]=collide_field_s[Q_LBM*idx_block+17];
//		collide_field_d[Q_LBM*idx+18]=collide_field_s[Q_LBM*idx_block+18];
	}
}


//__device__ void TreatBoundary(int *flag_field_d, int x, int y, int z){
__global__ void TreatBoundary(int *flag_field_d){
	int x = threadIdx.x+blockIdx.x*blockDim.x;
	int y = threadIdx.y+blockIdx.y*blockDim.y;
	int z = threadIdx.z+blockIdx.z*blockDim.z;

	int step=xlength_d+2, idx=x+y*step+z*step*step, nx, ny, nz, i, boundary_side=0, boundary_idx=100500;
	float density, dot_prod;

	if(idx<num_cells_d) {
		if(flag_field_d[idx] == BOTTOM_BOUNDARY) {
			boundary_side = BOTTOM_BOUNDARY;
			boundary_idx = BOTTOM_BOUNDARY_IDX;
		} else if (flag_field_d[idx] == LEFT_BOUNDARY) {
			boundary_side = LEFT_BOUNDARY;
			boundary_idx = LEFT_BOUNDARY_IDX;
		} else if (flag_field_d[idx] == RIGHT_BOUNDARY) {
			boundary_side = RIGHT_BOUNDARY;
			boundary_idx = RIGHT_BOUNDARY_IDX;
		} else if (flag_field_d[idx] == BACK_BOUNDARY) {
			boundary_side = BACK_BOUNDARY;
			boundary_idx = BACK_BOUNDARY_IDX;
		} else if (flag_field_d[idx] == FRONT_BOUNDARY) {
			boundary_side = FRONT_BOUNDARY;
			boundary_idx = FRONT_BOUNDARY_IDX;
		} else if (flag_field_d[idx] == LEFT_BOTTOM_EDGE) {
			boundary_side = LEFT_BOTTOM_EDGE;
			boundary_idx = 13;
		} else if (flag_field_d[idx] == RIGHT_BOTTOM_EDGE) {
			boundary_side = RIGHT_BOTTOM_EDGE;
			boundary_idx = 11;
		} else if (flag_field_d[idx] == BACK_BOTTOM_EDGE) {
			boundary_side = BACK_BOTTOM_EDGE;
			boundary_idx = 18;
		} else if (flag_field_d[idx] == FRONT_BOTTOM_EDGE) {
			boundary_side = FRONT_BOTTOM_EDGE;
			boundary_idx = 4;
		} else if (flag_field_d[idx] == LEFT_BACK_EDGE) {
			boundary_side = LEFT_BACK_EDGE;
			boundary_idx = 17;
		} else if (flag_field_d[idx] == LEFT_FRONT_EDGE) {
			boundary_side = LEFT_FRONT_EDGE;
			boundary_idx = 3;
		} else if (flag_field_d[idx] == RIGHT_BACK_EDGE) {
			boundary_side = RIGHT_BACK_EDGE;
			boundary_idx = 15;
		} else if (flag_field_d[idx] == RIGHT_FRONT_EDGE) {
			boundary_side = RIGHT_FRONT_EDGE;
			boundary_idx = 1;
		} else if (flag_field_d[idx] == LEFT_UPPER_EDGE) {
			boundary_side = LEFT_UPPER_EDGE;
			boundary_idx = 7;
		} else if (flag_field_d[idx] == RIGHT_UPPER_EDGE) {
			boundary_side = RIGHT_UPPER_EDGE;
			boundary_idx = 5;
		} else if (flag_field_d[idx] == BACK_UPPER_EDGE) {
			boundary_side = BACK_UPPER_EDGE;
			boundary_idx = 14;
		} else if (flag_field_d[idx] == FRONT_UPPER_EDGE) {
			boundary_side = FRONT_UPPER_EDGE;
			boundary_idx = 0;
		} else if (flag_field_d[idx] == TOP_BOUNDARY) {
			boundary_side = TOP_BOUNDARY;
			boundary_idx = TOP_BOUNDARY_IDX;
		}

		if( boundary_side==LEFT_BOUNDARY || boundary_side==RIGHT_BOUNDARY ||
				boundary_side==BOTTOM_BOUNDARY ||
				boundary_side==BACK_BOUNDARY || boundary_side==FRONT_BOUNDARY) {
			i = treat_boundary_indeces[boundary_idx][0];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			collide_field_d[Q_LBM*(x+y*step+z*step*step)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)];
			i = treat_boundary_indeces[boundary_idx][1];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			collide_field_d[Q_LBM*(x+y*step+z*step*step)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)];
			i = treat_boundary_indeces[boundary_idx][2];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			collide_field_d[Q_LBM*(x+y*step+z*step*step)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)];
			i = treat_boundary_indeces[boundary_idx][3];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			collide_field_d[Q_LBM*(x+y*step+z*step*step)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)];
			i = treat_boundary_indeces[boundary_idx][4];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			collide_field_d[Q_LBM*(x+y*step+z*step*step)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)];
		} else if (boundary_side == LEFT_BOTTOM_EDGE || boundary_side == RIGHT_BOTTOM_EDGE ||
				boundary_side == BACK_BOTTOM_EDGE || boundary_side == FRONT_BOTTOM_EDGE ||
				boundary_side == LEFT_BACK_EDGE || boundary_side == LEFT_FRONT_EDGE ||
				boundary_side == RIGHT_BACK_EDGE || boundary_side == RIGHT_FRONT_EDGE) {
			nx=x+LATTICE_VELOCITIES_D[boundary_idx][0];
			ny=y+LATTICE_VELOCITIES_D[boundary_idx][1];
			nz=z+LATTICE_VELOCITIES_D[boundary_idx][2];
			collide_field_d[Q_LBM*(x+y*step+z*step*step)+boundary_idx]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(boundary_idx)];
		} else if(boundary_side == LEFT_UPPER_EDGE || boundary_side == RIGHT_UPPER_EDGE ||
				boundary_side == BACK_UPPER_EDGE || boundary_side == FRONT_UPPER_EDGE) {
			i = boundary_idx;
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			/* Compute density in the neighbour cell */
			ComputeDensityGpu(&collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)],&density);
			/* Compute dot product */
			dot_prod=LATTICE_VELOCITIES_D[i][0]*wall_velocity_d[0]+
					LATTICE_VELOCITIES_D[i][1]*wall_velocity_d[1]+
					LATTICE_VELOCITIES_D[i][2]*wall_velocity_d[2];
			/* Assign the boudary cell value */
			collide_field_d[Q_LBM*(idx)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)]+
					2*LATTICE_WEIGHTS_D[i]*density*C_S_POW2_INV*dot_prod;
		} else if(boundary_side == TOP_BOUNDARY) {
			i = treat_boundary_indeces[boundary_idx][0];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			ComputeDensityGpu(&collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)],&density);
			dot_prod=LATTICE_VELOCITIES_D[i][0]*wall_velocity_d[0]+
					LATTICE_VELOCITIES_D[i][1]*wall_velocity_d[1]+
					LATTICE_VELOCITIES_D[i][2]*wall_velocity_d[2];
			collide_field_d[Q_LBM*(idx)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)]+
					2*LATTICE_WEIGHTS_D[i]*density*C_S_POW2_INV*dot_prod;
			i = treat_boundary_indeces[boundary_idx][1];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			ComputeDensityGpu(&collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)],&density);
			dot_prod=LATTICE_VELOCITIES_D[i][0]*wall_velocity_d[0]+
					LATTICE_VELOCITIES_D[i][1]*wall_velocity_d[1]+
					LATTICE_VELOCITIES_D[i][2]*wall_velocity_d[2];
			collide_field_d[Q_LBM*(idx)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)]+
					2*LATTICE_WEIGHTS_D[i]*density*C_S_POW2_INV*dot_prod;
			i = treat_boundary_indeces[boundary_idx][2];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			ComputeDensityGpu(&collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)],&density);
			dot_prod=LATTICE_VELOCITIES_D[i][0]*wall_velocity_d[0]+
					LATTICE_VELOCITIES_D[i][1]*wall_velocity_d[1]+
					LATTICE_VELOCITIES_D[i][2]*wall_velocity_d[2];
			collide_field_d[Q_LBM*(idx)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)]+
					2*LATTICE_WEIGHTS_D[i]*density*C_S_POW2_INV*dot_prod;
			i = treat_boundary_indeces[boundary_idx][3];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			ComputeDensityGpu(&collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)],&density);
			dot_prod=LATTICE_VELOCITIES_D[i][0]*wall_velocity_d[0]+
					LATTICE_VELOCITIES_D[i][1]*wall_velocity_d[1]+
					LATTICE_VELOCITIES_D[i][2]*wall_velocity_d[2];
			collide_field_d[Q_LBM*(idx)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)]+
					2*LATTICE_WEIGHTS_D[i]*density*C_S_POW2_INV*dot_prod;
			i = treat_boundary_indeces[boundary_idx][4];
			nx=x+LATTICE_VELOCITIES_D[i][0];
			ny=y+LATTICE_VELOCITIES_D[i][1];
			nz=z+LATTICE_VELOCITIES_D[i][2];
			ComputeDensityGpu(&collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)],&density);
			dot_prod=LATTICE_VELOCITIES_D[i][0]*wall_velocity_d[0]+
					LATTICE_VELOCITIES_D[i][1]*wall_velocity_d[1]+
					LATTICE_VELOCITIES_D[i][2]*wall_velocity_d[2];
			collide_field_d[Q_LBM*(idx)+i]=
					collide_field_d[Q_LBM*(nx+ny*step+nz*step*step)+inv2(i)]+
					2*LATTICE_WEIGHTS_D[i]*density*C_S_POW2_INV*dot_prod;
		}
	}
}

__global__ void DoSwap(){
	int x = threadIdx.x+blockIdx.x*blockDim.x;
	int y = threadIdx.y+blockIdx.y*blockDim.y;
	int z = threadIdx.z+blockIdx.z*blockDim.z;
    int step=xlength_d+2, idx=x+y*step+z*step*step;
    float *swap = NULL;
	if(idx==0) {
		swap=collide_field_d; collide_field_d=stream_field_d; stream_field_d=swap;
	}
}

/**
 * Performs streaming, collision and treatment step.
 */
__global__ void StreamCollideTreat(int *flag_field_d){
	int x = threadIdx.x+blockIdx.x*blockDim.x;
	int y = threadIdx.y+blockIdx.y*blockDim.y;
	int z = threadIdx.z+blockIdx.z*blockDim.z;

	DoStreaming(flag_field_d, x, y, z);
	__syncthreads();
//	DoSwap(x,y,z);
//	__syncthreads();
//	DoCollision(flag_field_d, x, y, z);
//	__syncthreads();
//	TreatBoundary(flag_field_d, x, y, z);

}


void DoIteration(float *collide_field, float *stream_field, int *flag_field, float tau,
		float *wall_velocity, int xlength, float **collide_field_dd, float **stream_field_dd,
		int **flag_field_d, float *mlups_sum){
	int num_cells = pow(xlength+2, D_LBM);
	float *swap=NULL;
	size_t computational_field_size = Q_LBM*num_cells*sizeof(float);
	clock_t mlups_time;

	//initialize constant data
	cudaErrorCheck(cudaMemcpyToSymbol(xlength_d, &xlength, sizeof(int), 0, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemcpyToSymbol(num_cells_d, &num_cells, sizeof(int), 0, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemcpyToSymbol(tau_d, &tau, sizeof(float), 0, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemcpyToSymbol(wall_velocity_d, wall_velocity, D_LBM*sizeof(float), 0, cudaMemcpyHostToDevice));

	cudaErrorCheck(cudaMemcpyToSymbol(collide_field_d, collide_field_dd, sizeof(*collide_field_dd), 0, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemcpyToSymbol(stream_field_d, stream_field_dd, sizeof(*stream_field_dd), 0, cudaMemcpyHostToDevice));

	//define grid structure
	//NOTE:redundant threads for boundary cells are not accounted for
	dim3 block(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
	dim3 grid((xlength+2+block.x-1)/block.x, (xlength+2+block.y-1)/block.y, (xlength+2+block.z-1)/block.z);

	mlups_time = clock();

	//perform streaming
//	DoStreaming<<<grid,block>>>(*flag_field_d);
//	cudaErrorCheck(cudaPeekAtLastError());
//	cudaErrorCheck(cudaThreadSynchronize());

	StreamCollideTreat<<<grid,block>>>(*flag_field_d);
	cudaErrorCheck(cudaPeekAtLastError());
	cudaErrorCheck(cudaThreadSynchronize());

	/* Perform the swapping of collide and stream fields */
//	swap=*collide_field_dd; *collide_field_dd=*stream_field_dd; *stream_field_dd=swap;
	DoSwap<<<grid,block>>>();
	cudaErrorCheck(cudaThreadSynchronize());

	//perform collision
	DoCollision<<<grid,block>>>(*flag_field_d);
	cudaErrorCheck(cudaPeekAtLastError());
	cudaErrorCheck(cudaThreadSynchronize());

	TreatBoundary<<<grid,block>>>(*flag_field_d);
	cudaErrorCheck(cudaPeekAtLastError());
	cudaErrorCheck(cudaThreadSynchronize());

	mlups_time = clock()-mlups_time;

	*mlups_sum += num_cells/(MLUPS_EXPONENT*(float)mlups_time/CLOCKS_PER_SEC);
	if(VERBOSE)
		printf("MLUPS: %f\n", num_cells/(MLUPS_EXPONENT*(float)mlups_time/CLOCKS_PER_SEC));

	//TODO:do we need to copy it every time?
	//copy data back to host
	cudaErrorCheck(cudaMemcpyFromSymbol(collide_field_dd, collide_field_d, sizeof(*collide_field_dd), 0, cudaMemcpyDeviceToHost));
	cudaErrorCheck(cudaMemcpyFromSymbol(stream_field_dd, stream_field_d, sizeof(*stream_field_dd), 0, cudaMemcpyDeviceToHost));

	cudaErrorCheck(cudaMemcpy(collide_field, *collide_field_dd, computational_field_size, cudaMemcpyDeviceToHost));
	cudaErrorCheck(cudaMemcpy(stream_field, *stream_field_dd, computational_field_size, cudaMemcpyDeviceToHost));
}
