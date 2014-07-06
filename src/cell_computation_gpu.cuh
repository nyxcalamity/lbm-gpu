#ifndef _CELL_COMPUTATION_GPU_CUH_
#define _CELL_COMPUTATION_GPU_CUH_


/**
 * Computes the density from the particle distribution functions stored at currentCell.
 * currentCell thus denotes the address of the first particle distribution function of the
 * respective cell. The result is stored in density.
 */
__device__ void ComputeDensityGpu(float *current_cell, float *density){
    *density=0;
    *density+=current_cell[0];
    *density+=current_cell[1];
    *density+=current_cell[2];
    *density+=current_cell[3];
    *density+=current_cell[4];
    *density+=current_cell[5];
    *density+=current_cell[6];
    *density+=current_cell[7];
    *density+=current_cell[8];
    *density+=current_cell[9];
    *density+=current_cell[10];
    *density+=current_cell[11];
    *density+=current_cell[12];
    *density+=current_cell[13];
    *density+=current_cell[14];
    *density+=current_cell[15];
    *density+=current_cell[16];
    *density+=current_cell[17];
    *density+=current_cell[18];
}


/**
 * Computes the velocity within currentCell and stores the result in velocity
 */
__device__ void ComputeVelocityGpu(float *current_cell, float *density, float *velocity){
    velocity[0]=0;
	velocity[0]+=current_cell[0]*LATTICE_VELOCITIES_D[0][0];
	velocity[0]+=current_cell[1]*LATTICE_VELOCITIES_D[1][0];
	velocity[0]+=current_cell[2]*LATTICE_VELOCITIES_D[2][0];
	velocity[0]+=current_cell[3]*LATTICE_VELOCITIES_D[3][0];
	velocity[0]+=current_cell[4]*LATTICE_VELOCITIES_D[4][0];
	velocity[0]+=current_cell[5]*LATTICE_VELOCITIES_D[5][0];
	velocity[0]+=current_cell[6]*LATTICE_VELOCITIES_D[6][0];
	velocity[0]+=current_cell[7]*LATTICE_VELOCITIES_D[7][0];
	velocity[0]+=current_cell[8]*LATTICE_VELOCITIES_D[8][0];
	velocity[0]+=current_cell[9]*LATTICE_VELOCITIES_D[9][0];
	velocity[0]+=current_cell[10]*LATTICE_VELOCITIES_D[10][0];
	velocity[0]+=current_cell[11]*LATTICE_VELOCITIES_D[11][0];
	velocity[0]+=current_cell[12]*LATTICE_VELOCITIES_D[12][0];
	velocity[0]+=current_cell[13]*LATTICE_VELOCITIES_D[13][0];
	velocity[0]+=current_cell[14]*LATTICE_VELOCITIES_D[14][0];
	velocity[0]+=current_cell[15]*LATTICE_VELOCITIES_D[15][0];
	velocity[0]+=current_cell[16]*LATTICE_VELOCITIES_D[16][0];
	velocity[0]+=current_cell[17]*LATTICE_VELOCITIES_D[17][0];
	velocity[0]+=current_cell[18]*LATTICE_VELOCITIES_D[18][0];
	velocity[0]/=*density;

	velocity[1]=0;
	velocity[1]+=current_cell[0]*LATTICE_VELOCITIES_D[0][1];
	velocity[1]+=current_cell[1]*LATTICE_VELOCITIES_D[1][1];
	velocity[1]+=current_cell[2]*LATTICE_VELOCITIES_D[2][1];
	velocity[1]+=current_cell[3]*LATTICE_VELOCITIES_D[3][1];
	velocity[1]+=current_cell[4]*LATTICE_VELOCITIES_D[4][1];
	velocity[1]+=current_cell[5]*LATTICE_VELOCITIES_D[5][1];
	velocity[1]+=current_cell[6]*LATTICE_VELOCITIES_D[6][1];
	velocity[1]+=current_cell[7]*LATTICE_VELOCITIES_D[7][1];
	velocity[1]+=current_cell[8]*LATTICE_VELOCITIES_D[8][1];
	velocity[1]+=current_cell[9]*LATTICE_VELOCITIES_D[9][1];
	velocity[1]+=current_cell[10]*LATTICE_VELOCITIES_D[10][1];
	velocity[1]+=current_cell[11]*LATTICE_VELOCITIES_D[11][1];
	velocity[1]+=current_cell[12]*LATTICE_VELOCITIES_D[12][1];
	velocity[1]+=current_cell[13]*LATTICE_VELOCITIES_D[13][1];
	velocity[1]+=current_cell[14]*LATTICE_VELOCITIES_D[14][1];
	velocity[1]+=current_cell[15]*LATTICE_VELOCITIES_D[15][1];
	velocity[1]+=current_cell[16]*LATTICE_VELOCITIES_D[16][1];
	velocity[1]+=current_cell[17]*LATTICE_VELOCITIES_D[17][1];
	velocity[1]+=current_cell[18]*LATTICE_VELOCITIES_D[18][1];
	velocity[1]/=*density;

	velocity[2]=0;
	velocity[2]+=current_cell[0]*LATTICE_VELOCITIES_D[0][2];
	velocity[2]+=current_cell[1]*LATTICE_VELOCITIES_D[1][2];
	velocity[2]+=current_cell[2]*LATTICE_VELOCITIES_D[2][2];
	velocity[2]+=current_cell[3]*LATTICE_VELOCITIES_D[3][2];
	velocity[2]+=current_cell[4]*LATTICE_VELOCITIES_D[4][2];
	velocity[2]+=current_cell[5]*LATTICE_VELOCITIES_D[5][2];
	velocity[2]+=current_cell[6]*LATTICE_VELOCITIES_D[6][2];
	velocity[2]+=current_cell[7]*LATTICE_VELOCITIES_D[7][2];
	velocity[2]+=current_cell[8]*LATTICE_VELOCITIES_D[8][2];
	velocity[2]+=current_cell[9]*LATTICE_VELOCITIES_D[9][2];
	velocity[2]+=current_cell[10]*LATTICE_VELOCITIES_D[10][2];
	velocity[2]+=current_cell[11]*LATTICE_VELOCITIES_D[11][2];
	velocity[2]+=current_cell[12]*LATTICE_VELOCITIES_D[12][2];
	velocity[2]+=current_cell[13]*LATTICE_VELOCITIES_D[13][2];
	velocity[2]+=current_cell[14]*LATTICE_VELOCITIES_D[14][2];
	velocity[2]+=current_cell[15]*LATTICE_VELOCITIES_D[15][2];
	velocity[2]+=current_cell[16]*LATTICE_VELOCITIES_D[16][2];
	velocity[2]+=current_cell[17]*LATTICE_VELOCITIES_D[17][2];
	velocity[2]+=current_cell[18]*LATTICE_VELOCITIES_D[18][2];
    velocity[2]/=*density;
}


/**
 * Computes the equilibrium distributions for all particle distribution functions of one
 * cell from density and velocity and stores the results in feq.
 */
__device__ void ComputeFeqGpu(float *density, float *velocity, float *feq){
    float s1, s2, s3=velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2];
	s1 = LATTICE_VELOCITIES_D[0][0]*velocity[0]+LATTICE_VELOCITIES_D[0][1]*velocity[1]+LATTICE_VELOCITIES_D[0][2]*velocity[2];
	s2 = s1*s1;
	feq[0]=LATTICE_WEIGHTS_D[0]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[1][0]*velocity[0]+LATTICE_VELOCITIES_D[1][1]*velocity[1]+LATTICE_VELOCITIES_D[1][2]*velocity[2];
	s2 = s1*s1;
	feq[1]=LATTICE_WEIGHTS_D[1]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[2][0]*velocity[0]+LATTICE_VELOCITIES_D[2][1]*velocity[1]+LATTICE_VELOCITIES_D[2][2]*velocity[2];
	s2 = s1*s1;
	feq[2]=LATTICE_WEIGHTS_D[2]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[3][0]*velocity[0]+LATTICE_VELOCITIES_D[3][1]*velocity[1]+LATTICE_VELOCITIES_D[3][2]*velocity[2];
	s2 = s1*s1;
	feq[3]=LATTICE_WEIGHTS_D[3]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[4][0]*velocity[0]+LATTICE_VELOCITIES_D[4][1]*velocity[1]+LATTICE_VELOCITIES_D[4][2]*velocity[2];
	s2 = s1*s1;
	feq[4]=LATTICE_WEIGHTS_D[4]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[5][0]*velocity[0]+LATTICE_VELOCITIES_D[5][1]*velocity[1]+LATTICE_VELOCITIES_D[5][2]*velocity[2];
	s2 = s1*s1;
	feq[5]=LATTICE_WEIGHTS_D[5]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[6][0]*velocity[0]+LATTICE_VELOCITIES_D[6][1]*velocity[1]+LATTICE_VELOCITIES_D[6][2]*velocity[2];
	s2 = s1*s1;
	feq[6]=LATTICE_WEIGHTS_D[6]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[7][0]*velocity[0]+LATTICE_VELOCITIES_D[7][1]*velocity[1]+LATTICE_VELOCITIES_D[7][2]*velocity[2];
	s2 = s1*s1;
	feq[7]=LATTICE_WEIGHTS_D[7]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[8][0]*velocity[0]+LATTICE_VELOCITIES_D[8][1]*velocity[1]+LATTICE_VELOCITIES_D[8][2]*velocity[2];
	s2 = s1*s1;
	feq[8]=LATTICE_WEIGHTS_D[8]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[9][0]*velocity[0]+LATTICE_VELOCITIES_D[9][1]*velocity[1]+LATTICE_VELOCITIES_D[9][2]*velocity[2];
	s2 = s1*s1;
	feq[9]=LATTICE_WEIGHTS_D[9]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[10][0]*velocity[0]+LATTICE_VELOCITIES_D[10][1]*velocity[1]+LATTICE_VELOCITIES_D[10][2]*velocity[2];
	s2 = s1*s1;
	feq[10]=LATTICE_WEIGHTS_D[10]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[11][0]*velocity[0]+LATTICE_VELOCITIES_D[11][1]*velocity[1]+LATTICE_VELOCITIES_D[11][2]*velocity[2];
	s2 = s1*s1;
	feq[11]=LATTICE_WEIGHTS_D[11]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[12][0]*velocity[0]+LATTICE_VELOCITIES_D[12][1]*velocity[1]+LATTICE_VELOCITIES_D[12][2]*velocity[2];
	s2 = s1*s1;
	feq[12]=LATTICE_WEIGHTS_D[12]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[13][0]*velocity[0]+LATTICE_VELOCITIES_D[13][1]*velocity[1]+LATTICE_VELOCITIES_D[13][2]*velocity[2];
	s2 = s1*s1;
	feq[13]=LATTICE_WEIGHTS_D[13]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[14][0]*velocity[0]+LATTICE_VELOCITIES_D[14][1]*velocity[1]+LATTICE_VELOCITIES_D[14][2]*velocity[2];
	s2 = s1*s1;
	feq[14]=LATTICE_WEIGHTS_D[14]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[15][0]*velocity[0]+LATTICE_VELOCITIES_D[15][1]*velocity[1]+LATTICE_VELOCITIES_D[15][2]*velocity[2];
	s2 = s1*s1;
	feq[15]=LATTICE_WEIGHTS_D[15]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[16][0]*velocity[0]+LATTICE_VELOCITIES_D[16][1]*velocity[1]+LATTICE_VELOCITIES_D[16][2]*velocity[2];
	s2 = s1*s1;
	feq[16]=LATTICE_WEIGHTS_D[16]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[17][0]*velocity[0]+LATTICE_VELOCITIES_D[17][1]*velocity[1]+LATTICE_VELOCITIES_D[17][2]*velocity[2];
	s2 = s1*s1;
	feq[17]=LATTICE_WEIGHTS_D[17]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);

	s1 = LATTICE_VELOCITIES_D[18][0]*velocity[0]+LATTICE_VELOCITIES_D[18][1]*velocity[1]+LATTICE_VELOCITIES_D[18][2]*velocity[2];
	s2 = s1*s1;
	feq[18]=LATTICE_WEIGHTS_D[18]*(*density)*(1+s1*C_S_POW2_INV+s2*C_S_POW4_INV/2.0-s3*C_S_POW2_INV/2.0);
}


/**
 * Finds an inverse probability distribution for the provided lattice index.
 */
__device__ int InvGpu(int i){
    return (Q_LBM-1)-i;
}

#endif
