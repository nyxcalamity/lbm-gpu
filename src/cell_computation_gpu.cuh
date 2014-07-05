#ifndef _CELL_COMPUTATION_GPU_CUH_
#define _CELL_COMPUTATION_GPU_CUH_


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


// TODO: rename in inv
__device__ int inv2(int i){
    return (Q_LBM-1)-i;
}


#endif
