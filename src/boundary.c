#include "boundary.h"
#include "lbm_model.h"
#include "cell_computation.h"

/**
 * Inverts the value of the lattice index in order to find the vector opposite to the provided one.
 * @param i
 *      index to inverse
 * @return 
 *      inversed index
 */
int inv(int i){
    return (Q_LBM-1)-i;
}


void TreatBoundary(float *collide_field, int* flag_field, float *wall_velocity, int xlength){
    int x,nx,y,ny,z,nz,i,step=xlength+2;
    float density,dot_prod;
    
    for(x=0;x<step;x++){
        for(y=0;y<step;y++){
            for(z=0;z<step;z++){
                if(flag_field[x+y*step+z*step*step]!=FLUID){
                    for(i=0;i<Q_LBM;i++){
                        nx=x+LATTICE_VELOCITIES[i][0];
                        ny=y+LATTICE_VELOCITIES[i][1];
                        nz=z+LATTICE_VELOCITIES[i][2];

                        /* We don't need the values outside of our extended domain */
                        if(0<nx && nx<step-1 && 0<ny && ny<step-1 && 0<nz && nz<step-1){
                            if (flag_field[x+y*step+z*step*step]==MOVING_WALL){
                                /* Compute density in the neighbour cell */
                                ComputeDensity(&collide_field[Q_LBM*(nx+ny*step+nz*step*step)],&density);
                                /* Compute dot product */
                                dot_prod=LATTICE_VELOCITIES[i][0]*wall_velocity[0]+
                                        LATTICE_VELOCITIES[i][1]*wall_velocity[1]+
                                        LATTICE_VELOCITIES[i][2]*wall_velocity[2];
                                /* Assign the boudary cell value */
                                collide_field[Q_LBM*(x+y*step+z*step*step)+i]=
                                        collide_field[Q_LBM*(nx+ny*step+nz*step*step)+inv(i)]+
                                        2*LATTICE_WEIGHTS[i]*density*C_S_POW2_INV*dot_prod;
                            }else if(flag_field[x+y*step+z*step*step]==NO_SLIP){
                                collide_field[Q_LBM*(x+y*step+z*step*step)+i]=
                                        collide_field[Q_LBM*(nx+ny*step+nz*step*step)+inv(i)];
                            }
                        }
                    }
                }
            }
        }
    }
}
