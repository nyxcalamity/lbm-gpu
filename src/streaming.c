#include "streaming.h"
#include "lbm_model.h"

void DoStreaming(float *collide_field, float *stream_field, int *flag_field, int xlength){
    int x,nx,y,ny,z,nz,i,step=xlength+2;

    for(x=0;x<step;x++){
        for(y=0;y<step;y++){
            for(z=0;z<step;z++){
                if(flag_field[x+y*step+z*step*step]==FLUID){
                    for(i=0;i<Q_LBM;i++){
                        nx=x-LATTICE_VELOCITIES[i][0];
                        ny=y-LATTICE_VELOCITIES[i][1];
                        nz=z-LATTICE_VELOCITIES[i][2];
                        
                        stream_field[Q_LBM*(x+y*step+z*step*step)+i]=
                                collide_field[Q_LBM*(nx+ny*step+nz*step*step)+i];
                    }
                }
            }
        }
    }
}
