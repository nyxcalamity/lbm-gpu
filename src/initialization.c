#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "initialization.h"
#include "lbm_model.h"
#include "utils.h"

/**
 * Prints help message that includes program usage instructions and control flags.
 */
void PrintHelpMessage(){
	printf("List of control flags:\n");
	printf("\t -gpu             all computations are to be performed on gpu\n");
	printf("\t -cpu             all computations are to be performed on cpu\n");
	printf("\t -gpu-streaming   computes everything except streaming on cpu\n");
	printf("\t -gpu-collision   computes everything except collision on cpu\n");
	printf("\t -gpu-boundaries  computes everything except boundaries on cpu\n");
	printf("\t -help            prints this help message\n");
	printf("NOTE: Control flags are mutually exclusive and only one flag at a time is allowed\n");
	printf("Example program usage:\n");
	printf("\t./lbm-sim ./data/lbm.dat -gpu\n");
	exit(1);
}


void ReadParameters(int *xlength, float *tau, float *velocity_wall, int *timesteps,
		int *timesteps_per_plotting, int argc, char *argv[], int *gpu_enabled, int *gpu_streaming,
        int *gpu_collision, int *gpu_boundaries){
    float *velocity_wall_1, *velocity_wall_2, *velocity_wall_3;
    
    if(argc<3) PrintHelpMessage();
    if(!strcmp(argv[1], "-help") || !strcmp(argv[2], "-help")) PrintHelpMessage();
    if(access(argv[1], R_OK) != 0)
        ERROR("Provided configuration file path either doesn't exist or can not be read.");
    
    if(!strcmp(argv[2], "-gpu")) *gpu_enabled=1; else *gpu_enabled=0;
    if(!strcmp(argv[2], "-gpu-streaming")) *gpu_streaming=1; else *gpu_streaming=0;
    if(!strcmp(argv[2], "-gpu-collision")) *gpu_collision=1; else *gpu_collision=0;
    if(!strcmp(argv[2], "-gpu-boundaries")) *gpu_boundaries=1; else *gpu_boundaries=0;

    READ_FLOAT(argv[1], *tau);

    velocity_wall_1=&velocity_wall[0];
    velocity_wall_2=&velocity_wall[1];
    velocity_wall_3=&velocity_wall[2];

    READ_FLOAT(argv[1], *velocity_wall_1);
    READ_FLOAT(argv[1], *velocity_wall_2);
    READ_FLOAT(argv[1], *velocity_wall_3);
    
    READ_INT(argv[1], *xlength);
    READ_INT(argv[1], *timesteps);
    READ_INT(argv[1], *timesteps_per_plotting);
}


void InitialiseFields(float *collide_field, float *stream_field, int *flag_field, int xlength){
    int x,y,z,i,step=xlength+2;
    
    /* NOTE: We use z=xlength+1 as the moving wall */
    for(x=0;x<step;x++){
        for(y=0;y<step;y++){
            for(z=0;z<step;z++){
                /* Initializing flags */
                if(z == xlength+1)
                    flag_field[x+y*step+z*step*step]=MOVING_WALL;
                else if(x == 0 || x == xlength+1 || y == 0 || y == xlength+1 || z == 0)
                    flag_field[x+y*step+z*step*step]=NO_SLIP;
                else
                    flag_field[x+y*step+z*step*step]=FLUID;
                
                /* Initializing distributions for stream and collide fields */
                for(i=0;i<Q_LBM;i++){
                    /* NOTE:Stream field is initilized to 0s because that helps to 
                     * track down mistakes and has no impact whatsoever to on the 
                     * computation further on.
                     */
                    stream_field[Q_LBM*(x+y*step+z*step*step)+i]=0;
                    collide_field[Q_LBM*(x+y*step+z*step*step)+i]=LATTICE_WEIGHTS[i];
                }
            }
        }
    }
}
