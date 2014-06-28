#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_


/** Reads the parameters for the lid driven cavity scenario from a configuration file */
void ReadParameters(
    int *xlength,                      /* reads domain size */
    float *tau,                        /* relaxation parameter tau */
    float *velocity_wall,              /* velocity of the lid */
    int *timesteps,                    /* number of time steps */
    int *timesteps_per_plotting,       /* time steps between subsequent VTK plots */
    int argc,                          /* number of arguments. Should equal 2 (program + name of config file */
    char *argv[],                      /* argv[1] shall contain the path to the config file */
    int *gpu_enabled,
    int *gpu_streaming,
    int *gpu_collision,
    int *gpu_boundaries
);


/* Initializes the particle distribution functions and the flag field */
void InitialiseFields(float *collide_field, float *stream_field,int *flag_field, int xlength);

#endif
