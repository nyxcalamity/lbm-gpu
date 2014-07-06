#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_


/**
 * Reads the parameters for the lid driven cavity scenario from a configuration file
 */
void ReadParameters(int *xlength, float *tau, float *velocity_wall, int *timesteps,
		int *timesteps_per_plotting, int argc, char *argv[], int *gpu_enabled);


/**
 * Initializes the particle distribution functions and the flag field
 */
void InitialiseFields(float *collide_field, float *stream_field,int *flag_field, int xlength,
		int gpu_enabled);

#endif
