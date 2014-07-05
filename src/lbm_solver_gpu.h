#ifndef _COLLISION_GPU_H_
#define _COLLISION_GPU_H_

/**
 * Performs streaming, collision and boundary treatment.
 */
void DoIteration(float *collide_field, float *stream_field, int *flag_field, float tau,
		float *wall_velocity, int xlength, float **collide_field_d, float **stream_field_d,
		int **flag_field_d, float *mlups_sum);

#endif
