#ifndef _COLLISION_GPU_H_
#define _COLLISION_GPU_H_

/** Carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
//void DoCollisionGpu(float *collide_field, int *flag_field, float tau, int xlength);

/**
 * Handles the boundaries in our simulation setup
 */
//void TreatBoundaryGpu(float *collide_field, int* flag_field, float *wall_velocity, int xlength);

/**
 * Carries out the streaming step and writes the respective distribution functions from
 * collideField to streamField.
 */
//void DoStreamingGpu(float *collide_field, float *stream_field, int *flag_field, int xlength);

/**
 * Performs streaming, collision and boundary treatment.
 */
void DoIteration(float *collide_field, float *stream_field, int *flag_field, float tau,
		float *wall_velocity, int xlength);

#endif
