#ifndef _BOUNDARY_GPU_H_
#define _BOUNDARY_GPU_H_


/**
 * Handles the boundaries in our simulation setup
 */
void TreatBoundaryGpu(float *collide_field, int* flag_field, float *wall_velocity, int xlength);

#endif
