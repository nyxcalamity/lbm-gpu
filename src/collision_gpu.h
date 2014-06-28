#ifndef _COLLISION_GPU_H_
#define _COLLISION_GPU_H_

/** Carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void DoCollisionGpu(float *collide_field, int *flag_field, float tau, int xlength);

#endif
