#ifndef _COLLISION_H_
#define _COLLISION_H_

/** carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void DoCollision(float *collide_field, int *flag_field, float tau, int xlength);

#endif
