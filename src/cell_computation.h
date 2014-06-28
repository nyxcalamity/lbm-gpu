#ifndef _CELL_COMPUTATIONS_H_
#define _CELL_COMPUTATIONS_H_

/**
 * Computes the density from the particle distribution functions stored at currentCell.
 * currentCell thus denotes the address of the first particle distribution function of the
 * respective cell. The result is stored in density.
 */
void ComputeDensity(const float *const current_cell, float *density);

/**
 * Computes the velocity within currentCell and stores the result in velocity
 */
void ComputeVelocity(const float * const current_cell, const float * const density, float *velocity);

/**
 * Computes the equilibrium distributions for all particle distribution functions of one
 * cell from density and velocity and stores the results in feq.
 */
void ComputeFeq(const float * const density, const float * const velocity, float *feq);

#endif
