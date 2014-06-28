#ifndef _VISUALIZATION_H_
#define _VISUALIZATION_H_

/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. You can re-use parts of the code
 *  from visual.c (VTK output for Navier-Stokes solver) and modify it for 3D datasets.
 */
void WriteVtkOutput(const float * const collide_field, const int * const flag_field,
        const char * filename, unsigned int t, int xlength);

/**
 * Function that prints out the point by point values of the provided field (4D).
 * @param field
 *          linerized 4D array, with (x,y,z,i)=Q*(x+y*(ncell+2)+z*(ncell+2)*(ncell+2))+i
 * @param ncell
 *          number of inner cells, the ones which are there before adding a boundary layer
 */
void PrintField(float *field, int ncell);

void WriteField(const float * const field, const char * filename, unsigned int t, const int xlength,
		const int rank);

#endif
