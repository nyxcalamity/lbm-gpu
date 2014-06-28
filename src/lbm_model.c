#include <math.h>
#include <stdio.h>

#include "lbm_model.h"
#include "utils.h"


void ValidateModel(float wall_velocity[D_LBM], int domain_size, float tau){
	float u_wall_length, mach_number, reynolds_number;
	/* Compute Mach number and Reynolds number */
	u_wall_length = sqrt(wall_velocity[0]*wall_velocity[0]+
			wall_velocity[1]*wall_velocity[1]+
			wall_velocity[2]*wall_velocity[2]);
	mach_number = u_wall_length*SQRT3;
	reynolds_number = u_wall_length*domain_size*C_S_POW2_INV/(tau-0.5);
	printf("Computed Mach number: %f\n", mach_number);
	printf("Computed Reynolds number: %f\n", reynolds_number);

	/* Check if characteristic numbers are correct */
	if (mach_number >= 1)
		ERROR("Computed Mach number is too large.");
	if (reynolds_number > 500)
		ERROR("Computed Reynolds number is too large for simulation to be run on a laptop/pc.");
}
