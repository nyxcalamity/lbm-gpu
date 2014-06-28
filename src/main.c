#ifndef _MAIN_C_
#define _MAIN_C_

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "lbm_model.h"
#include "initialization.h"
#include "streaming.h"
#include "streaming_gpu.h"
#include "collision.h"
#include "collision_gpu.h"
#include "boundary.h"
#include "boundary_gpu.h"
#include "visualization.h"
#include "utils.h"

int main(int argc, char *argv[]) {
	float *collide_field=NULL, *stream_field=NULL, *swap=NULL, tau, wall_velocity[D_LBM], num_cells,
			mlups_sum;
	int *flag_field=NULL, xlength, t, timesteps, timesteps_per_plotting,
			gpu_enabled, gpu_streaming, gpu_collision, gpu_boundaries;
	clock_t mlups_time;
	size_t field_size;

	//process parameters
	ReadParameters(&xlength, &tau, wall_velocity, &timesteps, &timesteps_per_plotting, argc, argv,
			&gpu_enabled, &gpu_streaming, &gpu_collision, &gpu_boundaries);

	//check if provided parameters are legitimate
	ValidateModel(wall_velocity, xlength, tau);

	//pre-computing constants
	num_cells = pow(xlength + 2, D_LBM);
	field_size = Q_LBM*num_cells*sizeof(float);

	//initializing fields
	collide_field = (float*) malloc(field_size);
	stream_field = (float*) malloc(field_size);
	flag_field = (int*) malloc(num_cells*sizeof(int));
	InitialiseFields(collide_field, stream_field, flag_field, xlength);

	for (t = 0; t < timesteps; t++) {
		mlups_time = clock();
		/* Copy pdfs from neighbouring cells into collide field */
		if (gpu_enabled || gpu_streaming)
			DoStreamingGpu(collide_field, stream_field, flag_field, xlength);
		else
			DoStreaming(collide_field, stream_field, flag_field, xlength);
		/* Perform the swapping of collide and stream fields */
		swap = collide_field; collide_field = stream_field; stream_field = swap;
		/* Compute post collision distributions */
		if (gpu_enabled || gpu_collision)
			DoCollisionGpu(collide_field, flag_field, tau, xlength);
		else
			DoCollision(collide_field, flag_field, tau, xlength);
		/* Treat boundaries */
		if (gpu_enabled || gpu_boundaries)
			TreatBoundaryGpu(collide_field, flag_field, wall_velocity, xlength);
		else
			TreatBoundary(collide_field, flag_field, wall_velocity, xlength);
		/* Print out the MLUPS value */
		mlups_time = clock()-mlups_time;
		printf("Time step: #%d\n", t);
		mlups_sum += num_cells/(MLUPS_EXPONENT*(float)mlups_time/CLOCKS_PER_SEC);
		/* Print out vtk output if needed */
		if (t % timesteps_per_plotting == 0)
			WriteVtkOutput(collide_field, flag_field, "img/lbm-img", t, xlength);
	}

	printf("Average MLUPS: %f\n", mlups_sum/(t+1));

	if (VERBOSE)
		WriteField(collide_field, "img/collide-field", t, xlength, gpu_enabled);

	/* Free memory */
	free(collide_field);
	free(stream_field);
	free(flag_field);

	printf("Simulation complete.\n");
	return 0;
}
#endif
