#include "initialization_gpu.h"
#include "lbm_model.h"
#include "utils_gpu.h"

void InitialiseDeviceFields(float *collide_field, float *stream_field,int *flag_field, int xlength, float **collide_field_d, float **stream_field_d,int **flag_field_d){
	int num_cells = pow(xlength+2, D_LBM);
	size_t computational_field_size = Q_LBM*num_cells*sizeof(float);
	size_t flag_field_size = num_cells*sizeof(int);

	cudaErrorCheck(cudaMalloc(collide_field_d, computational_field_size));
	cudaErrorCheck(cudaMemcpy(*collide_field_d, collide_field, computational_field_size, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMalloc(stream_field_d, computational_field_size));
	cudaErrorCheck(cudaMemcpy(*stream_field_d, stream_field, computational_field_size, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMalloc(flag_field_d, flag_field_size));
	cudaErrorCheck(cudaMemcpy(*flag_field_d, flag_field, flag_field_size, cudaMemcpyHostToDevice));
}


void FreeDeviceFields(float **collide_field_d, float **stream_field_d,int **flag_field_d){
	cudaErrorCheck(cudaFree(*collide_field_d));
	cudaErrorCheck(cudaFree(*stream_field_d));
	cudaErrorCheck(cudaFree(*flag_field_d));
}
