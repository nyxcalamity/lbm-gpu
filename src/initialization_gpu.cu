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

    int x,y,z,i,step=xlength+2;

    /* NOTE: We use z=xlength+1 as the moving wall */
    for(x=0;x<step;x++){
        for(y=0;y<step;y++){
            for(z=0;z<step;z++){

                /* Initializing distributions for stream and collide fields */
                for(i=0;i<Q_LBM;i++){
                    /* NOTE:Stream field is initilized to 0s because that helps to
                     * track down mistakes and has no impact whatsoever to on the
                     * computation further on.
                     */
                    stream_field[Q_LBM*(x+y*step+z*step*step)+i]=0;
                    collide_field[Q_LBM*(x+y*step+z*step*step)+i]=0;
                }
            }
        }
    }


}


void FreeDeviceFields(float **collide_field_d, float **stream_field_d,int **flag_field_d){
	cudaErrorCheck(cudaFree(*collide_field_d));
	cudaErrorCheck(cudaFree(*stream_field_d));
	cudaErrorCheck(cudaFree(*flag_field_d));
}


void CheckGPU(float *collide_field, float *stream_field,int *flag_field, int xlength, float *collide_field_d, float *stream_field_d,int *flag_field_d){
	int num_cells = pow(xlength+2, D_LBM);
	size_t computational_field_size = Q_LBM*num_cells*sizeof(float);
	size_t flag_field_size = num_cells*sizeof(int);

	cudaErrorCheck(cudaMemcpy(collide_field, collide_field_d, computational_field_size, cudaMemcpyDeviceToHost));
	cudaErrorCheck(cudaMemcpy(stream_field, stream_field_d, computational_field_size, cudaMemcpyDeviceToHost));
}
