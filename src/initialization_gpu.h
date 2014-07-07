#ifndef _INITIALIZATION_GPU_H_
#define _INITIALIZATION_GPU_H_

/**
 * Copy everything from fields on DRAM in GPU
 */
void InitialiseDeviceFields(float *collide_field, float *stream_field,int *flag_field, int xlength,
		float **collide_field_d, float **stream_field_d,int **flag_field_d);

/**
 * Free GPU memory
 */
void FreeDeviceFields(float **collide_field_d, float **stream_field_d,int **flag_field_d);

/**
 * Copy data from device in global memory.
 */
void CopyFieldsFromDevice(float *collide_field, float *stream_field, int xlength,
		float **collide_field_dd, float **stream_field_dd);

#endif
