#ifndef _LBM_MODEL_GPU_CUH_
#define _LBM_MODEL_GPU_CUH_

/* Restricts 3D blocks to have 512 threads (limits: 512 CC<2.x; 1024 CC>2.x) */
#define BLOCK_SIZE 8

/* Arbitrary boundary identifiers */
#define LEFT_BOUNDARY_IDX 0
#define RIGHT_BOUNDARY_IDX 1
#define BOTTOM_BOUNDARY_IDX 2
#define TOP_BOUNDARY_IDX 3
#define BACK_BOUNDARY_IDX 4
#define FRONT_BOUNDARY_IDX 5


__constant__ static const int LATTICE_VELOCITIES_D[19][3] = {
    {0,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},{0,1,-1},{-1,-1,0},{0,-1,0},{1,-1,0},
    {-1,0,0}, {0,0,0},  {1,0,0}, {-1,1,0},{0,1,0}, {1,1,0},  {0,-1,1},{-1,0,1},
    {0,0,1},  {1,0,1},  {0,1,1}
};

__constant__ static const float LATTICE_WEIGHTS_D[19] = {
    1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0,
    2.0/36.0, 12.0/36.0,2.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
    2.0/36.0, 1.0/36.0, 1.0/36.0
};

/**
 * Stores number of probability distributions which we need to copy in boundary treatment step.
 */
__constant__ static const float TREAT_BOUNDARY_INDECES[6][5] = {
	{3,7,10,13,17},
	{1,5,8,11,15},
	{4,11,12,13,18},
	{0,5,6,7,14},
	{14,15,16,17,18},
	{0,1,2,3,4}
};

#endif
