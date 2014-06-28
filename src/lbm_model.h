#ifndef _LBM_MODEL_H_
#define _LBM_MODEL_H_

#define VERBOSE 0
#define MLUPS_MIN_CELLS 300
#define MLUPS_EXPONENT 1000000
//simulation parameters
#define D_LBM 3
#define Q_LBM 19
//arbitrary cell identifiers
#define FLUID 0
#define NO_SLIP 1
#define MOVING_WALL 2
//computation enhancing values
#define EPS 0.05
static const float C_S_POW2_INV = 3.0;
static const float C_S_POW4_INV = 9.0;
//reference values, not used in actual computation
#define SQRT3 1.73205080756887729
static const float C_S = 1.0/SQRT3;
//predefined simulation parameters
static const int LATTICE_VELOCITIES[19][3] = {
    {0,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},{0,1,-1},{-1,-1,0},{0,-1,0},{1,-1,0},
    {-1,0,0}, {0,0,0},  {1,0,0}, {-1,1,0},{0,1,0}, {1,1,0},  {0,-1,1},{-1,0,1},
    {0,0,1},  {1,0,1},  {0,1,1}
};
static const float LATTICE_WEIGHTS[19] = {
    1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 
    2.0/36.0, 12.0/36.0,2.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 
    2.0/36.0, 1.0/36.0, 1.0/36.0
};


/**
 * Validates the configured physical model by calculating characteristic numbers
 */
void ValidateModel(float wall_velocity[D_LBM], int domain_size, float tau);

#endif
