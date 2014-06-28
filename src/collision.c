#include "collision.h"
#include "cell_computation.h"
#include "lbm_model.h"
#include "utils.h"


/** Computes the post-collision distribution functions according to the BGK update rule and
 *  stores the results again at the same position.
 */
void ComputePostCollisionDistributions(float *current_cell, float tau, const float *const feq){
    int i;
    for(i=0;i<Q_LBM;i++){
        current_cell[i]=current_cell[i]-(current_cell[i]-feq[i])/tau;
        
        /* Probability distribution function can not be less than 0 */
        if (current_cell[i] < 0)
            ERROR("Probability distribution function can not be negative.");
    }
}


void DoCollision(float *collide_field, int *flag_field, float tau, int xlength){
    float density, velocity[3], feq[Q_LBM], *currentCell;
    int x,y,z,step=xlength+2;
    
    for(x=1;x<step-1;x++){
        for(y=1;y<step-1;y++){
            for(z=1;z<step-1;z++){
                currentCell=&collide_field[Q_LBM*(x+y*step+z*step*step)];
                
                ComputeDensity(currentCell,&density);
                ComputeVelocity(currentCell,&density,velocity);
                ComputeFeq(&density,velocity,feq);
                ComputePostCollisionDistributions(currentCell,tau,feq);
            }
        }
    }
}
