#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "utils.h"
#define true 0 == 0

void LM(double *seed, double mu, double alpha, double beta, unsigned long iters){
    double x_vecs[iters + 1][2];
    double lambda = mu;
    for(int i = 0; i < 2; ++i)
        x_vecs[0][i] = seed[i];
    for(int k  = 0; k < iters; ++k){
        //Now we are going to solve the linear system  "" (H - mu * I) delta = -grad "" to get delta in order to take the next step
        //First we need to compute the inverse of the lhs matrix:  
        // So we define the lhs matrix is 
        
        // Now we invert it
        
        //Done! Now we are ready to solve the system to get delta
        double gradient[2];
        eval_gradient(x_vecs[k], gradient);
        double delta[2];
    
        double Hessian_plus_muI[2][2] = { //checked, it's correct
                                    {-400*(x_vecs[k][1] - x_vecs[k][0]*x_vecs[k][0]) + 800*x_vecs[k][0]*x_vecs[k][0] + 2 + lambda, -400*x_vecs[k][0]}, 
                                    {-400*x_vecs[k][0], 200 + lambda}
                                    };
        double inverted_HmuI[2][2];
        invert(Hessian_plus_muI, inverted_HmuI, 2);
        solve_for_delta(inverted_HmuI, gradient, delta, 2);
        
        if(rho(x_vecs[k], delta) <= 0){
            lambda *= alpha;
        }else{
            lambda *= beta; 
        }
        for(int i = 0; i < 2; ++i)
            x_vecs[k + 1][i] = x_vecs[k][i] +  delta[i];
    }
    printf("Convergence point: (%f, %f)\n", x_vecs[iters][0], x_vecs[iters][1]);
}
int main(int argc, char const *argv[])
{
    double alpha = 2;
    double beta = 0.5;
    double mu  = 0.005;
    unsigned long iters = atoi(argv[1]);
    double seed[] = {-1.5, -1};
    LM(seed, mu, alpha, beta, iters);
    
    return 0;
}
