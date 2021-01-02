#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "utils.h"
    
void do_displacement(double sigma, double rho, double *x, double *new_x, double *d, int len){
    // set alpha to one
    double alpha = 1;
    double init_value_f = rosenbrock(x);
    // evaluate gradient at point x
    double grad_f[2];
    eval_gradient(x, grad_f);
    // now change x -> x + d
    for(int i = 0; i < len; ++i)
        new_x[i] = x[i] + alpha * d[i];
    // evaluate f at this new point
    double new_value_f = rosenbrock(new_x);
    // printf("Alpha: %f\n", alpha);
    while (new_value_f > init_value_f + sigma * alpha * scalar_product(grad_f, d, len))
    {
        alpha = alpha * rho;
        // printf("Alpha: %.20g\n", alpha);
        for(int i = 0; i < len; ++i)
            new_x[i] = x[i] + alpha * d[i];
        new_value_f = rosenbrock(new_x);
        if(alpha < 1e-4) break;
    }
}
double calc_beta(double *new_grad, double *d, int len){
    double beta;
    double numerator = 0;
    double denominator = 0;
    for(int i = 0; i < len; ++i)
        numerator += new_grad[i]*d[i];
    for(int i = 0; i < len; ++i)
        denominator += d[i]*d[i];
    // if(denominator == 0) puts("Denominator is 0");
    beta = numerator/denominator;
    return beta;
}
double CGD(long unsigned int iters, double sigma, double rho, double seed_x1, double seed_x2){
    double **grad_vecs;
    grad_vecs = (double**)malloc((iters+1) * sizeof(double));
    for(int i = 0; i < iters + 1; ++i) grad_vecs[i] = (double*)malloc(2*sizeof(double));
    double **d_vecs;
    d_vecs = (double**)malloc((iters+1) * sizeof(double));
    for(int i = 0; i < iters + 1; ++i) d_vecs[i] = (double*)malloc(2*sizeof(double));
    double **x_vecs;
    x_vecs = (double**)malloc((iters+1) * sizeof(double));
    for(int i = 0; i < iters + 1; ++i) x_vecs[i] = (double*)malloc(2*sizeof(double));
    double beta;
    double alpha;
    
    //initial conditions
    x_vecs[0][0] = seed_x1;
    x_vecs[0][1] = seed_x2;
    eval_gradient(x_vecs[0],grad_vecs[0]);
    for(int i = 0; i < 2; ++i)
        d_vecs[0][i] =  - grad_vecs[0][i]; // d0 = grad0
    //iterative process
    for(int k = 0; k < iters; ++k){
        do_displacement(sigma, rho, x_vecs[k], x_vecs[k + 1], d_vecs[k], 2);
        // we moved to the next point. Now we must compute the next direction we will jump in
        // for that we have to compute the beta_k, but to do that we need to know the -gradient (r) value at this point we just moved to
        eval_gradient(x_vecs[k + 1], grad_vecs[k + 1]);
        // now we compute beta
        beta = calc_beta(grad_vecs[k + 1], d_vecs[k], 2);
        // printf("beta: %.20f\n", beta);
        for(int i = 0; i < 2; ++i)
            d_vecs[k + 1][i] = -grad_vecs[k + 1][i] + beta * d_vecs[k][i];
        //printf("Convergence point: (%f, %f)\n", x_vecs[iters][0], x_vecs[iters][1]);

    }
    printf("Convergence point: (%.30f, %.30f)\n", x_vecs[iters][0], x_vecs[iters][1]);
    for(int i  = 0; i < iters + 1; ++i){
        free(grad_vecs[i]);
        free(d_vecs[i]);
        free(x_vecs[i]);
    }
    free(grad_vecs); 
    free(d_vecs); 
    free(x_vecs);
    return 0.4345;
}

int main(int argc, char const *argv[])
{
    unsigned long iters = atoi(argv[1]);
    double sigma = 1e-4;
    double rho = 0.5;
    double x1 = -1.5;
    double x2 = -1;
    CGD(iters, sigma, rho, x1, x2);
    return 0;
}