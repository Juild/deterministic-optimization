#include <stdlib.h>
#include <math.h>
#include <stdio.h>
static double  A[][2] = {{202, -200}, {-200, 202}};
static double b[] = {2,0};

void calc_r(double x1, double x2, double *r, int len){
    double x[] = {x1, x2};
    double Ax[] = {0, 0};
    for(int i; i < len; ++i)
        for(int j; j < len; ++j)
            Ax[i] = A[i][j]*x[j];
    for(int i  = 0; i < len; ++i)
        r[i] =  - Ax[i] + b[i];
}
double rosenbrock(double x1, double x2){
    return 100*pow(x2 -x1, 2) + pow(1 - x1, 2);
}
double calc_alpha(double *r, double *d, int len){
    double numerator = 0;
    double denominator = 0;
    for(int i = 0; i < len; ++i)
        numerator += r[i]*d[i];
    for(int i = 0; i < len ; ++i)
        for(int j  = 0; j < len; ++j)
            denominator += d[i]*A[i][j]*d[j];
    double alpha = numerator/denominator;
    return alpha;
}
double calc_beta(double *old_r, double *new_r, int len){
    double beta;
    double numerator = 0;
    double denominator = 0;
    for(int i = 0; i < len; ++i)
        numerator += new_r[i]*new_r[i];
    for(int i = 0; i < len; ++i)
        denominator += old_r[i]*old_r[i];
    beta = numerator/denominator;
    return beta;
}
double CGD(long unsigned int iters, double seed_x1, double seed_x2){
    double r_vecs[iters + 1][2];
    double d_vecs[iters + 1][2];
    double x_vecs[iters + 1][2];
    double beta;
    double alpha;
    
    //initial conditions
    x_vecs[0][0] = seed_x1;
    x_vecs[0][1] = seed_x2;
    calc_r(seed_x1, seed_x2, r_vecs[0], 2);
    for(int i = 0; i < 2; ++i)
        d_vecs[0][i] = r_vecs[0][i]; // d0 = r0
    //iterative process
    for(int k = 0; k != iters; ++k){
        alpha = calc_alpha(r_vecs[k], d_vecs[k], 2);
        for(int i = 0; i < 2; ++i)
            x_vecs[k + 1][i] = x_vecs[k][i] + alpha * d_vecs[k][i];
        // we moved to the next point. Now we must compute the next direction we will jump in
        // for that we have to compute the beta_k, but to do that we need to know the -gradient (r) value at this point we just moved to
        calc_r(x_vecs[k][0], x_vecs[k][1], r_vecs[k + 1], 2);
        // now we compute beta
        beta = calc_beta(r[k], r[k + 1], 2);
        for(int i = 0; i < 2: ++i)
            d_vecs[k + 1][i] = r_vecs[k][i] + beta * d_vecs[k][i];

    }
    
    
}

int main(int argc, char const *argv[])
{
    /* 
    Set seed
    Compute d(k) = - Grad_f(x(k)) + beta(k)d(k-1) such that d(k)*d(k-1) = 0
    Compute alpha(k); as it is a quadratic function just use minimization of phi(alpha(k)) = f(x(k) + alpha(k)d(k)) 
    move to the next point x(k+1) = x(k) + alpha(k)d(k)
    Compute error 
    if error it's below a certain threshold, then stop the iterative process and return the last point as the one corresponding to the global minimum
    */   

    return 0;
}