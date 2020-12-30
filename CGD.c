#include <stdlib.h>
#include <math.h>
#include <stdio.h>
double vector_mult(double *in1, double *in2, int len){
    double out = 0;
    for(int i = 0; i < len; ++i)
        out += in1[i]*in2[i];
    return out;
}
double rosenbrock(double *x){
    return 100*pow(x[1] -x[0], 2) + pow(1 - x[0], 2);
}
void eval_gradient(double *x, double *out){
    out[0] = 2*(200*pow(x[0], 3) - 200*x[0]*x[1] + x[0] - 1);
    out[1] = 200*(x[1]- pow(x[0], 2));

}
    
double calc_alpha(double sigma, double rho, double *x, double *d, int len){
    // set alpha to one
    double alpha = 1;
    double init_value_f = rosenbrock(x);
    // evaluate gradient at point x
    double grad_f[2];
    eval_gradient(x, grad_f);
    // store current x
    double xx[2];
    for(int i = 0; i < len ; ++i)
        xx[i] = x[i];
    // now change x -> x + d
    for(int i = 0; i < len; ++i)
        x[i] = x[i] + alpha * d[i];
    // evaluate f at this new point
    double new_value_f = rosenbrock(x);
    while (new_value_f > init_value_f + sigma * alpha * vector_mult(grad_f, d, len))
    {
        alpha = alpha * rho;
        for(int i = 0; i < len; ++i)
            x[i] = xx[i] + alpha * d;
        new_value_f = rosenbrock(x);
    }
    
    return alpha;
}
double calc_beta(double *new_grad, double *d, int len){
    double beta;
    double numerator = 0;
    double denominator = 0;
    for(int i = 0; i < len; ++i)
        numerator += new_grad[i]*d[i];
    for(int i = 0; i < len; ++i)
        denominator += d[i]*d[i];
    beta = numerator/denominator;
    return beta;
}
double CGD(long unsigned int iters, double sigma, double rho, double seed_x1, double seed_x2){
    double grad_vecs[iters + 1][2];
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
    for(int k = 0; k < iters; ++k){
        alpha = calc_alpha(sigma, rho, x_vecs[k], d_vecs[k], 2);
        for(int i = 0; i < 2; ++i)
            x_vecs[k + 1][i] = x_vecs[k][i] + alpha * d_vecs[k][i];
        // we moved to the next point. Now we must compute the next direction we will jump in
        // for that we have to compute the beta_k, but to do that we need to know the -gradient (r) value at this point we just moved to
        eval_gradient(x_vecs[k], grad_vecs[k]);
        // now we compute beta
        beta = calc_beta(grad_vecs[k + 1], d_vecs[k], 2);
        for(int i = 0; i < 2; ++i)
            d_vecs[k + 1][i] = -grad_vecs[k][i] + beta * d_vecs[k][i];
        //printf("Convergence point: (%f, %f)\n", x_vecs[iters][0], x_vecs[iters][1]);

    }
    printf("Convergence point: (%f, %f)\n", x_vecs[iters][0], x_vecs[iters][1]);
    return 0.4345;
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

   unsigned long iters = 1e5;
   double sigma = 1e-4;
   double rho = 0.5;
   double x1 = -1.5;
   double x2 = -1;
   CGD(iters, sigma, rho, x1, x2);
   return 0;
}