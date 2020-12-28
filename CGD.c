#include <stdlib.h>
#include <math.h>
#include <stdio.h>
static double  A[][2] = {{202, -200}, {-200, 202}};
static double b[] = {2,0};


double calc_r(double x1, double x2, double *r){
    double x[] = {x1, x2};
    double Ax[] = {0, 0};
    for(int i; i < 2; ++i)
        for(int j; j < 2; ++j)
            Ax[i] = A[i][j]*x[j];
    for(int i  = 0; i < 2; ++i)
        r[i] =  - Ax[i] + b[i];
}
double rosenbrock( double x1, double x2){
    return 100*pow(x2 -x1, 2) + pow(1 - x1, 2);
}

double CGD_step(){
    double r;
    
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