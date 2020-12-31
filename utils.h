#ifndef utils
#define utils
#include <stdlib.h>
#include <math.h>

double rosenbrock(double *x){
    return 100*pow(x[1] -x[0], 2) + pow(1 - x[0], 2);
}

double scalar_product(double *in1, double *in2, int len){
    double out = 0;
    for(int i = 0; i < len; ++i)
        out += in1[i]*in2[i];
    return out;
}

void eval_gradient(double *x, double *out){
    out[0] = 2*(200*pow(x[0], 3) - 200*x[0]*x[1] + x[0] - 1);
    out[1] = 200*(x[1]- pow(x[0], 2));
}
#endif