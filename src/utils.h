#ifndef utils
#define utils
#include <stdlib.h>
#include <math.h>
#define true  0 == 0
double sum_vect(double *in1, double *in2, double *out, int len){
    for(int i = 0; i < len; ++i)
        out[i] = in1[i] + in2[i];
}
double rosenbrock(double *x){
    return 100*pow(x[1] -x[0]*x[0], 2) + pow(1 - x[0], 2);
}

double scalar_product(double *in1, double *in2, int len){
    double out = 0;
    for(int i = 0; i < len; ++i)
        out += in1[i]*in2[i];
    return out;
}
double matrix_vector_mult(double (*A)[2], double *v, double *out, int len){
    for(int i = 0; i < len; ++i)
        for(int j = 0; j < len; ++j)
            out[i] += A[i][j] * v[j];
}

void eval_gradient(double *x, double *out){ // checked, it's correct
    out[0] = 2*(200*pow(x[0], 3) - 200*x[0]*x[1] + x[0] - 1);
    out[1] = 200*(x[1]- pow(x[0], 2));
}
double taylor_expansion_rosenbrock(double *x, double *delta){
    double gradient[2];
    double hessian[][2] = {{-400*(x[1] - x[0]*x[0]) + 800*x[0]*x[0] + 2, -400*x[0]}, {-400*x[0], 200}};
    double image = rosenbrock(x);
    double hessian_times_delta[] = {0, 0};
    eval_gradient(x, gradient);
    matrix_vector_mult(hessian, delta, hessian_times_delta, 2);
    return image + scalar_product(gradient, delta, 2) + 0.5 *scalar_product(delta, hessian_times_delta, 2);    
}
void invert(double (*in)[2], double (*out)[2], int len){ // checked, it's correct.
    double det = in[0][0] * in[1][1] - in[0][1] * in[1][0];
    out[0][0] = in[1][1] /det;
    out[0][1] = - in[0][1] /det;
    out[1][0] = - in[1][0] / det;
    out[1][1] = in[0][0] / det;
}
void solve_for_delta(double (*A)[2], double *gradient, double *delta, double len){
    double solution[] = {0, 0};

    for(int i = 0; i < len; ++i)
        for(int j = 0; j < len; ++j)
            solution[i] += A[i][j] * gradient[j];
    
    for(int i = 0; i < len; ++i)
        delta[i] = -solution[i];
}
double rho(double *x, double *delta){
    double x_plus_delta[] = {0,0};
    sum_vect(x, delta, x_plus_delta, 2);
    return (rosenbrock(x) - rosenbrock(x_plus_delta))/(rosenbrock(x) -taylor_expansion_rosenbrock(x, delta));
}
#endif