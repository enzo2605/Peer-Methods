/**
 * Title: Peer methods
 * Author: Vincenzo Iannucci
 * **/
#include <math.h>
#include <cblas.h>
#include <string.h>
#include "CLab.h"
#include "utilities.h"

#define a  1.5
#define B1 0.45
#define B2 0.3611 
#define F  0.802 
#define H  0.802 
#define S  0.0002 
#define d  500 
#define D  0.802

int main(int argc, char *argv[]) {
    // intervals
    double t_start, t_end, x_start, x_end;

    // Checking arguments passed by the user
    if (argc != 5) {
        printf("Usage: %s t_start t_end x_start x_end\n", argv[0]);
        printf("Using default parameters...\n\n");
        t_start = 0.0f;
        t_end   = 50.0f;
        x_start = -50.0f;
        x_end   = 50.0f;
    }
    else {
        t_start = atof(argv[1]);
        t_end   = atof(argv[2]);
        x_start = atof(argv[3]);
        x_end   = atof(argv[4]);
    }

    // Printing the interval used
    printf("Interval used:\ntime span: [%f, %f]\nspace span: [%f, %f]\n", t_start, t_end, x_start, x_end);

    // Time initialization
    double t_span[2] = { t_start, t_end };
    double Delta_t = 1.0f / pow(2.0f, 11.0f);

    double *t_int;
    int n_points_t;
    t_int = setDVector(t_int, t_start, t_end, Delta_t, &n_points_t);
    int N = (t_span[1] - t_span[0]) / Delta_t;

    // Space initialization
    int M = 64;
    double x_span[2] = { x_start, x_end };
    double Delta_x = (x_span[1] - x_span[0]) / M;

    double *x_int;
    int n_points_x;
    x_int = setDVector(x_int, x_start, x_end, Delta_x, &n_points_x);

    /*** Define initial conditions ***/
    // u10_time
    // Allocation of the vector
    double *tempVecU10 = (double *)Calloc(M, sizeof(double));
    // Initialization with random values between 0 and 1
    initializeRandomVector(tempVecU10, M);
    // vector by scalar product
    cblas_dscal(M, 0.7, tempVecU10, 1);
    // Add the scalar alpha at every element of the array tempVecU10
    double *u10_time = sumScalarByVector(tempVecU10, M, 1.4f);
    printDVector(u10_time, M);

    // u20_time
    // Allocation of the vector
    double *tempVecU20 = (double *)Calloc(M, sizeof(double));
    // Initialization with random values between 0 and 1
    initializeRandomVector(tempVecU20, M);
    // vector by scalar product
    cblas_dscal(M, 0.7, tempVecU20, 1);
    // Add the scalar alpha at every element of the array tempVecU20
    double *u20_time = sumScalarByVector(tempVecU20, M, 1.4f);
    printDVector(u20_time, M);

    // w0_time
    // Allocation of the vector
    double *tempVecW0 = (double *)Calloc(M, sizeof(double));
    // Initialization with random values between 0 and 1
    initializeRandomVector(tempVecW0, M);
    // vector by scalar product
    cblas_dscal(M, 0.07, tempVecW0, 1);
    // Add the scalar alpha at every element of the array tempVecW0
    double *w0_time = sumScalarByVector(tempVecW0, M, 0.14f);
    printDVector(w0_time, M);

    // Finite differences of order two and periodic boundary conditions
    /*
        Ldiff = (1/Delta_x^2)*(-2*eye(M)+diag(ones(M-1,1),1)+diag(ones(M-1,1),-1));
        Ldiff(1,M) = 1/Delta_x^2; Ldiff(M,1) = 1/Delta_x^2;
        L = blkdiag(Ldiff,D*Ldiff,d*Ldiff);
    */

    exit(0);
}