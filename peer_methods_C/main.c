/**
 * Title: Peer methods
 * Author: Vincenzo Iannucci
 * **/
#include <math.h>
#include <cblas.h>
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
#define NDIM 2

int main(int argc, char *argv[]) {
    // Random initialization for the seed
    srand((unsigned int)time(NULL));

    // Intervals
    double t_start, t_end, x_start, x_end;

    // Checking arguments passed by the user
    if (argc != 5) {
        fprintf(stdout, "Usage: %s t_start t_end x_start x_end\n", argv[0]);
        fprintf(stdout, "Using default parameters...\n\n");
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
    fprintf(stdout, "Interval used:\ntime span: [%f, %f]\nspace span: [%f, %f]\n", t_start, t_end, x_start, x_end);

    /*********************************************** 
     *          Time initialization 
     * *********************************************/
    double t_span[2] = { t_start, t_end };
    double Delta_t = 1.0f / pow(2.0f, 11.0f);

    double *t_int;
    int n_points_t;
    t_int = intervalDiscretization(t_int, t_start, t_end, Delta_t, &n_points_t);
    int N = (t_span[1] - t_span[0]) / Delta_t;

    /*********************************************** 
     *          Space initialization 
     * *********************************************/
    int M = 64;
    double x_span[2] = { x_start, x_end };
    double Delta_x = (x_span[1] - x_span[0]) / M;
    printf("Delta_x: %f\n",Delta_x);

    double *x_int;
    int n_points_x;
    x_int = intervalDiscretization(x_int, x_start, x_end, Delta_x, &n_points_x);

    /*********************************************** 
     *          Define initial conditions 
     * *********************************************/
    /**** u10_time *****/
    // Allocation of the vector
    double *u10_time = (double *)Calloc(M, sizeof(double));
    // Initialization with random values between 0 and 1
    initializeRandomVector(u10_time, M);
    // vector by scalar product
    cblas_dscal(M, 0.7, u10_time, 1);
    // Add the scalar alpha at every element of the array u10_time
    sumScalarByVector(u10_time, M, 1.4f);

    /**** u20_time *****/
    // Allocation of the vector
    double *u20_time = (double *)Calloc(M, sizeof(double));
    // Initialization with random values between 0 and 1
    initializeRandomVector(u20_time, M);
    // vector by scalar product
    cblas_dscal(M, 0.7, u20_time, 1);
    // Add the scalar alpha at every element of the array u20_time
    sumScalarByVector(u20_time, M, 1.4f);

    /**** w0_time *****/
    // Allocation of the vector
    double *w0_time = (double *)Calloc(M, sizeof(double));
    // Initialization with random values between 0 and 1
    initializeRandomVector(w0_time, M);
    // vector by scalar product
    cblas_dscal(M, 0.07, w0_time, 1);
    // Add the scalar alpha at every element of the array w0_time
    sumScalarByVector(w0_time, M, 0.14f);

    /************************************************************** 
     *  Create vector y0 = [U10;U20;W0] with initial conditions 
     * ***********************************************************/
    int y0Dimension;
    double *y0 = packThreeVectors(M, u10_time, u20_time, w0_time, &y0Dimension);

    /********************************************************************
     * Finite differences of order two and periodic boundary conditions
     * *****************************************************************/
    // Calculate Ldiff
    int sizeTempDiagOne, sizeTempDiagMinusOne;
    double *eyeM = eyeD(eyeM, M);
    eyeM = scalarByMatrix(eyeM, M, M, -2.0f);
    double *onesVector = onesD(onesVector, M - 1);
    double *tempDiagOne = diagD(onesVector, M - 1, 1, &sizeTempDiagOne);
    double *tempDiagMinusOne = diagD(onesVector, M - 1, -1, &sizeTempDiagMinusOne);
    double *addend1 = sumPuntSquareMatrices(eyeM, tempDiagOne, M);
    double *Ldiff = sumPuntSquareMatrices(addend1, tempDiagMinusOne, M);
    
    Ldiff = scalarByMatrix(Ldiff, M, M, 1.0f / (Delta_x * Delta_x));
    Ldiff[M - 1] = 1.0f / (Delta_x * Delta_x);
    Ldiff[(M - 1) * M] = 1.0f / (Delta_x * Delta_x);

    // Calculate the matrix L
    int LSize;
    double *DLdiff = scalarByMatrix(Ldiff, M, M, D);
    double *dLdiff = scalarByMatrix(Ldiff, M, M, d);
    double *L = threeBlockDiagD(M, Ldiff, DLdiff, dLdiff, &LSize);
    printDMatrix(L, LSize, LSize, "L");

    // Free all the memory dynamically allocated
    freeEverything(u10_time, u20_time, w0_time, y0, eyeM, onesVector, tempDiagOne, tempDiagMinusOne, addend1, Ldiff, (void *)0);

    exit(0);
}