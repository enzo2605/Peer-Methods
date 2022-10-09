/**
 * Title: Peer methods
 * Author: Vincenzo Iannucci
 * **/
#include "peerMethods.h"
#include <assert.h>

int main(int argc, char *argv[]) {
    // Random initialization for the seed
    srand((unsigned int)time(NULL));

    printf("FLT_EVAL_METHOD: %d\n", FLT_EVAL_METHOD);
    //assert(FLT_EVAL_METHOD == 2);

    // Intervals
    double_t t_start, t_end, x_start, x_end;

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
    double_t t_span[2] = { t_start, t_end };
    // for test
    double_t Delta_t = 1.0f / pow(2.0f, 3.0f);
    fprintf(stdout, "Delta_t: %f\n", Delta_t);

    double_t *t_int;
    int n_points_t;
    t_int = intervalDiscretization(t_int, t_start, t_end, Delta_t, &n_points_t);
    int N = (t_span[1] - t_span[0]) / Delta_t;
    fprintf(stdout, "N: %d\n", N);

    /*********************************************** 
     *          Space initialization 
     * *********************************************/
    double_t x_span[2] = { x_start, x_end };
    double_t Delta_x = (x_span[1] - x_span[0]) / M;
    printf("Delta_x: %f\n",Delta_x);

    double_t *x_int;
    int n_points_x;
    x_int = intervalDiscretization(x_int, x_start, x_end, Delta_x, &n_points_x);
    assert(n_points_x == M + 1);

    /*********************************************** 
     *          Define initial conditions 
     * *********************************************/

    // Test values
    double_t test_u10[] = { 0.9528, 0.7041, 0.9539, 0.5982, 0.8407, 0.4428, 0.8368, 0.5187 };
    double_t test_u20[] = { 0.0222, 0.3759, 0.8986, 0.4290, 0.1996, 0.3031, 0.5383, 0.9102 };
    double_t test_w0[]  = { 0.5253, 0.3068, 0.0345, 0.7153, 0.7687, 0.0595, 0.6271, 0.2652 };

    /**** u10_time *****/
    // Allocation of the vector
    double_t *u10_time = (double_t *)Calloc(M, sizeof(double_t));
    // Initialization with random values between 0 and 1
    initVectorWAnotherVector(u10_time, test_u10, M);
    //initializeRandomVector(u10_time, M);
    // vector by scalar product
    cblas_dscal(M, 0.7, u10_time, 1);
    // Add the scalar alpha at every element of the array u10_time
    sumScalarByVector(u10_time, M, 1.4f);
    //printDVector(u10_time, M, "u10_time");

    /**** u20_time *****/
    // Allocation of the vector
    double_t *u20_time = (double_t *)Calloc(M, sizeof(double_t));
    // Initialization with random values between 0 and 1
    initVectorWAnotherVector(u20_time, test_u20, M);
    //initializeRandomVector(u20_time, M);
    // vector by scalar product
    cblas_dscal(M, 0.7, u20_time, 1);
    // Add the scalar alpha at every element of the array u20_time
    sumScalarByVector(u20_time, M, 1.4f);
    //printDVector(u20_time, M, "u20_time");

    /**** w0_time *****/
    // Allocation of the vector
    double_t *w0_time = (double_t *)Calloc(M, sizeof(double_t));
    // Initialization with random values between 0 and 1
    initVectorWAnotherVector(w0_time, test_w0, M);
    //initializeRandomVector(w0_time, M);
    // vector by scalar product
    cblas_dscal(M, 0.07, w0_time, 1);
    // Add the scalar alpha at every element of the array w0_time
    sumScalarByVector(w0_time, M, 0.14f);
    //printDVector(w0_time, M, "w0_time");

    /************************************************************** 
     *  Create vector y0 = [U10;U20;W0] with initial conditions 
     * ***********************************************************/
    int y0Dimension;
    double_t *y0 = packThreeVectors(M, u10_time, u20_time, w0_time, &y0Dimension);
    //printDVector(y0, y0Dimension, "y0");

    /********************************************************************
     * Finite differences of order two and periodic boundary conditions
     * *****************************************************************/
    // Calculate Ldiff
    int sizeTempDiagOne, sizeTempDiagMinusOne;
    double_t *eyeM = eyeD(eyeM, M);
    eyeM = scalarByMatrix(eyeM, M, M, -2.0f);
    double_t *onesVector = onesD(onesVector, M - 1);
    double_t *tempDiagOne = diagD(onesVector, M - 1, 1, &sizeTempDiagOne);
    double_t *tempDiagMinusOne = diagD(onesVector, M - 1, -1, &sizeTempDiagMinusOne);
    double_t *addend1 = sumPuntSquareMatrices(eyeM, tempDiagOne, M);
    double_t *Ldiff = sumPuntSquareMatrices(addend1, tempDiagMinusOne, M);
    
    Ldiff = scalarByMatrix(Ldiff, M, M, 1.0f / (Delta_x * Delta_x));
    Ldiff[M - 1] = 1.0f / (Delta_x * Delta_x);
    Ldiff[(M - 1) * M] = 1.0f / (Delta_x * Delta_x);
    //printDMatrix(Ldiff, M, M, "Ldiff");

    // Calculate the matrix L
    int LSize;
    double_t *DLdiff = scalarByMatrix(Ldiff, M, M, D);
    double_t *dLdiff = scalarByMatrix(Ldiff, M, M, d);
    double_t *L = threeBlockDiagD(M, Ldiff, DLdiff, dLdiff, &LSize);
    //printDMatrix(L, LSize, LSize, "L");

    /*
    */
    int ySize;
    double_t *y = RungeKutta4th(2.0f, 0.0f, y0, y0Dimension, L, LSize, &ySize);
    // printDVector(y, ySize, "y");
    
    return_values result;
    // Initialize the returning pointers to NULL
    initReturnStruct(&result);
    // Compute using peer methods with stage = 2
    fPeerClassic_twoStages(N, t_span, 2, L, LSize, y0, y0Dimension, &result);

    // Using assertions to check if the pointers still NULL after the function invocation
    assert(result.t != NULL);
    assert(result.y != NULL);
    assert(result.yT != NULL);

    // Printing values
    printDVector(result.t, result.t_size, "t");
    printDVector(result.yT, result.yT_size, "yT_ClPeer");
    printDMatrix(result.y, result.y_rows, result.y_cols, "y_ClPeer");

    // Free all the memory dynamically allocated
    freeEverything(u10_time, u20_time, w0_time, y0, eyeM, onesVector, tempDiagOne, tempDiagMinusOne, 
                    addend1, Ldiff, DLdiff, dLdiff, L, result.t, result.y, result.yT, (void *)0);

    exit(0);
}