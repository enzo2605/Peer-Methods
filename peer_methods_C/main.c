/**
 * Title: Peer methods
 * Author: Vincenzo Iannucci
 * **/
#include "peerMethods.h"
#include <assert.h>

int main(int argc, char *argv[]) {
    // Random initialization for the seed
    srand((unsigned int)time(NULL));

    // Measure time taken by the whole algorithm
    clock_t start, end;

    // Intervals
    double t_start, t_end, x_start, x_end;
    char inputFileName[MAX_FILE_NAME_CHAR];
    char outputFileName[MAX_FILE_NAME_CHAR];

    // Checking arguments passed by the user
    if (argc != 7) {
        fprintf(stdout, "Usage: %s t_start t_end x_start x_end file_name\n", argv[0]);
        fprintf(stdout, "Using default parameters...\n\n");
        t_start = 0.0f;
        t_end   = 50.0f;
        x_start = -50.0f;
        x_end   = 50.0f;
        strcpy(inputFileName, "inputC.txt");
        strcpy(outputFileName, "outputC.txt");
    }
    else {
        t_start = atof(argv[1]);
        t_end   = atof(argv[2]);
        x_start = atof(argv[3]);
        x_end   = atof(argv[4]);
        strcpy(inputFileName, argv[5]);
        strcpy(outputFileName, argv[6]);
    }

    // Printing the interval used
    fprintf(stdout, "Interval used:\ntime span: [%f, %f]\nspace span: [%f, %f]\n", t_start, t_end, x_start, x_end);

    // Timer starting here
    start = clock();

    /*********************************************** 
     *          Time initialization 
     * *********************************************/
    double t_span[2] = { t_start, t_end };
    // for test
    double Delta_t = 1.0f / pow(2.0f, 8.0f);
    fprintf(stdout, "Delta_t: %f\n", Delta_t);

    double *t_int;
    int n_points_t;
    t_int = intervalDiscretization(t_int, t_start, t_end, Delta_t, &n_points_t);
    int N = (t_span[1] - t_span[0]) / Delta_t;
    fprintf(stdout, "N: %d\n", N);

    /*********************************************** 
     *          Space initialization 
     * *********************************************/
    double x_span[2] = { x_start, x_end };
    double Delta_x = (x_span[1] - x_span[0]) / M;
    printf("Delta_x: %f\n",Delta_x);

    double *x_int;
    int n_points_x;
    x_int = intervalDiscretization(x_int, x_start, x_end, Delta_x, &n_points_x);
    assert(n_points_x == M + 1);

    /*********************************************** 
     *          Define initial conditions 
     * *********************************************/

    /**** u10_time *****/
    // Dynamic allocation of the vectors
    double *u10_time = (double *)Calloc(M, sizeof(double));
    double *u20_time = (double *)Calloc(M, sizeof(double));
    double *w0_time  = (double *)Calloc(M, sizeof(double));

    // Initialize vectors by reading data from a file
    if (initInputVectors(inputFileName, u10_time, u20_time, w0_time, M) == 0) {
        fprintf(stdout, "\nVectors have been successfully initialized by the file %s.\n", inputFileName);
    }
    else {
        fprintf(stdout, "\nAn error occurred while initializing vectors.\n");
    }

    /************************************************************** 
     *  Create vector y0 = [U10;U20;W0] with initial conditions 
     * ***********************************************************/
    int y0Dimension;
    double *y0 = packThreeVectors(M, u10_time, u20_time, w0_time, &y0Dimension);
    //printDVector(y0, y0Dimension, "y0");

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
    //printDMatrix(Ldiff, M, M, "Ldiff");

    // Calculate the matrix L
    int LSize;
    double *DLdiff = scalarByMatrix(Ldiff, M, M, D);
    double *dLdiff = scalarByMatrix(Ldiff, M, M, d);
    double *L = threeBlockDiagD(M, Ldiff, DLdiff, dLdiff, &LSize);
    //printDMatrix(L, LSize, LSize, "L");
    
    return_values result;
    // Initialize the returning pointers to NULL
    initReturnStruct(&result);
    // Compute using peer methods with stage = 2
    fPeerClassic_twoStages(N, t_span, 2, L, LSize, y0, y0Dimension, &result);

    // Using assertions to check if the pointers still NULL after the function invocation
    assert(result.t != NULL);
    assert(result.y != NULL);
    assert(result.yT != NULL);

    // Saving data into a file
    if (saveInFile(outputFileName, result) == 0) {
        fprintf(stdout, "\nData have been successfully saved in the file %s.\n", outputFileName);
    }
    else {
        fprintf(stdout, "\nAn error occoured while saving data in file %s.\n", outputFileName);
    }

    // Timer finishing here
    end = clock();

    // Compute the time taken
    double cpuTimeUsed = ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(stdout, "\nTotal time: %lf\n", cpuTimeUsed);

    // Free all the memory dynamically allocated
    freeEverything(u10_time, u20_time, w0_time, y0, eyeM, onesVector, tempDiagOne, tempDiagMinusOne, 
                    addend1, Ldiff, DLdiff, dLdiff, L, result.t, result.y, result.yT, (void *)0);

    exit(0);
}