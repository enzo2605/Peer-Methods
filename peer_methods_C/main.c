/**
 * Title: Peer methods
 * Author: Vincenzo Iannucci
 * **/
#include <stdio.h>
#include <stdlib.h>
#include "helper.h"
#include "peerMethods.h"

double a  = 1.5;
double B1 = 0.45;
double B2 = 0.3611;
double F  = 0.802;
double H  = 0.802;
double S  = 0.0002;
double d  = 500;
double D  = 0.802;
int M     = 4;

int main(int argc, char *argv[]) {
    // Random initialization for the seed
    srand((unsigned int)time(NULL));

    // Intervals
    double t_start, t_end, x_start, x_end, val;
    char inputFileName[MAX_FILE_NAME_CHAR];
    char outputFileName[MAX_FILE_NAME_CHAR];
    t_start = 0.0f;
    t_end   = 50.0f;
    x_start = -50.0f;
    x_end   = 50.0f;
    strcpy(inputFileName, "inputC.txt");
    strcpy(outputFileName, "../file_analyzer/outputC.txt");

    if (argc != 2) {
        fprintf(stdout, "Usage: %s val\n", argv[0]);
        val = 11.0f;
        //fprintf(stdout, "Using default parameters...\nDelta_t: 1 / 2^%f\n", val);
    }
    else {
        val = atof(argv[1]);
        //fprintf(stdout, "Delta_t: 1 / 2^%f\n", val);
    }

    /*
    // Checking arguments passed by the user
    if (argc != 7) {
        fprintf(stdout, "Usage: %s t_start t_end x_start x_end file_name\n", argv[0]);
        fprintf(stdout, "Using default parameters...\n\n");
        t_start = 0.0f;
        t_end   = 50.0f;
        x_start = -50.0f;
        x_end   = 50.0f;
        strcpy(inputFileName, "inputC.txt");
        strcpy(outputFileName, "../file_analyzer/outputC.txt");
    }
    else {
        t_start = atof(argv[1]);
        t_end   = atof(argv[2]);
        x_start = atof(argv[3]);
        x_end   = atof(argv[4]);
        strcpy(inputFileName, argv[5]);
        strcpy(outputFileName, argv[6]);
    }
    */

    // Printing the interval used
    fprintf(stdout, "\nInterval used:\ntime span: [%f, %f]\nspace span: [%f, %f]\n", t_start, t_end, x_start, x_end);
    fprintf(stdout, "\nM: %d\n", M);

    /*********************************************** 
     *          Time initialization 
     * *********************************************/
    double t_span[2] = { t_start, t_end };
    // for test
    double Delta_t = 1.0f / pow(2.0f, val);
    fprintf(stdout, "\nDelta_t: 1.0f / 2.0f ^ %f = %f\n", val, Delta_t);

    double *t_int;
    int n_points_t;
    t_int = intervalDiscretization(t_start, t_end, Delta_t, &n_points_t);
    int N = (t_span[1] - t_span[0]) / Delta_t;
    fprintf(stdout, "\nN: %d\n", N);

    /*********************************************** 
     *          Space initialization 
     * *********************************************/
    double x_span[2] = { x_start, x_end };
    double Delta_x = (x_span[1] - x_span[0]) / M;
    printf("\nDelta_x: %lf\n",Delta_x);

    double *x_int;
    int n_points_x;
    x_int = intervalDiscretization(x_start, x_end, Delta_x, &n_points_x);

    /*********************************************** 
     *          Define initial conditions 
     * *********************************************/

    /**** u10_time *****/
    // Dynamic allocation of the vectors
    double *u10_time = (double *)Calloc(M, sizeof(double));
    double *u20_time = (double *)Calloc(M, sizeof(double));
    double *w0_time  = (double *)Calloc(M, sizeof(double));
    //printDVector(u10_time, M, "u10_time");
    //printDVector(u20_time, M, "u20_time");
    //printDVector(w0_time, M, "w0_time");

    // Initialize vectors by reading data from a file
    if (initInputVectors(inputFileName, u10_time, u20_time, w0_time, M) == 0) {
        fprintf(stdout, "\nVectors have been successfully initialized by the file %s.\n", inputFileName);
    }
    else {
        fprintf(stdout, "\nAn error occurred while initializing vectors.\n");
    }

    // Measure the time taken by the whole software
    struct timeval startTotal, endTotal;

    // ------------------------ Total timer starting here ----------------- 
    gettimeofday(&startTotal, NULL);

    /************************************************************** 
     *  Create vector y0 = [U10;U20;W0] with initial conditions 
     * ***********************************************************/
    int y0Dimension;
    double *y0 = packThreeVectors(M, u10_time, u20_time, w0_time, &y0Dimension);
    //printDVector(y0, y0Dimension, "y0");

    /********************************************************************
     * Finite differences of order two and periodic boundary conditions
     * *****************************************************************/
    double *L;
    int LSize;
    double LMatrixCPUTimeUsed = 0.0f;
    EVALUATE_EXECUTION_TIME(computeLMatrix(&L, &LSize, Delta_x), LMatrixCPUTimeUsed);
    fprintf(stdout, "\nTime taken by computeLMatrix module: %lf\n", LMatrixCPUTimeUsed);
    //saveMatrixInFile("../parallel/good.txt", L, LSize, LSize);
    //exit(0);

    /*
    double *sherrattResult;
    int sherrattSize;
    double sherrattCPUTimeUsed = 0.0f;
    EVALUATE_EXECUTION_TIME(sherrattResult = Sherratt(y0, y0Dimension, L, LSize, &sherrattSize), sherrattCPUTimeUsed);

    double *rungeKuttaResult;
    int rungeKuttaSize;
    double RungeKutta4thCPUTimeUsed = 0.0f;
    double h = (t_span[1] - t_span[0]) / N;
    printf("h: %lf\n", h);
    EVALUATE_EXECUTION_TIME(rungeKuttaResult = RungeKutta4th(h, 0, y0, y0Dimension, L, LSize, &rungeKuttaSize), RungeKutta4thCPUTimeUsed);
    fprintf(stdout, "\nTime taken by RungeKutta4th module: %lf\n", RungeKutta4thCPUTimeUsed);
    saveVectorsInFile("../parallel/good.txt", 3, y0, y0Dimension, L, LSize * LSize, rungeKuttaResult, rungeKuttaSize, (void *)0);
    exit(0);
    */
    
    return_values result;
    // Initialize the returning pointers to NULL
    initReturnStruct(&result);
    double fPeerClassicCPUTimeUsed = 0.0f;
    // Compute using peer methods with stage = 2
    EVALUATE_EXECUTION_TIME(fPeerClassic_twoStages(N, t_span, 2, L, LSize, y0, y0Dimension, &result), fPeerClassicCPUTimeUsed);
    fprintf(stdout, "\nTime taken by fPeerClassic_twoStages module: %lf\n", fPeerClassicCPUTimeUsed);

    // ------------------------ Total timer ending here ----------------- 
    gettimeofday(&endTotal, NULL);
    double totalExecutionTime = getTimeSpent(startTotal, endTotal);
    printf("Total execution time: %lf\n", totalExecutionTime);

    /*
    // Saving data into a file
    if (saveInFile(outputFileName, result) == 0) {
        fprintf(stdout, "\nData have been successfully saved in the file %s.\n", outputFileName);
    }
    else {
        fprintf(stdout, "\nAn error occoured while saving data in the file %s.\n", outputFileName);
    }
    */

    // Free all the memory dynamically allocated
    freeEverything(u10_time, u20_time, w0_time, y0, L, result.t, result.y, result.yT, (void *)0);

    exit(0);
}