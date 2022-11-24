#ifndef helper_h
#define helper_h

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>

#define MAX_FILE_NAME_CHAR 256

#define EVALUATE_EXECUTION_TIME(call, time)         \
{                                                   \
    struct timeval startModule, endModule;          \
    gettimeofday(&startModule, NULL);               \
    (call);                                         \
    gettimeofday(&endModule, NULL);                 \
    (time) = getTimeSpent(startModule, endModule);  \
}                                                   \

#ifdef __cplusplus
extern "C" {
#endif

/*************************************************
 *          Displaying stuff functions
 ***********************************************/
void printDMatrix(double *matrix, int M, int N, const char *string);
void printDMatrixRowMajor(double *matrix, int M, int N, const char *string);
void printDMatrixStacked(double *matrix, int M, int N, const char *string);
void printDVector(double *vector, int N, const char *string);

/*************************************************
 *          Saving file function
 ***********************************************/
int saveVectorsInFile(const char *fileName, int elements, double *arr1, int dim1, ...);
int initInputVectors(const char *fileName, double *u10_time, double *u20_time, double *w0_time, int dimension);
int saveMatrixInFile(const char *fileName, double *matrix, int matrix_rows, int matrix_cols);

/* Timing utility function */
double getTimeSpent(struct timeval start, struct timeval end);

#ifdef __cplusplus
}
#endif
#endif