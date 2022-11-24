/**
 * Author: Vincenzo Iannucci
 * Purpose: The library provides an implementation for useful functions
 * like wrapper functions, functions used to dispay results and so on.
 * **/
#ifndef utilities_h
#define utilities_h

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>

#define MAX_FILE_NAME_CHAR 256

#define EVALUATE_EXECUTION_TIME(call, time)         \
{                                                   \
    struct timeval startModule, endModule;          \
    gettimeofday(&startModule, NULL);               \
    (call);                                         \
    gettimeofday(&endModule, NULL);                 \
    (time) = getTimeSpent(startModule, endModule);  \
}                                                   \

/*************************************************
 *          Wrapper functions
 ***********************************************/
void *Malloc(size_t size);
void *Calloc(size_t nmemb, size_t size);

/*************************************************
 *          Displaying stuff functions
 ***********************************************/
void printDMatrix(double *matrix, int M, int N, const char *string);
void printDMatrixRowMajor(double *matrix, int M, int N, const char *string);
void printDMatrixStacked(double *matrix, int M, int N, const char *string);
void printDVector(double *vector, int N, const char *string);

/*************************************************
 *          Initialization functions
 ***********************************************/
void initializeRandomVector(double *vector, int N);
void initializeRandomMatrix(double *matrix, int M, int N);
int initMatrixByRowWithValuesFromVector(double *matrix, int M, int N, double *vector, int vector_size);
void initVectorWAnotherVector(double *newVector, double *oldVector, int n);
double getTimeSpent(struct timeval start, struct timeval end);
int saveVectorsInFile(const char *fileName, int elements, double *arr1, int dim1, ...);
int initInputVectors(const char *fileName, double *u10_time, double *u20_time, double *w0_time, int dimension);
int saveMatrixInFile(const char *fileName, double *matrix, int matrix_rows, int matrix_cols);

/*************************************************
 *  Free all the memory dynamically allocated
 ***********************************************/
void freeEverything(void *arg1, ...);

#endif // !utilities_h