/**
 * Author: Vincenzo Iannucci
 * Purpose: The library provides an implementation for the main function
 * for solving peer method.
 * **/
#ifndef peerMethods_h
#define peerMethods_h

/*
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include "external_libs/CBLAS/include/cblas.h"
*/

#define STAGES 2

#ifdef __cplusplus
extern "C" {
#endif

extern double a, B1, B2, F, H, S, d, D, L;
extern int M;

/*********************************************************************************
 * This struct has been created with the only purpose to return the value
 * obtained by the fPeerMethod function.
 ********************************************************************************/
typedef struct {
    double *yT;
    int yT_size;
    double *y;
    int y_rows;
    int y_cols;
    double *t;
    int t_size;
} return_values;


/*********************************************************************************
 *                          Utilities functions
 ********************************************************************************/
void initReturnStruct(return_values *rv);
int saveResultsInFile(const char* fileName, return_values result);


/*********************************************************************************
 * Given a double pointer to matrix L, a pointer to its size and Delta_x, the 
 * function compute the matrix L of size LSize * LSize.
 ********************************************************************************/
void computeLMatrix(double **L, int *LSize, double Delta_x);

/*********************************************************************************
 * Given the array y0 with y0Size elements and the matrix L with Lsize elements
 * return the pointer to an array of size sherratSize.
 ********************************************************************************/
double *Sherratt(const double *y0, int y0Size, const double *L, int Lsize, int *sherrattSize);

/*********************************************************************************
 * Given the increase h, the initial time t0 and the vector y0 with ySize elements
 * and the matrix L with Lsize elements, return the vector y with ySize element.
 ********************************************************************************/
double *RungeKutta4th(double h, double t0, const double *y0, int y0Size, const double *L, int Lsize, int *ySize);

/*********************************************************************************
 * Given in input the time span, the L matrix and y0 vector with their relative
 * sizes, returns a struct called "return_values" that contains the result
 * of the computation.
 ********************************************************************************/
void fPeerClassic_twoStages(int N, double *t_span, int t_span_size, const double *L, int Lsize, const double *y0, int y0_size, return_values *collect_result);

/*************************************************
 *          Wrapper functions
 ***********************************************/
void *Malloc(size_t size);
void *Calloc(size_t nmemb, size_t size);

/*************************************************
 *          Initialization functions
 ***********************************************/
void initializeRandomVector(double *vector, int N);
void initializeRandomMatrix(double *matrix, int M, int N);
int initMatrixByRowWithValuesFromVector(double *matrix, int M, int N, double *vector, int vector_size);
void initVectorWAnotherVector(double *newVector, double *oldVector, int n);

/*************************************************
 *  Free all the memory dynamically allocated
 ***********************************************/
void freeEverything(void *arg1, ...);

/*******************************************************
 * Return a vector in which the first element is 
 * first, the last element is last with a step that is 
 * represented by step parameters. N represent the
 * size of the array at the end.
 *****************************************************/
double *intervalDiscretization(double first, double last, double step, int *N);

/*******************************************************
 * Return a square matrix of size N in which the  
 * elemets at the position i == j assume value 1
 *****************************************************/
double *eyeD(int N);

/*******************************************************
 * Return a vector of size N in which every element
 * has value 1.
 *****************************************************/
double *onesD(int N);

/*******************************************************
 * Return a vector of size N in which every element
 * has value 0.
 *****************************************************/
double *zerosD(int N);

/*******************************************************
 * Return a matrix of size M * N in which every element
 * has value 0.
 *****************************************************/
double *zerosMatrixD(int M, int N);

/*******************************************************
 * Return a matrix of matrix_size x matrix_size elements 
 * in which each element of the k-th diagonal is an 
 * element of the vector passed by parameter.
 *****************************************************/
double *diagD(double *vector, int size, int k, int *matrix_size);

/*******************************************************
 * Packing three matrices side by side into one. It's
 * an utility function for the threeBlockDiagD function
 *****************************************************/
double *packThreeMatrices(int n, double *A, double *B, double *C);

/*******************************************************
 * Return a block diag matrix of three matrices passed
 * by arguments of size n x n.
 *****************************************************/
double *threeBlockDiagD(int n, double *A, double *B, double *C, int *blckSize);

/*******************************************************
 * Packing three vectors side by side into one.
 *****************************************************/
double *packThreeVectors(int n, double *A, double *B, double *C, int *newDimension);

/****************************************************
 * Generates n points. The spacing between the points 
 * is (x2-x1)/(n-1)
 ***************************************************/
double *linspace(double x1, double x2, int n);

#ifdef __cplusplus
}
#endif
#endif // !peerMethods_h