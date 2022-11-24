/**
 * Author: Vincenzo Iannucci
 * Purpose: The library provides an implementation for the main function
 * for solving peer method.
 * **/
#ifndef peerMethods_h
#define peerMethods_h

#include "utilities.h"
#include "CLab.h"
#include "external_libs/CBLAS/include/cblas.h"

#define a  1.5
#define B1 0.45
#define B2 0.3611 
#define F  0.802 
#define H  0.802 
#define S  0.0002 
#define d  500 
#define D  0.802
#define M  4
#define STAGES 2

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
int saveVectorsInFile(const char *fileName, int elements, double *arr1, int dim1, ...);
int saveResultsInFile(const char* fileName, return_values result);
int initInputVectors(const char *fileName, double *u10_time, double *u20_time, double *w0_time, int dimension);
int saveMatrixInFile(const char *fileName, double *matrix, int matrix_rows, int matrix_cols);


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

#endif // !peerMethods_h