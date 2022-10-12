/**
 * Author: Vincenzo Iannucci
 * Purpose: The library provides an implementation for the main function
 * for solving peer method.
 * **/
#ifndef peerMethods_h
#define peerMethods_h

#include "utilities.h"
#include "CLab.h"
#include <cblas.h>

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
    double_t *yT;
    int yT_size;
    double_t *y;
    int y_rows;
    int y_cols;
    double_t *t;
    int t_size;
} return_values;


/*********************************************************************************
 *                          Utilities functions
 ********************************************************************************/
int initInputVectors(const char *fileName, double_t *u10_time, double_t *u20_time, double_t *w0_time, int dimension);
int saveInFile(const char* fileName, return_values result);
void initReturnStruct(return_values *rv);

/*********************************************************************************
 * Given the array y0 with y0Size elements and the matrix L with Lsize elements
 * return the pointer to an array of size sherratSize.
 ********************************************************************************/
double_t *Sherratt(double_t *y0, int y0Size, double_t *L, int Lsize, int *sherrattSize);

/*********************************************************************************
 * Given the increase h, the initial time t0 and the vector y0 with ySize elements
 * and the matrix L with Lsize elements, return the vector y with ySize element.
 ********************************************************************************/
double_t *RungeKutta4th(double_t h, double_t t0, double_t *y0, int y0Size, double_t *L, int Lsize, int *ySize);

/*********************************************************************************
 * Given in input the time span, the L matrix and y0 vector with their relative
 * sizes, returns a struct called "return_values" that contains the result
 * of the computation.
 ********************************************************************************/
void fPeerClassic_twoStages(int N, double_t *t_span, int t_span_size, double_t *L, int Lsize, double_t *y0, int y0_size, return_values *collect_result);

#endif // !peerMethods_h