/**
 * Author: Vincenzo Iannucci
 * Purpose: The library provides an implementation for some MatLab routines
 * that are useful for the project.
 * **/

#ifndef CLab_h
#define CLab_h

#include "utilities.h"

/*******************************************************
 * Return a vector in which the first element is 
 * first, the last element is last with a step that is 
 * represented by step parameters. N represent the
 * size of the array at the end.
 *****************************************************/
double *intervalDiscretization(double *vector, double first, double last, double step, int *N);

/*******************************************************
 * Return a square matrix of size N in which the  
 * elemets at the position i == j assume value 1
 *****************************************************/
double *eyeD(double *a, int N);

/*******************************************************
 * Return a vector of size N in which every element
 * has value 1.
 *****************************************************/
double *onesD(double *a, int N);

/*******************************************************
 * Return a vector of size N in which every element
 * has value 0.
 *****************************************************/
double *zerosD(double *a, int N);

/*******************************************************
 * Return a matrix of size M * N in which every element
 * has value 0.
 *****************************************************/
double *zerosMatrixD(int M, int N);

/*******************************************************
 * Sum each element of vector in by the scalar alpha.
 * Return the pointer at the resulting array. 
 *****************************************************/
void sumScalarByVector(double *in, int N, double alpha);

/*******************************************************
 * Multiply each element of vector in by the scalar alpha.
 * Return the pointer at the resulting array. 
 *****************************************************/
void scalarByVector(double *in, int N, double alpha);

/*******************************************************
 * Sum each element of vector v1 and v2 in the same 
 * position and store the result in res vector. 
 *****************************************************/
double *sumPuntVectors(double *v1, double *v2, int size);

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

/*******************************************************
 * Sum element by element of two square matrices
 * of the same dimension.
 *****************************************************/
double *sumPuntSquareMatrices(double *matrix1, double *matrix2, int size);

/*******************************************************
 * Multiply a scalar value for a matrix.
 *****************************************************/
double *scalarByMatrix(double *matrix, int M, int N, double alpha);

/****************************************************
 * Generates n points. The spacing between the points 
 * is (x2-x1)/(n-1)
 ***************************************************/
double *linspace(double x1, double x2, int n);

#endif // !CLab_h