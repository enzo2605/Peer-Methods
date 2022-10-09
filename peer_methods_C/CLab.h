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
double_t *intervalDiscretization(double_t *vector, double_t first, double_t last, double_t step, int *N);

/*******************************************************
 * Return a square matrix of size N in which the  
 * elemets at the position i == j assume value 1
 *****************************************************/
double_t *eyeD(double_t *a, int N);

/*******************************************************
 * Return a vector of size N in which every element
 * has value 1.
 *****************************************************/
double_t *onesD(double_t *a, int N);

/*******************************************************
 * Return a vector of size N in which every element
 * has value 0.
 *****************************************************/
double_t *zerosD(int N);

/*******************************************************
 * Return a matrix of size M * N in which every element
 * has value 0.
 *****************************************************/
double_t *zerosMatrixD(int M, int N);

/*******************************************************
 * Sum each element of vector in by the scalar alpha.
 * Return the pointer at the resulting array. 
 *****************************************************/
void sumScalarByVector(double_t *in, int N, double_t alpha);

/*******************************************************
 * Multiply each element of vector in by the scalar alpha.
 * Return the pointer at the resulting array. 
 *****************************************************/
void scalarByVector(double_t *in, int N, double_t alpha);

/*******************************************************
 * Sum each element of vector v1 and v2 in the same 
 * position and store the result in res vector. 
 *****************************************************/
double_t *sumPuntVectors(double_t *v1, double_t *v2, int size);

/*******************************************************
 * Return a matrix of matrix_size x matrix_size elements 
 * in which each element of the k-th diagonal is an 
 * element of the vector passed by parameter.
 *****************************************************/
double_t *diagD(double_t *vector, int size, int k, int *matrix_size);

/*******************************************************
 * Packing three matrices side by side into one. It's
 * an utility function for the threeBlockDiagD function
 *****************************************************/
double_t *packThreeMatrices(int n, double_t *A, double_t *B, double_t *C);

/*******************************************************
 * Return a block diag matrix of three matrices passed
 * by arguments of size n x n.
 *****************************************************/
double_t *threeBlockDiagD(int n, double_t *A, double_t *B, double_t *C, int *blckSize);

/*******************************************************
 * Packing three vectors side by side into one.
 *****************************************************/
double_t *packThreeVectors(int n, double_t *A, double_t *B, double_t *C, int *newDimension);

/*******************************************************
 * Sum element by element of two square matrices
 * of the same dimension.
 *****************************************************/
double_t *sumPuntSquareMatrices(double_t *matrix1, double_t *matrix2, int size);

/*******************************************************
 * Multiply a scalar value for a matrix.
 *****************************************************/
double_t *scalarByMatrix(double_t *matrix, int M, int N, double_t alpha);

/****************************************************
 * Generates n points. The spacing between the points 
 * is (x2-x1)/(n-1)
 ***************************************************/
double_t *linspace(double_t x1, double_t x2, int n);

#endif // !CLab_h