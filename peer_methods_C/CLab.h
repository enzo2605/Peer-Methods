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
double *setDVector(double *vector, double first, double last, double step, int *N);

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
 * Sum each element of vector in by the scalar alpha.
 * Return the pointer at the resulting array. 
 *****************************************************/
double *sumScalarByVector(double *in, int N, double alpha);

/*******************************************************
 * Return a matrix matrix_size x matrix_size elements in 
 * which each element of the k-th diagonal is an element 
 * of the vector passed by parameter.
 *****************************************************/
double *diagD(double *vector, int size, int k, int *matrix_size);

/*******************************************************
 * Packing three matrices side by side into one. It's
 * an utility function for the threeBlockDiagD function
 *****************************************************/
double *packMatrices(int n, double *A, double *B, double *C);

/*******************************************************
 * Return a block diag matrix of three matrices passed
 * by arguments of size n x n.
 *****************************************************/
double *threeBlockDiagD(int n, double *A, double *B, double *C);

#endif // !CLab_h