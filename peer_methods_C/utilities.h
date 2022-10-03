/**
 * Author: Vincenzo Iannucci
 * Purpose: The library provides an implementation for useful functions
 * like wrapper functions, functions used to dispay results and so on.
 * **/

#ifndef utilities_h
#define utilities_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>

/*************************************************
 *          Wrapper functions
 ***********************************************/
void *Malloc(size_t size);
void *Calloc(size_t nmemb, size_t size);

/*************************************************
 *          Displaying stuff functions
 ***********************************************/
void printDMatrix(double *matrix, int M, int N, const char *string);
void printDVector(double *vector, int N);

/*************************************************
 *          Initialization functions
 ***********************************************/
void initializeRandomVector(double *vector, int N);
void initializeRandomMatrix(double *matrix, int M, int N);

/*************************************************
 *  Free all the memory dynamically allocated
 ***********************************************/
void freeEverything(void *arg1, ...);

#endif // !utilities_h