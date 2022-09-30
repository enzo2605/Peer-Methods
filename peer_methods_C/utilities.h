#ifndef utilities_h
#define utilities_h

#include <stdio.h>
#include <stdlib.h>

/*************************************************
 *          Wrapper functions
 ***********************************************/
void *Malloc(size_t size);
void *Calloc(size_t nmemb, size_t size);

/*************************************************
 *          Displaying stuff functions
 ***********************************************/
void printDMatrix(double *matrix, int M, int N);
void printDVector(double *vector, int N);

#endif // !utilities_h