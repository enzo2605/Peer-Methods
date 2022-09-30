/**
 * Author: Vincenzo Iannucci
 * Purpose: The library provides an implementation for some MatLab routines
 * that are useful for the project.
 * **/

#ifndef CLab_h
#define CLab_h

#include "utilities.h"

/*******************************************************
 * Auto generate a vector in which the first element is 
 * first, the last element is last with a step that is 
 * represented by step parameters.
 *****************************************************/
double *setDVector(double *vector, double first, double last, double step, int *N);

#endif // !CLab_h