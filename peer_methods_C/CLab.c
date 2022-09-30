#include "CLab.h"

double *setDVector(double *vector, double first, double last, double step, int *N) {
    // Number of elements
    int size = ((int)(last - first) / step) + 1;
    // Allocate the array using calloc
    vector = (double *)Calloc(size, sizeof(double));
    // Fill the array
    for (int i = 0; i < size; i++) {
        *(vector + i) = first;
        first += step;
    }
    *N = size;
    return vector;
}