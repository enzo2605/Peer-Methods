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

double *eyeD(double *a, int N) {
    // Allocate the matrix, rembering that we need to allocate
    // the matrix by columns and start from 1
    a = (double *)Calloc(N * N + 1, sizeof(double));
    int k = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            a[k++] = (i == j) ? 1 : 0;
        }
    }
    return a;
}

double *onesD(double *a, int N) {
    // Allocate the array using calloc
    a = (double *)Calloc(N, sizeof(double));
    // Fill the array
    for (int i = 0; i < N; i++) {
        *(a + i) = 1.0f;
    }
    return a;
}

double *zerosD(double *a, int N) {
    // Allocate the array using calloc
    a = (double *)Calloc(N, sizeof(double));
    // Fill the array
    for (int i = 0; i < N; i++) {
        *(a + i) = 0.0f;
    }
    return a;
}

double *sumScalarByVector(double *in, int N, double alpha) {
    int i;
    for (i = 0; i < N; i++) {
        in[i] += alpha;
    }
    return in;
}

double *diagD(double *vector, int size, int k, int *matrix_size) {
    // The dimension of the final matrix will be N x N
    // with N = k + 1
    int N = size + abs(k);
    double *matrix = (double *)Calloc(N * N, sizeof(double));

    int j = abs(k);
    for (int i = 0; i < size; i++) {
        if (k > 0) {
            matrix[k++ * N + i] = vector[i];
        }
        else if (k < 0) {
            matrix[i * N + j++] = vector[i];
        }
        else {
            matrix[i * N + i] = vector[i];
        }
    }

    *matrix_size = N;
    return matrix;
}

double *threeBlockDiagD(int n, double *A, double *B, double *C) {
    
}