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
    fprintf(stdout, "N: %d\n", N);

    // Allocate dynamically the new matrix
    double *matrix = (double *)Calloc(N * N, sizeof(double));

    // Decide where to place arguments based on the absoulute value of k
    int j = abs(k);
    fprintf(stdout, "j: %d\n", j);
    for (int i = 0; i < size; i++) {
        if (k > 0) {
            matrix[j++ * N + i] = vector[i];
        }
        else if (k < 0) {
            matrix[i * N + j++] = vector[i];
        }
        else {
            matrix[i * N + i] = vector[i];
        }
    }

    // Saving the new matrix size
    *matrix_size = N;

    return matrix;
}

double *packThreeMatrices(int n, double *A, double *B, double *C) {
    // Allocates the new matrix
    double *pack = (double *)Calloc(n * n * 3, sizeof(double));

    // Copy the three matrix one side by another
    memcpy(pack, A, n * n * sizeof(double));
    memcpy(pack + n * n, B, n * n * sizeof(double));
    memcpy(pack + n * n + n * n, C, n * n * sizeof(double));

    return pack;
}

double *threeBlockDiagD(int n, double *A, double *B, double *C) {
    // The new size will be the old one multiplied by 3
    int newSize = n * 3;

    // Allocates the final matrix
    double *blockMatrix = (double *)Calloc(newSize * newSize, sizeof(double));
    // Pack the three matrix into one using an apposite function
    double *pack = packThreeMatrices(n, A, B, C);

    // Copy row by row the element of pack into blockMatrix
    for (int i = 0; i < newSize; i += n) {
        for (int j = 0; j < n; j++) {
            memcpy(blockMatrix + (i * newSize + i) + j * newSize, pack + (i * n) + j * n, n * sizeof(double));
        }
    }

    return blockMatrix;
}

double *packThreeVectors(int n, double *A, double *B, double *C, int *newDimension) {
    int newSize = n * 3;
    
    // Allocates the new vector
    double *pack = (double *)Calloc(newSize, sizeof(double));

    // Copy the three vector one side by another
    memcpy(pack, A, n * sizeof(double));
    memcpy(pack + n, B, n * sizeof(double));
    memcpy(pack + 2 * n, C, n * sizeof(double));

    // Saving the new dimension of the vector
    *newDimension = newSize;

    return pack;
}

double *sumPuntSquareMatrices(double *matrix1, double *matrix2, int size) {
    int i, j;

    // Allocate the resulting matrix
    double *result = (double *)Calloc(size * size, sizeof(double));

    // Sum element by element of the matrices
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            result[i * size + j] = matrix1[i * size + j] + matrix2[i * size + j];
        }
    }

    return result;
}

void scalarByMatrix(double *matrix, int M, int N, double alpha) {
    int i, j;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            matrix[j * M + i] *= alpha;
        }
    } 
}