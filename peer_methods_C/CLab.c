#include "CLab.h"

double_t *intervalDiscretization(double_t *vector, double_t first, double_t last, double_t step, int *N) {
    // Number of elements
    int size = ((int)(last - first) / step) + 1;
    // Allocate the array using calloc
    vector = (double_t *)Calloc(size, sizeof(double_t));
    // Fill the array
    for (int i = 0; i < size; i++) {
        *(vector + i) = first;
        first += step;
    }
    *N = size;
    return vector;
}

double_t *eyeD(double_t *a, int N) {
    // Allocate the matrix, rembering that we need to allocate
    // the matrix by columns and start from 1
    a = (double_t *)Calloc(N * N, sizeof(double_t));
    int k = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            a[i * N + j] = (i == j) ? 1 : 0;
        }
    }
    return a;
}

double_t *onesD(double_t *a, int N) {
    // Allocate the array using calloc
    a = (double_t *)Calloc(N, sizeof(double_t));
    // Fill the array
    for (int i = 0; i < N; i++) {
        *(a + i) = 1.0f;
    }
    return a;
}

double_t *zerosD(int N) {
    // Allocate the array using calloc and initialize automatically
    // every element with 0
    double_t *a = (double_t *)Calloc(N, sizeof(double_t));
    return a;
}

double_t *zerosMatrixD(int M, int N) {
    // Allocate the array using calloc and initialize automatically
    // every element with 0
    double_t *a = (double_t *)Calloc(M * N, sizeof(double_t));
    return a;
}

void sumScalarByVector(double_t *in, int N, double_t alpha) {
    int i;
    for (i = 0; i < N; i++) {
        in[i] += alpha;
    }
}

void scalarByVector(double_t *in, int N, double_t alpha) {
    int i;
    for (i = 0; i < N; i++) {
        in[i] *= alpha;
    }
}

double_t *sumPuntVectors(double_t *v1, double_t *v2, int size) {
    double_t *res = (double_t *)Calloc(size, sizeof(double_t));

    for (int i = 0; i < size; i++) {
        res[i] = v1[i] + v2[i];
    }

    return res;
}

double_t *diagD(double_t *vector, int size, int k, int *matrix_size) {
    // The dimension of the final matrix will be N x N
    // with N = k + 1
    int N = size + abs(k);

    // Allocate dynamically the new matrix
    double_t *matrix = (double_t *)Calloc(N * N, sizeof(double_t));

    // Decide where to place arguments based on the absoulute value of k
    int j = abs(k);
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

double_t *packThreeMatrices(int n, double_t *A, double_t *B, double_t *C) {
    // Allocates the new matrix
    double_t *pack = (double_t *)Calloc(n * n * 3, sizeof(double_t));

    // Copy the three matrix one side by another
    memcpy(pack, A, n * n * sizeof(double_t));
    memcpy(pack + n * n, B, n * n * sizeof(double_t));
    memcpy(pack + n * n + n * n, C, n * n * sizeof(double_t));

    return pack;
}

double_t *threeBlockDiagD(int n, double_t *A, double_t *B, double_t *C, int *blckSize) {
    // The new size will be the old one multiplied by 3
    int newSize = n * 3;

    // Allocates the final matrix
    double_t *blockMatrix = (double_t *)Calloc(newSize * newSize, sizeof(double_t));
    // Pack the three matrix into one using an apposite function
    double_t *pack = packThreeMatrices(n, A, B, C);

    // Copy row by row the element of pack into blockMatrix
    for (int i = 0; i < newSize; i += n) {
        for (int j = 0; j < n; j++) {
            memcpy(blockMatrix + (i * newSize + i) + j * newSize, pack + (i * n) + j * n, n * sizeof(double_t));
        }
    }
    *blckSize = newSize;
    return blockMatrix;
}

double_t *packThreeVectors(int n, double_t *A, double_t *B, double_t *C, int *newDimension) {
    int newSize = n * 3;
    
    // Allocates the new vector
    double_t *pack = (double_t *)Calloc(newSize, sizeof(double_t));

    // Copy the three vector one side by another
    memcpy(pack, A, n * sizeof(double_t));
    memcpy(pack + n, B, n * sizeof(double_t));
    memcpy(pack + 2 * n, C, n * sizeof(double_t));

    // Saving the new dimension of the vector
    *newDimension = newSize;

    return pack;
}

double_t *sumPuntSquareMatrices(double_t *matrix1, double_t *matrix2, int size) {
    int i, j;

    // Allocate the resulting matrix
    double_t *result = (double_t *)Calloc(size * size, sizeof(double_t));

    // Sum element by element of the matrices
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            result[i * size + j] = matrix1[i * size + j] + matrix2[i * size + j];
        }
    }

    return result;
}

double_t *scalarByMatrix(double_t *matrix, int M, int N, double_t alpha) {
    int i, j;
    double_t *newMatrix = (double_t *)Calloc(M * N, sizeof(double_t));

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            newMatrix[i * N + j] = matrix[i * N + j] * alpha;
        }
    }

    return newMatrix;
}

double_t *linspace(double_t x1, double_t x2, int n) {
    // Allocate the vector of size n dynamically
    double_t *v = (double_t *)Calloc(n, sizeof(double_t));
    // Define the spacing between eache element of the vector
    double_t step = (x2 - x1) / (n - 1);
    // Generate the values of the vector
    v[0] = x1;
    for (int i = 1; i < n; i++) {
        v[i] = v[i - 1] + step;
    }

    return v;
}