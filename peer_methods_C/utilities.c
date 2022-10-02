#include "utilities.h"

void *Malloc(size_t size) {
    void *res = malloc(size);
    if (res == NULL) {
        perror("malloc");
        exit(1);
    }
    return res;
}
       
void *Calloc(size_t nmemb, size_t size) {
    void *res = calloc(nmemb, size);
    if (res == NULL) {
        perror("calloc");
        exit(1);
    }
    return res;
}

void printDMatrix(double *matrix, int M, int N) {
    int i, j;
    fprintf(stdout, "Matrix size: %d x %d\n", M, N);
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            fprintf(stdout, "%8.2f", matrix[j * M + i]);
        }
        fprintf(stdout, "\n");
    }
}

void printDVector(double *vector, int N) {
    int i;
    fprintf(stdout, "Vector size: %d\n", N);
    for (i = 0; i < N; i++) {
        fprintf(stdout, "%8.2lf\n", vector[i]);
    }
    fprintf(stdout, "\n");
}

void initializeRandomVector(double *vector, int N) {
    int i;
    double randomNumber;

    for (i = 0; i < N; i++) {
        randomNumber = (double)rand() / ((double)RAND_MAX);
        vector[i] = randomNumber;
    }
}

void initializeRandomMatrix(double *matrix, int M, int N) {
    int k = 0;
    int i, j;
    // generation by columns starting by the index 1
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            matrix[k++] = (double)rand() / ((double)RAND_MAX);
        }
    }
}