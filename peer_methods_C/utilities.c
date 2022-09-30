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
    printf("Matrix size: %d x %d\n", M, N);
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            printf("%8.2f", matrix[j * M + i]);
        }
        printf("\n");
    }
}

void printDVector(double *vector, int N) {
    int i;
    printf("Vector size: %d\n", N);
    for (i = 0; i < N; i++) {
        printf("%8.2lf\n", vector[i]);
    }
    printf("\n");
}