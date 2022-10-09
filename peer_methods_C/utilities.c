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

void printDMatrix(double_t *matrix, int M, int N, const char *string) {
    int i, j;
    fprintf(stdout, "\n%s\nMatrix size: %d x %d\n", string, M, N);
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            fprintf(stdout, "%10.4f", matrix[j * M + i]);
        }
        fprintf(stdout, "\n");
    }
}

void printDVector(double_t *vector, int N, const char *string) {
    int i;
    fprintf(stdout, "\n%s\nVector size: %d\n", string, N);
    for (i = 0; i < N; i++) {
        fprintf(stdout, "%10.5lf\n", vector[i]);
    }
    fprintf(stdout, "\n");
}

void initializeRandomVector(double_t *vector, int N) {
    int i;
    double_t randomNumber;

    for (i = 0; i < N; i++) {
        randomNumber = (double_t)rand() / ((double_t)RAND_MAX);
        vector[i] = randomNumber;
    }
}

void initializeRandomMatrix(double_t *matrix, int M, int N) {
    int k = 0;
    int i, j;
    // generation by columns starting by the index 1
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            matrix[i * N + j] = (double_t)rand() / ((double_t)RAND_MAX);
        }
    }
}

int initMatrixByRowWithValuesFromVector(double_t *matrix, int M, int N, double_t *vector, int vector_size) {
    if (M * N != vector_size) {
        fprintf(stdout, "\nIt has been impossible to initialize the matrix with the vector passed.");
        fprintf(stdout, "\nThe size of matrix and the vector are incompatible.\nReturned -1\n");
        return -1;
    }

    int k = 0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            matrix[j * M + i] = vector[k++];
        }
    }
    return 0;
}

void initVectorWAnotherVector(double_t *newVector, double_t *oldVector, int n) {
    int i;
    for (i = 0; i < n; i++) {
        newVector[i] = oldVector[i];
    }
}

void freeEverything(void *arg1, ...) {
    va_list args; // list of arguments
    void *vp; // pointer to the i-th arguments
    // Free the first argument
    free(arg1);
    va_start(args, arg1);
    // Until you reach the end of the arguments, free the argument
    // pointed by vp
    while ((vp = va_arg(args, void *)) != 0) {
        free(vp);
    }
    va_end(args);
}