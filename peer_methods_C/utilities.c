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

void printDMatrix(double *matrix, int M, int N, const char *string) {
    int i, j;
    fprintf(stdout, "\n%s\nMatrix size: %d x %d\n", string, M, N);
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            fprintf(stdout, " %.15lf", matrix[j * M + i]);
        }
        fprintf(stdout, "\n");
    }
}

void printDMatrixRowMajor(double *matrix, int M, int N, const char *string) {
    int i, j;
    fprintf(stdout, "\n%s\nMatrix size: %d x %d\n", string, M, N);
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            fprintf(stdout, " %.15f", matrix[i * N + j]);
        }
        fprintf(stdout, "\n");
    }
}

void printDMatrixStacked(double *matrix, int M, int N, const char *string) {
    int i, j;
    fprintf(stdout, "\n%s\nMatrix size: %d x %d\n", string, M, N);
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            fprintf(stdout, "%.15lf\n", matrix[i * N + j]);
        }
        //fprintf(stdout, "\n");
    }
}

void printDVector(double *vector, int N, const char *string) {
    int i;
    fprintf(stdout, "\n%s\nVector size: %d\n", string, N);
    for (i = 0; i < N; i++) {
        // 10.5lf
        fprintf(stdout, "%.15lf\n", vector[i]);
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
    int i, j;
    // generation by columns starting by the index 1
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            matrix[i * N + j] = (double)rand() / ((double)RAND_MAX);
        }
    }
}

int initMatrixByRowWithValuesFromVector(double *matrix, int M, int N, double *vector, int vector_size) {
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

void initVectorWAnotherVector(double *newVector, double *oldVector, int n) {
    int i;
    for (i = 0; i < n; i++) {
        newVector[i] = oldVector[i];
    }
}

inline double getTimeSpent(struct timeval start, struct timeval end) {
    return ((end.tv_sec - start.tv_sec) + ((end.tv_usec - start.tv_usec) / 1e6));
}

int saveMatrixInFile(const char *fileName, double *matrix, int matrix_rows, int matrix_cols) {
    // Open the file
    FILE *filePtr;
    filePtr = fopen(fileName, "w+");
    // Check possible errors
    if (filePtr == NULL) {
        perror("output file error");
        return 1;
    }
    int numArrays = 1;
    fprintf(filePtr, "%d\n", numArrays);
    fprintf(filePtr, "%d\n", matrix_rows * matrix_cols);
    for (int i = 0; i < matrix_rows; i++) {
        for (int j = 0; j < matrix_cols; j++) {
            fprintf(filePtr, "%.15f\n", matrix[j * matrix_rows + i]);
        }
    }
    fclose(filePtr);
    return 0;
}

int initInputVectors(const char *fileName, double *u10_time, double *u20_time, double *w0_time, int dimension) {
    FILE *filePtr;
    int readedDimension;
    // Open the file
    filePtr = fopen(fileName, "r+");
    // Check possible errors
    if (filePtr == NULL) {
        perror("input file error");
        return 1;
    }
    // Check if the file is empty
    fseek(filePtr, 0, SEEK_END);
    long size = ftell(filePtr);
    if (size == 0) {
        fprintf(stdout, "\nFile %s is empty. Please provide to fill opportunely the file with values.\n", fileName);
        return 1;
    }
    fseek(filePtr, 0, SEEK_SET);
    // Vector u10_time
    fscanf(filePtr, "%d", &readedDimension);
    for (int i = 0; i < dimension; i++) {
        fscanf(filePtr, "%lf", (u10_time + i));
    }
    // Vector u10_time
    fscanf(filePtr, "%d", &readedDimension);
    for (int i = 0; i < dimension; i++) {
        fscanf(filePtr, "%lf", (u20_time + i));
    }
    // Vector u10_time
    fscanf(filePtr, "%d", &readedDimension);
    for (int i = 0; i < dimension; i++) {
        fscanf(filePtr, "%lf", (w0_time + i));
    }
    fclose(filePtr);
    return 0;
}

int saveVectorsInFile(const char *fileName, int elements, double *arr1, int dim1, ...) {
    va_list ptr;
    FILE *filePtr;
    // Check the number of elements to write on the file
    if (elements == 0) {
        return 1;
    }

    // Initializing argument to the
    // list pointer
    va_start(ptr, dim1);

    // Open the file
    filePtr = fopen(fileName, "w+");
    // Check possible errors
    if (filePtr == NULL) {
        perror("File opening error");
        return 1;
    }

    // Write the number of arrays
    fprintf(filePtr, "%d\n", elements);

    // Write the number of elements of the array
    fprintf(filePtr, "%d\n", dim1);
    // Save the array in the file
    for (int i = 0; i < dim1; i++) {
        fprintf(filePtr, "%.15lf\n", arr1[i]);
    }

    // Write the the other vectors
    for (int i = 0; i < elements - 1; i++) {
        // Accessing current variable
        // and pointing to next one
        double *tempArray = va_arg(ptr, double *);
        int tempDim = va_arg(ptr, int);
        
        // Write the number of elements of the array
        fprintf(filePtr, "%d\n", tempDim);
        // Save the array in the file
        for (int i = 0; i < tempDim; i++) {
            fprintf(filePtr, "%.15lf\n", tempArray[i]);
        }
    }
 
    // Ending argument list traversal
    va_end(ptr);
    // Close the file
    fclose(filePtr);
    return 0;
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