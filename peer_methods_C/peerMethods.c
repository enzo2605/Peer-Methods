#include "peerMethods.h"

void initReturnStruct(return_values *rv) {
    rv->t = NULL;
    rv->y = NULL;
    rv->yT = NULL;
}

void saveInFile(const char* fileName, return_values result) {
    // Open the file
    FILE *filePtr;
    filePtr = fopen(fileName, "w");
    // Check possible errors
    if (filePtr == NULL) {
        fprintf(stdout, "\nError while opening %s.\n", fileName);
        exit(1);
    }
    int numArrays = 3;
    fprintf(filePtr, "%d\n", numArrays);
    // Write y_T
    fprintf(filePtr, "%d\n", result.yT_size);
    for (int i = 0; i < result.yT_size; i++) {
        fprintf(filePtr, "%.4f\n", result.yT[i]);
    }
    // Write y
    fprintf(filePtr, "%d\n", result.y_rows * result.y_cols);
    for (int i = 0; i < result.y_rows; i++) {
        for (int j = 0; j < result.y_cols; j++) {
            fprintf(filePtr, "%.4f\n", result.y[i * result.y_cols + j]);
        }
    }
    // Write t
    // Write y_T
    fprintf(filePtr, "%d\n", result.t_size);
    for (int i = 0; i < result.t_size; i++) {
        fprintf(filePtr, "%.4f\n", result.t[i]);
    }
    // Close the file
    fclose(filePtr);
}

double_t *Sherratt(double_t *y0, int y0Size, double_t *L, int Lsize, int *sherrattSize) {
    // Declaring three dynamic array representing the three function contained in y0
    double_t *U = (double_t *)Calloc(M, sizeof(double_t));
    double_t *V = (double_t *)Calloc(M, sizeof(double_t));
    double_t *W = (double_t *)Calloc(M, sizeof(double_t));
    for (int i = 0; i < M; i++) {
        U[i] = y0[i];
        V[i] = y0[i + M];
        W[i] = y0[i + 2 * M];
    }

    int newSize = M * 3;
    double_t *fun = (double_t *)Calloc(newSize, sizeof(double_t));
    for (int i = 0; i < M; i++) {
        fun[i] = W[i] * U[i] * (U[i] + H * V[i]) - B1 * U[i] - S * U[i] * V[i];
        fun[i + M] = F * W[i] * V[i] * (U[i] + H * V[i]) - B2 * V[i];
        fun[i + 2 * M] = a - W[i] - W[i] * (U[i] + V[i]) * (U[i] + H * V[i]);
    }

    // Matrix by vector product
    double_t *Ly = (double_t *)Calloc(newSize, sizeof(double_t));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, Lsize, Lsize, 1, L, Lsize, y0, 1, 0, Ly, 1);

    // Punctual sum of two vectors
    double_t *dydt = sumPuntVectors(Ly, fun, Lsize);
    *sherrattSize = Lsize;

    // Free all the useless memory allocated
    freeEverything(U, V, W, fun, Ly, (void *)0);

    return dydt;
}

double_t *RungeKutta4th(double_t h, double_t t0, double_t *y0, int y0Size, double_t *L, int Lsize, int *ySize) {
    double_t c[4] = { 0.0f, 1.0f / 2.0f, 1.0f / 2.0f, 1.0f };
    double_t A[4][4] = { 
                        { 0.0f, 0.0f, 0.0f, 0.0f }, 
                        { 1.0f / 2.0f, 0.0f, 0.0f, 0.0f }, 
                        { 0.0f, 1.0f / 2.0f, 0.0f, 0.0f }, 
                        { 0.0f, 0.0f, 1.0f, 0.0f }
                    };
    double_t b[4] = { 1.0f / 6.0f, 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 6.0f };

    // Compute Y1
    double_t *Y1 = y0;
    //printDVector(Y1, y0Size, "Y1");

    // Compute Y2
    int fY1_size;
    double_t *fY1 = Sherratt(Y1, y0Size, L, Lsize, &fY1_size);
    scalarByVector(fY1, y0Size, h * A[1][0]);
    double_t *Y2 = sumPuntVectors(y0, fY1, y0Size);
    //printDVector(Y2, y0Size, "Y2");

    // Compute Y3
    int fY2_size;
    double_t *fY2 = Sherratt(Y2, y0Size, L, Lsize, &fY2_size);
    scalarByVector(fY2, fY2_size, h * A[2][1]);
    double_t *Y3 = sumPuntVectors(y0, fY2, fY2_size);
    //printDVector(Y3, y0Size, "Y3");

    // Compute Y4
    int fY3_size;
    double_t *fY3 = Sherratt(Y3, y0Size, L, Lsize, &fY3_size);
    scalarByVector(fY3, fY3_size, h * A[3][2]);
    double_t *Y4 = sumPuntVectors(y0, fY3, fY3_size);
    //printDVector(Y4, y0Size, "Y4");

    // Compute y
    int fY4_size;
    double_t *fY4;
    fY1 = Sherratt(Y1, y0Size, L, Lsize, &fY1_size);
    fY2 = Sherratt(Y2, y0Size, L, Lsize, &fY2_size);
    fY3 = Sherratt(Y3, y0Size, L, Lsize, &fY3_size);
    fY4 = Sherratt(Y4, y0Size, L, Lsize, &fY4_size);

    scalarByVector(fY1, fY1_size, h * b[0]);
    scalarByVector(fY2, fY2_size, h * b[1]);
    scalarByVector(fY3, fY3_size, h * b[2]);
    scalarByVector(fY4, fY4_size, h * b[3]);

    double_t *y = (double_t *)Calloc(y0Size, sizeof(double_t));
    y = sumPuntVectors(fY1, fY2, y0Size);
    y = sumPuntVectors(fY3, y, y0Size);
    y = sumPuntVectors(fY4, y, y0Size);
    y = sumPuntVectors(y, y0, y0Size);
    *ySize = y0Size;
    //printDVector(y, *ySize, "y");

    return y;
}

void fPeerClassic_twoStages(int N, double_t *t_span, int t_span_size, double_t *L, int Lsize, double_t *y0, int y0_size, return_values *collect_result) {
    /******************************* 
     * Fixing method coefficients 
     * ****************************/
    double_t b11 = 0.0f, b21 = 0.0f, b12 = 1.0f, b22 = 1.0f, c1 = 0.0f, c2 = 1.0f, r21 = 0.0f;

    double_t a11 = -((b11 - 2 * b11 * c1 - pow(c1, 2) + b11 * pow(c1, 2)) / (2 * (-1 + c1)));
    double_t a12 = -((b11 + 2 * c1 - 2 * b11 * c1 - pow(c1, 2) + b11 * pow(c1, 2)) / (2 * (-1 + c1)));
    double_t a21 = -((-1 + b21 - 2 * b21 * c1 + b21 * pow(c1, 2) + 2 * c1 * r21) / (2 * (-1 + c1)));
    double_t a22 = -((3 + b21 - 2 * c1 - 2 * b21 * c1 + b21 * pow(c1, 2) - 2 * r21) / (2 * (-1 + c1)));
    //fprintf(stdout, "a11: %f\na12: %f\na21: %f\na22: %f\n", a11, a12, a21, a22);

    double_t c[STAGES] = { c1, c2 };
    //printDVector(c, STAGES, "c");

    double_t *A = (double_t *)Calloc(STAGES * STAGES, sizeof(double_t));
    double_t tempA[STAGES * STAGES] = { a11, a12, a21, a22 };
    initMatrixByRowWithValuesFromVector(A, STAGES, STAGES, tempA, STAGES * STAGES);
    //printDMatrix(A, STAGES, STAGES, "A");

    double_t *B = (double_t *)Calloc(STAGES * STAGES, sizeof(double_t));
    double_t tempB[STAGES * STAGES] = { b11, b12, b21, b22 };
    initMatrixByRowWithValuesFromVector(B, STAGES, STAGES, tempB, STAGES * STAGES);
    //printDMatrix(B, STAGES, STAGES, "B");

    double_t *R = (double_t *)Calloc(STAGES * STAGES, sizeof(double_t));
    double_t tempR[STAGES * STAGES] = { 0.0f, 0.0f, r21, 0.0f };
    initMatrixByRowWithValuesFromVector(R, STAGES, STAGES, tempR, STAGES * STAGES);
    //printDMatrix(R, STAGES, STAGES, "R");

    /******************************* 
     *  Compute the solution
     * ****************************/
    double_t h = (t_span[1] - t_span[0]) / N;
    double_t *t = linspace(t_span[0], t_span[1], N + 1);
    int t_size = N + 1;
    int n = 1;

    int s = STAGES; // Number of stages
    int d1 = y0_size; // Dimension of the problem
    int Y_rows = s * d1, Y_cols = N;
    double_t *Y = zerosMatrixD(Y_rows, Y_cols);
    //fprintf(stdout, "\nh: %lf\ns: %d\nd1: %d\n", h, s, d1);
    //printDVector(t, N + 1, "t");
    //printDMatrix(Y, Y_rows, Y_cols, "Y");

    /************************************************
     * Runge-Kutta of order four to initialize stages 
     * **********************************************/
    double_t *FYiRK;
    int FYiRK_size;
    for (int i = 0; i < s; i++) {
        FYiRK = RungeKutta4th(c[i] * h, t[0], y0, y0_size, L, Lsize, &FYiRK_size);
        for (int k = 0; k < d1; k++) {
            Y[(n - 1) * Y_rows + (i * d1 + k)] = FYiRK[k];
        }
    }
    //printDMatrix(Y, Y_rows, Y_cols, "Y");
    //printDVector(FYiRK, FYiRK_size, "FYiRK");

    int y_rows = d1, y_cols = N + 1;
    double_t *y = zerosMatrixD(y_rows, y_cols);

    for (int k = 0; k < d1; k++) {
        y[0 * y_rows + k] = y0[k];
    }

    for (int k = 0; k < d1; k++) {
        y[(n) * y_rows + k] = Y[(n - 1) * Y_rows + ((s - 1) * d1 + k)];
    }
    //printDMatrix(y, y_rows, y_cols, "y");

    int Fnm1_size = s * d1;
    double_t *Fnm1 = zerosD(Fnm1_size);
    double_t *Yi = zerosD(d1);

    for (int i = 0; i < s; i++) {
        for (int k = 0; k < d1; k++) {
            Yi[k] = Y[(n - 1) * Y_rows + ((i) * d1 + k)];
        }
        int FYi_size;
        double_t *FYi = Sherratt(Yi, d1, L, Lsize, &FYi_size);
        for (int k = 0; k < d1; k++) {
            Fnm1[(i) * d1 + k] = FYi[k];
        }
    }
    //printDVector(Yi, d1, "Yi");
    //printDVector(Fnm1, Fnm1_size, "Fnm1");


    for (n = 1; n < N; n++) {
        //fprintf(stdout, "\nn: %d \n", n);
        for (int i = 0; i < s; i++) {
            for (int k = 0; k < d1; k++) {
                for (int j = 0; j < s; j++) {
                    Y[(n) * Y_rows + ((i) * d1 + k)] += B[j * STAGES + i] * Y[(n - 1) * Y_rows + ((j) * d1 + k)] + h * A[j * STAGES + i] * Fnm1[j * d1 + k];
                }
            }
        }
        //printDMatrix(Y, Y_rows, Y_cols, "Y");

        Fnm1 = zerosD(Fnm1_size);
        Yi = zerosD(d1);
        for (int i = 0; i < s; i++) {
            for (int k = 0; k < d1; k++) {
                Yi[k] = Y[(n) * Y_rows + ((i) * d1 + k)];
            }
            int FYi_size;
            double_t *FYi = Sherratt(Yi, d1, L, Lsize, &FYi_size);
            for (int k = 0; k < d1; k++) {
                Fnm1[(i) * d1 + k] = FYi[k];
            }
        }
        //printDVector(Yi, d1, "Yi");
        //printDVector(Fnm1, Fnm1_size, "Fnm1");

        for (int k = 0; k < d1; k++) {
            y[(n + 1) * y_rows + k] = Y[(n) * Y_rows + ((s - 1) * d1 + k)];
        }
    }
    //printDMatrix(Y, Y_rows, Y_cols, "Y");
    //printDMatrix(y, y_rows, y_cols, "y");
    //exit(0);

    //fprintf(stdout, "Here\n");
    double_t *yT = zerosD(d1);
    int yT_size = d1;
    for (int i = 0; i < d1; i++) {
        //fprintf(stdout, "y[(N) * y_cols + i]: %f\n", y[(N) * y_rows + i]);
        yT[i] = y[(N) * y_rows + i];
    }
    //printDVector(yT, yT_size, "yT");

    // After all the calculation, collecting and return the results
    collect_result->y = y;
    collect_result->y_rows = y_rows;
    collect_result->y_cols = y_cols;

    collect_result->t = t;
    collect_result->t_size = t_size;

    collect_result->yT = yT;
    collect_result->yT_size = yT_size;
}