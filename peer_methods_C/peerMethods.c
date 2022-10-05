#include "peerMethods.h"

double *Sherratt(double *y0, double *U, double *V, double *W, int m, double *L, int Lsize) {
    int newSize = M * 3;
    double *fun = (double *)Calloc(newSize, sizeof(double));

    for (int i = 0; i < M; i++) {
        fun[i] = W[i] * U[i] * (U[i] + H * V[i]) - B1 * U[i] - S * U[i] * V[i];
        fun[i + M] = F * W[i] * V[i] * (U[i] + H * V[i]) - B2 * V[i];
        fun[i + 2 * M] = a - W[i] - W[i] * (U[i] + V[i]) * (U[i] + H * V[i]);
    }

    // Matrix by vector product
    double *Ly = (double *)Calloc(newSize, sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, Lsize, Lsize, 1, L, Lsize, y0, 1, 0, Ly, 1);
    // Sum of two matrices of the same size
    double *dydt = sumPuntSquareMatrices(Ly, fun, Lsize);

    return dydt;
}

/*
        
    b11 = 0; b21 = 0; b12 = 1; b22 = 1;
    c1 = 0; c2 = 1;
    r21 = 0;
    
    a11 = -((b11 - 2 * b11 * c1 - c1^2 + b11 * c1^2)/(2 * (-1 + c1)));
    a12 = -((b11 + 2 * c1 - 2 * b11 * c1 - c1^2 + b11 * c1^2)/(2 * (-1 + c1)));
    a21 = -((-1 + b21 - 2 * b21 * c1 + b21 * c1^2 + 2 * c1 * r21)/(2 * (-1 + c1)));
    a22 = -((3 + b21 - 2 * c1 - 2 * b21 * c1 + b21 * c1^2 - 2 * r21)/(2 * (-1 + c1)));
    
    c = [c1 c2];
    A = [a11 a12; a21 a22];
    B = [b11 b12; b21 b22];
    R = [0 0; r21 0];
*/
void fPeerClassic_twoStages(int N, double *t_span, int t_span_size, double *y0, int y0_size, double *yT_ClPeer, int yT_ClPeer_rows, int yT_ClPeer_cols, double *y_ClPeer, int y_ClPeer_size, double *t,  int t_size) {
    double b11 = 0.0f, b21 = 0.0f, b12 = 1.0f, b22 = 1.0f, c1 = 0.0f, c2 = 1.0f, r21 = 0.0f;

    double a11 = -((b11 - 2 * b11 * c1 - pow(c1, 2) + b11 * pow(c1, 2)) / (2 * (-1 + c1)));
    double a12 = -((b11 + 2 * c1 - 2 * b11 * c1 - pow(c1, 2) + b11 * pow(c1, 2)) / (2 * (-1 + c1)));
    double a21 = -((-1 + b21 - 2 * b21 * c1 + b21 * pow(c1, 2) + 2 * c1 * r21) / (2 * (-1 + c1)));
    double a22 = -((3 + b21 - 2 * c1 - 2 * b21 * c1 + b21 * pow(c1, 2) - 2 * r21) / (2 * (-1 + c1)));

    double c[STAGES] = { c1, c2 };

    double *A = (double *)Calloc(STAGES * STAGES, sizeof(double));
    double tempA[4] = { a11, a12, a21, a22 };
    initMatrixByRowWithValuesFromVector(A, STAGES, STAGES, tempA, STAGES * STAGES);

    double *B = (double *)Calloc(STAGES * STAGES, sizeof(double));
    double tempB[4] = { b11, b12, b21, b22 };
    initMatrixByRowWithValuesFromVector(B, STAGES, STAGES, tempB, STAGES * STAGES);

    double *R = (double *)Calloc(STAGES * STAGES, sizeof(double));
    double temp[4] = { 0.0f, 0.0f, r21, 0.0f };
    initMatrixByRowWithValuesFromVector(A, STAGES, STAGES, temp, STAGES * STAGES);
}