#include "peerMethods.h"

double *Sherratt(double *y0, int y0Size, double *L, int Lsize, int *sherrattSize) {
    // Declaring three dynamic array representing the three function contained in y0
    double *U = (double *)Calloc(M, sizeof(double));
    double *V = (double *)Calloc(M, sizeof(double));
    double *W = (double *)Calloc(M, sizeof(double));
    for (int i = 0; i < M; i++) {
        U[i] = y0[i];
        V[i] = y0[i + M];
        W[i] = y0[i + 2 * M];
    }

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
    *sherrattSize = Lsize;

    // Free all the useless memory allocated
    freeEverything(U, V, W, fun, Ly);

    return dydt;
}

/*
function y = RK4(funz,h,t0,y0)

c = [0;1/2;1/2;1];
A = [0 0 0 0;1/2 0 0 0;0 1/2 0 0;0 0 1 0];
b = [1/6,1/3,1/3,1/6];

Y1 = y0;
Y2 = y0 + h*A(2,1)*funz(t0,Y1);
Y3 = y0 + h*A(3,2)*funz(t0+c(2)*h,Y2);
Y4 = y0 + h*A(4,3)*funz(t0+c(3)*h,Y3);

y = y0 + h*(b(1)*funz(t0,Y1) + b(2)*funz(t0+c(2)*h,Y2) + ...
    b(3)*funz(t0+c(3)*h,Y3) + b(4)*funz(t0+c(4)*h,Y4));

end
*/
void RungeKutta4th(double h, double t0, double *y0, int y0Size, double *y, int *ySize) {
    
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

    %% Compute solution

    h = (tspan(2)-tspan(1))/N;
    t = linspace(tspan(1),tspan(2),N+1); 
    n = 1;

    s = length(c); % Number of stages
    d = length(y0); % Dimension of the problem
    Y = zeros(s*d,N);

    % Runge-Kutta of order four to initialize stages

    for i = 1:s
        FYiRK = RK4(funz,c(i)*h,t(1),y0);
        for k = 1:d
            Y((i-1)*d+k,n) = FYiRK(k);
        end
    end

    y = zeros(d,N+1);
    for k = 1:d
        y(k,1) = y0(k);
    end 
*/
void fPeerClassic_twoStages(int N, double *t_span, int t_span_size, double *L, int Lsize, double *y0, int y0_size, double *yT_ClPeer, int *yT_ClPeer_rows, int *yT_ClPeer_cols, double *y_ClPeer, int *y_ClPeer_size, double *t,  int *t_size) {
    /******************************* 
     * Fixing method coefficients 
     * ****************************/
    double b11 = 0.0f, b21 = 0.0f, b12 = 1.0f, b22 = 1.0f, c1 = 0.0f, c2 = 1.0f, r21 = 0.0f;

    double a11 = -((b11 - 2 * b11 * c1 - pow(c1, 2) + b11 * pow(c1, 2)) / (2 * (-1 + c1)));
    double a12 = -((b11 + 2 * c1 - 2 * b11 * c1 - pow(c1, 2) + b11 * pow(c1, 2)) / (2 * (-1 + c1)));
    double a21 = -((-1 + b21 - 2 * b21 * c1 + b21 * pow(c1, 2) + 2 * c1 * r21) / (2 * (-1 + c1)));
    double a22 = -((3 + b21 - 2 * c1 - 2 * b21 * c1 + b21 * pow(c1, 2) - 2 * r21) / (2 * (-1 + c1)));
    fprintf(stdout, "a11: %f\na12: %f\na21: %f\na22: %f\n", a11, a12, a21, a22);

    double c[STAGES] = { c1, c2 };
    printDVector(c, STAGES, "c");

    double *A = (double *)Calloc(STAGES * STAGES, sizeof(double));
    double tempA[STAGES * STAGES] = { a11, a12, a21, a22 };
    initMatrixByRowWithValuesFromVector(A, STAGES, STAGES, tempA, STAGES * STAGES);
    printDMatrix(A, STAGES, STAGES, "A");

    double *B = (double *)Calloc(STAGES * STAGES, sizeof(double));
    double tempB[STAGES * STAGES] = { b11, b12, b21, b22 };
    initMatrixByRowWithValuesFromVector(B, STAGES, STAGES, tempB, STAGES * STAGES);
    printDMatrix(B, STAGES, STAGES, "B");

    double *R = (double *)Calloc(STAGES * STAGES, sizeof(double));
    double tempR[STAGES * STAGES] = { 0.0f, 0.0f, r21, 0.0f };
    initMatrixByRowWithValuesFromVector(R, STAGES, STAGES, tempR, STAGES * STAGES);
    printDMatrix(R, STAGES, STAGES, "R");

    /******************************* 
     *  Compute the solution
     * ****************************/
    double h = (t_span[2] - t_span[1]) / N;
    t = linspace(t_span[0], t_span[1], N + 1);
    double n = 1.0f;

    int s = STAGES; // Number of stages
    int d1 = y0_size; // Dimension of the problem
    int Y_rows = s * d1, Y_cols = N;
    double *Y = zerosMatrixD(Y_rows, Y_cols);

    /************************************************
     * Runge-Kutta of order four to initialize stages 
     * **********************************************/


}