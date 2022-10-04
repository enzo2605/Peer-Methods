#include "peerMethods.h"

/**
 * function dydt = Sherratt(t,y)

    global a B1 B2 F H S d D L M

    %% Case of finite differences of order four and periodic boundary conditions

    % Separe the components of y in U, V and W

    for i = 1:M
        U(i) = y(i);
        V(i) = y(i+M);
        W(i) = y(i+2*M);
    end

    % Define Fun

    for i = 1:M
        Fun1(i) = W(i)*U(i)*(U(i)+H*V(i)) - B1*U(i) - S*U(i)*V(i);
        Fun2(i) = F*W(i)*V(i)*(U(i)+H*V(i)) - B2*V(i);
        Fun3(i) = a - W(i) - W(i)*(U(i)+V(i))*(U(i)+H*V(i));
    end
    
    Fun = [Fun1,Fun2,Fun3]';

    % Determine y(t_{n+1}) at all space grid points

    % To improve using the structure of L
    dydt = L*y + Fun;
end
*/

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