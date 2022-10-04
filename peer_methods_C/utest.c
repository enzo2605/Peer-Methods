#include <math.h>
#include <cblas.h>
#include <string.h>
#include "CLab.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    int rip = 1;
    while (rip == 1) {
        int M = 8;
        double *a = (double *)calloc(M, sizeof(double));
        for (int i = 0; i < M; i++) {
            *(a + i) = i + 1;
        }
        double *b = (double *)calloc(M, sizeof(double));
        for (int i = 0; i < M; i++) {
            *(b + i) = i + 1;
        }

        double res = cblas_ddot(M, a, 1, b, 1);
        fprintf(stdout, "res: %lf\n", res);
    }
    exit(0);
}