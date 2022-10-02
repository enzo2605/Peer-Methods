#include <math.h>
#include <cblas.h>
#include <string.h>
#include "CLab.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    int rip = 1;
    while (rip == 1) {
        int size, k, N_mat;
        printf("\nInserisci n: ");
        scanf("%d", &size);

        double *a = (double *)calloc(size * size, sizeof(double));
        initializeRandomMatrix(a, size, size);
        printDMatrix(a, size, size);
        double *b = (double *)calloc(size * size, sizeof(double));
        initializeRandomMatrix(b, size, size);
        printDMatrix(b, size, size);
        double *c = sumPuntSquareMatrices(a, b, size);
        printDMatrix(c, size, size);
    }
    exit(0);
}