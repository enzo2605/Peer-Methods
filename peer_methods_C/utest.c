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
        initializeMatrix(a, size, size);
        printDMatrix(a, size, size);
        double *b = (double *)calloc(size * size, sizeof(double));
        initializeMatrix(b, size, size);
        printDMatrix(b, size, size);
        double *c = (double *)calloc(size * size, sizeof(double));
        initializeMatrix(c, size, size);
        printDMatrix(c, size, size);

        double *pack = packMatrices(size, a, b, c);
        printDMatrix(pack, size, size * 3);

        double *blockMatrix = threeBlockDiagD(size, a, b, c);
        printDMatrix(blockMatrix, size * 3, size * 3);
    }
    exit(0);
}