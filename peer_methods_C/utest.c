#include <math.h>
#include <cblas.h>
#include <string.h>
#include "CLab.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    srand((unsigned int)time(NULL));
    int rip = 1;
    while (rip == 1) {
        int size, k, N_mat;
        printf("\nInserisci n: ");
        scanf("%d", &size);

        double *a = (double *)calloc(size, sizeof(double));
        initializeRandomVector(a, size);
        printDVector(a, size);
        double *b = (double *)calloc(size, sizeof(double));
        initializeRandomVector(b, size);
        printDVector(b, size);
        double *c = (double *)calloc(size, sizeof(double));
        initializeRandomVector(c, size);
        printDVector(c, size);

        int newDimension = 0;
        double *pack = packThreeVectors(size, a, b, c, &newDimension);
        printDVector(pack, newDimension);
    }
    exit(0);
}