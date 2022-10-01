#include <math.h>
#include <cblas.h>
#include <string.h>
#include "CLab.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    int rip = 1;
    while (rip == 1) {
        int size, k, N_mat;
        printf("\nInserisci n e k: ");
        scanf("%d %d", &size, &k);
        double *v = (double *)Calloc(size, sizeof(double));
        for (int i = 0; i < size; i++) {
            v[i] = i + 1.0f;
        }
        printDVector(v, size);
        double *mat = diagD(v, size, k, &N_mat);
        printDMatrix(mat, N_mat, N_mat);
    }
    exit(0);
}