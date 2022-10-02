#include <math.h>
#include <cblas.h>
#include <string.h>
#include "CLab.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    int rip = 1;
    while (rip == 1) {
        int M, k, N_mat;
        fprintf(stdout, "\nInserisci M: ");
        scanf("%d", &M);

        double *onesVector, *tempDiagOne, *tempDiagMinusOne; 
        int sizeTempDiagOne, sizeTempDiagMinusOne;
        onesVector = onesD(onesVector, M - 1);
        printDVector(onesVector, M - 1);
        fprintf(stdout, "hi\n");
        tempDiagOne      = diagD(onesVector, M - 1, 1, &sizeTempDiagOne);
        fprintf(stdout, "hi\n");
        printDMatrix(tempDiagOne, sizeTempDiagOne, sizeTempDiagOne);
    }
    exit(0);
}