#include <math.h>
#include <cblas.h>
#include <string.h>
#include "CLab.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
    int rip = 1;
    while (rip == 1) {
        double x1, x2;
        int n;
        fscanf(stdin, "%lf %lf %d", &x1, &x2, &n);
        double *res = linspace(x1, x2, n);
        printDVector(res, n, "res");
    }
    exit(0);
}