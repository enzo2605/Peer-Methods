/**
 * Title: Peer methods
 * Author: Vincenzo Iannucci
 * **/
#include <math.h>
#include <cblas.h>
#include "CLab.h"
#include "utilities.h"

#define a  1.5
#define B1 0.45
#define B2 0.3611 
#define F  0.802 
#define H  0.802 
#define S  0.0002 
#define d  500 
#define D  0.802

int main(int argc, char *argv[]) {
    // intervals
    double t_start, t_end, x_start, x_end;

    // Checking arguments passed by the user
    if (argc != 5) {
        printf("Usage: %s t_start t_end x_start x_end\n", argv[0]);
        printf("Using default parameters...\n\n");
        t_start = 0.0f;
        t_end   = 50.0f;
        x_start = -50.0f;
        x_end   = 50.0f;
    }
    else {
        t_start = atof(argv[1]);
        t_end   = atof(argv[2]);
        x_start = atof(argv[3]);
        x_end   = atof(argv[4]);
    }

    // Printing the interval used
    printf("Interval used:\ntime span: [%f, %f]\nspace span: [%f, %f]\n", t_start, t_end, x_start, x_end);

    // Time initialization
    double t_span[2] = { t_start, t_end };
    double Delta_t = 1.0f / pow(2.0f, 11.0f);

    double *t_int;
    int n_points_t;
    t_int = setDVector(t_int, t_start, t_end, Delta_t, &n_points_t);
    int N = (t_span[1] - t_span[0]) / Delta_t;

    // Space initialization
    int M = 64;
    double x_span[2] = { x_start, x_end };
    double Delta_x = (x_span[1] - x_span[0]) / M;

    double *x_int;
    int n_points_x;
    x_int = setDVector(x_int, x_start, x_end, Delta_x, &n_points_x);

    

    exit(0);
}