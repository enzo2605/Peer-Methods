/**
 * Author: Vincenzo Iannucci
 * Purpose: The library provides an implementation for the main function
 * for solving peer method.
 * **/
#ifndef peerMethods_h
#define peerMethods_h

#include "utilities.h"
#include <cblas.h>
#include "CLab.h"

#define a  1.5
#define B1 0.45
#define B2 0.3611 
#define F  0.802 
#define H  0.802 
#define S  0.0002 
#define d  500 
#define D  0.802
#define M  4
#define STAGES 2

double *Sherratt(double *y0, double *U, double *V, double *W, int m, double *L, int Lsize);
void fPeerClassic(int N, double *t_span, int t_span_size, double *y0, int y0_size, double *yT_ClPeer, int yT_ClPeer_rows, int yT_ClPeer_cols, double *y_ClPeer, int y_ClPeer_size, double *t,  int t_size);

#endif // !peerMethods_h