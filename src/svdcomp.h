#ifndef SVDCOMP_H
#define SVDCOMP_H

#include "mastor.h"

void svdcomp(double **a, int m, int n, double w[], double **v);

void svdcomp_u_and_v(double **a, int m, int n, double w[], double **u, double **v);

#endif
