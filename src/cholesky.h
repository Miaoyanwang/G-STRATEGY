#ifndef CHOLESKY_H
#define CHOLESKY_H

#include <stdio.h>
#include <math.h>

int cholesky(double **orig, int n, double **aug, int mcol,double **chol, double **cholaug, int ofs);

#endif
