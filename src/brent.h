#ifndef BRENT_H
#define BRENT_H

double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);

#endif
