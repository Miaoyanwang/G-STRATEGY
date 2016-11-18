#include "brent.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define MAX_ITER 10000

/* brent's algorithm, modified from GSL min/brent.c
*/


double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin) {
  /* state initialization */

  const double golden = 0.3819660;      /* golden = (3 - sqrt(5))/2 */
  double x_lower = (ax < cx) ? ax : cx;
  double f_lower = (*f)(x_lower);
  double x_minimum = bx;
  double f_minimum = (*f)(bx);
  double x_upper = (ax < cx) ? cx : ax;
  double f_upper = (*f)(x_upper);
  double z = x_minimum;
  double f_z = f_minimum;

  double tolerance = GSL_SQRT_DBL_EPSILON * fabs (z);

  double v = x_lower + golden * (x_upper - x_lower);
  double w = v;


  double f_v = (*f)(v);
  double f_w = f_v;

  double d = 0;
  double e = 0;
  int first = 1;

  /* loop variables */
  double tmp;
  double w_lower, w_upper,midpoint;
  double x_left, x_right;
  double u,f_u;
  double p,q,r;
  double t2;

  /* iteration of brent  */
  int i = 0;
  while (first != 0 || (fabs(e) > (tol* fabs(z) + GSL_DBL_EPSILON))) {
      
    z = x_minimum;
    f_z = f_minimum;
    i++;
    /*printf("iteration number %d: e=%.9f, tol=%.9f\n",i,fabs(e),tolerance); 
      printf("a,b,c (%.6f,%.6f,%.6f) -> fa,fb,fc (%.6f,%.6f,%.6f)\n",x_lower,x_minimum,x_upper,f_lower,f_minimum,f_upper);*/
    if (first != 0) {
      first = 0;
    }
    if (i > MAX_ITER) {
      break;
    }
    x_left = x_lower;
    x_right = x_upper;

    /* swap e and d */
    tmp = d;
    d = e;
    e = tmp;
  
    w_lower = (z - x_left);
    w_upper = (x_right - z);

    tolerance =  GSL_SQRT_DBL_EPSILON * fabs (z);

    p = 0, q = 0, r = 0;

    midpoint = 0.5 * (x_left + x_right);

    if (fabs(z-midpoint) <= 2*(tol*fabs(z) + GSL_DBL_EPSILON)) {
      break;
    }

    if (fabs (e) > tolerance) {
      /* fit parabola */
      r = (z - w) * (f_z - f_v);
      q = (z - v) * (f_z - f_w);
      p = (z - v) * q - (z - w) * r;
      q = 2 * (q - r);
      
      if (q > 0) {
	p = -p;
      } else {
	q = -q;
      }

      r = e;
      e = d;
    }

    if (fabs (p) < fabs (0.5 * q * r) && p < q * w_lower && p < q * w_upper) {
      t2 = 2 * tolerance ;
      
      d = p / q;
      u = z + d;

      if ((u - x_left) < t2 || (x_right - u) < t2) {
	d = (z < midpoint) ? tolerance : -tolerance ;
      }
    } else {
      e = (z < midpoint) ? x_right - z : -(z - x_left) ;
      d = golden * e;
    }

    
    if (fabs (d) >= tolerance) {
      u = z + d;
    } else {
      u = z + ((d > 0) ? tolerance : -tolerance) ;
    }

    f_u = (*f)(u);

    if (f_u <= f_z) {
      /* case a */
      if (u < z) {
	x_upper = z;
	f_upper = f_z;
      } else {
	x_lower = z;
	f_lower = f_z;
      }

      v = w;
      f_v = f_w;
      w = z;
      f_w = f_z;
      x_minimum = u;
      f_minimum = f_u;
    } else {
      /* case b */
      if (u < z) {
	x_lower = u;
	f_lower = f_u;
      } else {
	x_upper = u;
	f_upper = f_u;
      }

      if (f_u <= f_w || w == z) {
	/* case b.1 */
	v = w;
	f_v = f_w;
	w = u;
	f_w = f_u;
      } else if (f_u <= f_v || v == z || v == w) {
	/* case b.2 */
	v = u;
	f_v = f_u;
      }
    }
  }

  *xmin = x_minimum;
    
    return(0);
}


