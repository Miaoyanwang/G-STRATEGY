#include "brak.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define EVAL_MAX 10000

#include "nrutil.h"

/* initial_ax finds a starting value for alpha */
void initial_ax(struct FAMILY *family, struct MLE_PARAM *mle_initial, struct DATA_STRUCT data_struct) {
	
	int i,j,k,n=0,n_fam,unr=0;
	n_fam = data_struct.n_fam;
	
	double cov_dot_cov=0.0, cov_dot_trait=0.0, trait_adj_dot_trait_adj=0.0, outer_dot=0.0;
	double beta_tilda, err_var_tilda, add_var_tilda;
	
	int n_pheno=0;
	int *pheno_typed;
	double **cov;
	double *trait;
	double **phi;
	
	/* pull out variables from the struct family to use to calulate W^tW and W^tD
	 (i.e. we are assuming K = I and doing simple ordinary regression
	 w/o accounting for the dependence due to relatedness)
	 */
	for(i=1; i<= n_fam; i++) {
		n_pheno = family[i].n_pheno;
		cov = family[i].cov;
		trait = family[i].trait;
		
		for(j=1; j<=n_pheno; j++) {
			cov_dot_cov = cov_dot_cov + cov[j][1]*cov[j][1];
			cov_dot_trait = cov_dot_trait + cov[j][1]*trait[j];
		}
	}
	
	beta_tilda = cov_dot_trait/cov_dot_cov;
	
	/* estimate the error variance ignoring relatedness (D-W*beta_tilda)^t * (D-W*beta_tilda) */
	for(i=1; i<= n_fam; i++) {
		n_pheno = family[i].n_pheno;
		pheno_typed = family[i].pheno_typed;
		cov = family[i].cov;
		trait = family[i].trait;
		
		n = n + n_pheno;
		
		for(j=1; j<=n_pheno; j++) {
			trait_adj_dot_trait_adj = trait_adj_dot_trait_adj + (trait[j]-cov[j][1]*beta_tilda)*(trait[j]-cov[j][1]*beta_tilda);
		}
	}
	
	err_var_tilda = trait_adj_dot_trait_adj/n;
	
	/* estimate add_var_tilda as the average of the outer product of the adjuted trait
	 values: [ (D-W*beta_tilda)_i * (D-W*beta_tilda)_j ] / phi[i][j]
	 */
	for(i=1; i<= n_fam; i++) {
		n_pheno = family[i].n_pheno;
		
		pheno_typed = family[i].pheno_typed;
		cov = family[i].cov;
		trait = family[i].trait;
		phi = family[i].phi;
		
		for(j=2; j<=n_pheno;j++) {
			for(k=1;k<j;k++) {
				if(phi[j][k]>1e-10) {
					unr++;
					outer_dot = outer_dot + (trait[j]-cov[j][1]*beta_tilda)*(trait[k]-cov[k][1]*beta_tilda)/phi[j][k];
				}
			}
		}
	}
	
	add_var_tilda = outer_dot/unr;
	
	if (add_var_tilda <= 0.0) {
		add_var_tilda = 1e-6;
	}
	
	mle_initial->ax = add_var_tilda/err_var_tilda;
}

/* inital_brak assumes that initial_ax has already been called, this function calls the 
 brak function (from GLS min/bracketing.c) to bracket the minimum (so that brent works)
 */
void initial_brak(struct MLE_PARAM *mle_initial, double (*func)(double)) {
	
	double ax, bx, cx, f_ax, f_bx, f_cx;
	ax = mle_initial->ax;
	
	if (ax < 1) {
		bx = 10+1000.0*ax;
		ax = 0.0;;
	} else {
		bx = 0.0;
		ax = 1000.0*ax;
	}
	
	brak(&ax, &bx, &cx, &f_ax, &f_bx, &f_cx, func);
	
	mle_initial->ax = ax;
	mle_initial->bx = bx;
	mle_initial->cx = cx;	
}


/** brak function adapted from GSL min/bracketing.c */

void brak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double))
{

  /* The three following variables must be declared volatile to avoid storage
     in extended precision registers available on some architecture. The code
     relies on the ability to compare double values. As the values will be
     store in regular memory, the extended precision will then be lost and
     values that are different in extended precision might have equal
     representation in double precision. This behavior might break the
     algorithm. 
   */
  double x_left = (*ax < *bx) ? *ax : *bx;
  double x_right= (*ax < *bx) ? *bx : *ax;
  double x_center;

  volatile double f_left = (*func)(x_left);
  volatile double f_right = (*func)(x_right);
  volatile double f_center;

  const double golden = 0.3819660;      /* golden = (3 - sqrt(5))/2 */
  size_t nb_eval = 0;
  
  
  if (f_right >= f_left) {
    x_center = (x_right - x_left) * golden + x_left;
    nb_eval++;
    f_center = (*func)(x_center);
  } else {
    x_center = x_right;
    f_center = f_right;
    x_right = (x_center - x_left) / golden + x_left;
    nb_eval++;
    f_right = (*func)(x_right);
  }
  
  do {
    if (f_center < f_left ) {
      if (f_center < f_right) {
	*ax = x_left;
	*cx = x_right;
	*bx = x_center;
	*fa = f_left;
	*fc = f_right;
	*fb = f_center;
	break;
      } else if (f_center > f_right) {
	x_left = x_center;
	f_left = f_center;
	x_center = x_right;
	f_center = f_right;
	x_right = (x_center - x_left) / golden + x_left;
	nb_eval++;
	f_right = (*func)(x_right);
      } else { /* f_center == f_right */
	x_right = x_center;
	f_right = f_center;
	x_center = (x_right - x_left) * golden + x_left;
	nb_eval++;
	f_center = (*func)(x_center);
      }
    } else { /* f_center >= f_left */
      x_right = x_center;
      f_right = f_center;
      x_center = (x_right - x_left) * golden + x_left;
      nb_eval++;
      f_center = (*func)(x_center);
    }
  } while (nb_eval < EVAL_MAX
         && (x_right - x_left) > GSL_SQRT_DBL_EPSILON * ( (x_right + x_left) * 0.5 ) + GSL_SQRT_DBL_EPSILON);

  
  *ax = x_left;
  *cx = x_right;
  *bx = x_center;
  *fa = f_left;
  *fc = f_right;
  *fb = f_center;
}

