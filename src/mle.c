#include "mle.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nrutil.h"
#include "cholesky.h"

extern struct OPT_STRUCT f_opt_struct;

double f_neg_log_like_mle(double alpha) {
	return neg_log_like_mle(alpha, f_opt_struct.svd, f_opt_struct.data_struct);
}

/* neg_log_like_mle returns the negative log-likelihood when the formula for the MLE
 error variance (sigam_e^2) has been plugged in, this function is minimized as a
 function of alpha = sigma_a^2/sigma_e^2

 family, svd are the whole structures, i.e. pointer on the structure for all families
 (not only a single family)

 -log_like = n*log(err_var) + sum_i(log(alpha*lambda_i+1))
 */
double neg_log_like_mle(double alpha, struct FAMILY_SVD *svd, struct DATA_STRUCT data_struct) {

	int n=0, n_f=0, i,j, n_fam = data_struct.n_fam;
	double err_var_hat, value = 0.0;

	if (alpha >= 0) {
		err_var_hat = calc_err_var_hat(alpha, svd, data_struct);

		value = 0.0;

		for(i=1; i<= n_fam; i++) {
			n_f = svd[i].n_ind;
			n = n + n_f;
			for(j=1; j<=n_f; j++) {
				value = value + log(alpha*svd[i].lambda[j]+1.0);
			}
		}

		value = n*log(err_var_hat) + value;

	} else {
		err_var_hat = calc_err_var_hat(0.0, svd, data_struct);

		value = 0.0;

		for(i=1; i<= n_fam; i++) {
			n_f = svd[i].n_ind;
			n = n + n_f;
		}

		value = n*log(err_var_hat);

		value = value - 100.0*alpha;

	}

	return value;
}

/* calc_err_var_hat returns the estimate of the error variance as a function of alpha
 the formula when n_cov = 1 (W is a vector):

 (D^tK^{-1}D - D^tK^{-1}W(W^tK^{-1}W)^{-1}W^tK^{-1}D)/n =
 1/n(sum_i(trait_u_2)_i - sum_i(trait_cov_u)_i * (sum_i(cov_u_2)_i)^{-1} * sum_i(trait_cov_u)_i)

 where (U Lambda U^t = Phi):
 trait_u_2_i = 1/(alpha*lambda_i+1) * ([D^tU]_i)^{2}
 cov_u_2_i = 1/(alpha*lambda_i+1) * ([W^tU]_i)^{2}
 trait_cov_u_i = 1/(alpha*lambda_i+1) * ([D^tU]_i)*([W^tU]_i)

 [D^tU]_i = svd[fam].trait_u[i] already initialized in function fill_struct_svd
 [W^tU]_i = svd[fam].cov_u[i][1] already initialized in function fill_struct_svd

 */
double calc_err_var_hat(double alpha, struct FAMILY_SVD *svd, struct DATA_STRUCT data_struct) {

	int n = 0, n_f=0, i, j, i_cov, j_cov, n_fam, n_cov;
	n_fam = data_struct.n_fam;
	n_cov = data_struct.n_cov;


	double trait_u_2=0.0, *trait_cov_u, **cov_u_2;
	trait_cov_u = dvector(1,n_cov);
	cov_u_2 = dmatrix(1,n_cov,1,n_cov);
	for(i_cov=1; i_cov<=n_cov;i_cov++) {
		trait_cov_u[i_cov]=0.0;
		for(j_cov=1; j_cov<=n_cov;j_cov++) {
			cov_u_2[i_cov][j_cov] = 0.0;
		}
	}
	for(i=1; i<= n_fam; i++) {
		n_f = svd[i].n_ind;
		n = n + n_f;
		for(j=1; j<=n_f; j++) {
			trait_u_2 = trait_u_2 + (svd[i].trait_u[j]*svd[i].trait_u[j])/(alpha*svd[i].lambda[j]+1.0);
			for(i_cov=1; i_cov<=n_cov;i_cov++) {
				trait_cov_u[i_cov] = trait_cov_u[i_cov] + svd[i].trait_u[j]*svd[i].cov_u[j][i_cov]/(alpha*svd[i].lambda[j]+1.0);
				for(j_cov=1; j_cov<=n_cov;j_cov++) {
					cov_u_2[i_cov][j_cov] = cov_u_2[i_cov][j_cov] + svd[i].cov_u[j][i_cov]*svd[i].cov_u[j][j_cov]/(alpha*svd[i].lambda[j]+1.0);
				}
			}
		}
	}
	double **chol, **aug, **cholaug;
	chol = dmatrix(1, n_cov, 1, n_cov);
	aug = dmatrix(1, n_cov,1,1);
	cholaug = dmatrix(1,n_cov, 1,1);
	for(i_cov=1; i_cov<=n_cov; i_cov++) {
		aug[i_cov][1] = trait_cov_u[i_cov];
	}

    int posdef=1;
    posdef = cholesky(cov_u_2, n_cov, aug, 1, chol, cholaug, 1);
    if(posdef == 0) {
        printf("\nERROR: Exiting program!\n"
               "The matrix W^T Omega^-1 W is not positive semi-definite\n"
               "(W is the matrix of covariates)\n"
               "Likely this is due to a covarite having the same value for everyone\n"
               "used in the estimation of variance component and fixed effects\n");
        exit(1);
    }


	double mult = 0.0;
	for(i_cov=1; i_cov<=n_cov; i_cov++) {
		mult = mult + cholaug[i_cov][1]*cholaug[i_cov][1];
	}

    free_dvector(trait_cov_u,1,n_cov);
	free_dmatrix(cov_u_2,1,n_cov,1,n_cov);
	free_dmatrix(chol, 1, n_cov, 1,n_cov);
	free_dmatrix(aug, 1, n_cov, 1, 1);
	free_dmatrix(cholaug, 1, n_cov, 1, 1);

	double value;
	value = (trait_u_2-mult)/n;
	return value;
}

/*
 Fill in the data into the structure:

 struct FAMILY_SVD {
 double *trait_u;
 double **cov_u;
 double *lambda;
 };

 trait_u = D^t * U (vector of length n_ind for each family)
 cov_u = U * W^t (matrix, first dimenion for 1...n_ind, second dimension 1...n_cov)
 lambda = eigenvalues of phi_rr (length n_ind)

 n_ind could be n_pheno => subset=0
 n_ind could be n_geno_pheno (we use this for VC estimates for GTAM_est) => subset=1
 n_ind could be n_pheno_genorel => subset=2
 n_ind could be n_pheno_phenorel_w_genorel => subset=3

 */
void fill_struct_svd(struct FAMILY family, struct DATA_STRUCT data_struct, struct FAMILY_SVD *svd, int subset) {

	int i,j, i_cov;

	if (subset == 0) { /* subset 0 is those phenotyped independed of genotype status */

		/* pull out variables from the struct to use internally in the function */
		int n_pheno = family.n_pheno;
		int *pheno_typed = family.pheno_typed;
		double **cov = family.cov;
		double *trait = family.trait;
		double **phi = family.phi;
		int n_cov = data_struct.n_cov;

        if (n_pheno != 0) {

            /* allocate memory for the vectors/matrices of svd */
            svd->trait_u = dvector(1,n_pheno);
            svd->cov_u = dmatrix(1,n_pheno,1,n_cov);
            svd->lambda = dvector(1,n_pheno);

            /* set up Phi_rr (subset of phi for those pheno-/cov-typed) */
            double **phi_rr;
            phi_rr = dmatrix(1,n_pheno,1,n_pheno);
            for(i=1; i<=n_pheno; i++) {
                for(j=1; j<=n_pheno; j++) {
                    phi_rr[i][j] = phi[pheno_typed[i]][pheno_typed[j]];
                }
            }

            double **u;
            u = dmatrix(1,n_pheno,1,n_pheno);

            svdcomp(phi_rr, n_pheno, n_pheno, svd->lambda, u);


            for(i=1; i<=n_pheno; i++) {
                svd->trait_u[i] = 0.0;
                for(i_cov=1; i_cov<=n_cov; i_cov++) {
                    svd->cov_u[i][i_cov] = 0.0;
                }
                for(j=1; j<=n_pheno; j++) {
                    svd->trait_u[i] = svd->trait_u[i] + trait[j]*u[j][i];
                    for(i_cov=1; i_cov<=n_cov; i_cov++) {
                        svd->cov_u[i][i_cov] = svd->cov_u[i][i_cov] + cov[j][i_cov]*u[j][i];
                    }
                }
            }

            svd->n_ind = n_pheno;

            free_dmatrix(u,1,n_pheno,1,n_pheno);
            free_dmatrix(phi_rr,1,n_pheno,1,n_pheno);
        }  else {
			svd->n_ind = 0;
			svd->trait_u = 0;
			svd->cov_u = 0;
			svd->lambda = 0;
		}
	} /* end subset = 0 */

    if (subset == 1) {/* subset 1 is the subset of phenotyped individuals, who are genotyped */

		/* pull out variables from the struct to use internally in the function */
		int n_geno_pheno = family.n_geno_pheno;
		int **geno_pheno_typed = family.geno_pheno_typed;
		double **cov = family.cov;
		double *trait = family.trait;
		double **phi = family.phi;
		int n_cov = data_struct.n_cov;

		if (n_geno_pheno != 0) {
			/* allocate memory for the vectors/matrices of svd */
			svd->trait_u = dvector(1,n_geno_pheno);
			svd->cov_u = dmatrix(1,n_geno_pheno,1,n_cov);
			svd->lambda = dvector(1,n_geno_pheno);

			/* set up Phi_rr (subset of phi for those pheno-/cov-typed) */
			double **phi_rr;
			phi_rr = dmatrix(1,n_geno_pheno,1,n_geno_pheno);
			for(i=1; i<=n_geno_pheno; i++) {
				for(j=1; j<=n_geno_pheno; j++) {
					phi_rr[i][j] = phi[geno_pheno_typed[i][1]][geno_pheno_typed[j][1]];
				}
			}

			double **u;
			u = dmatrix(1,n_geno_pheno,1,n_geno_pheno);

			svdcomp(phi_rr, n_geno_pheno, n_geno_pheno, svd->lambda, u);

			for(i=1; i<=n_geno_pheno; i++) {
				svd->trait_u[i] = 0.0;
				for(i_cov=1; i_cov<=n_cov; i_cov++) {
					svd->cov_u[i][i_cov] = 0.0;
				}
				for(j=1; j<=n_geno_pheno; j++) {
					svd->trait_u[i] = svd->trait_u[i] + trait[geno_pheno_typed[j][2]]*u[j][i];
					for(i_cov=1; i_cov<=n_cov; i_cov++) {
						svd->cov_u[i][i_cov] = svd->cov_u[i][i_cov] + cov[geno_pheno_typed[j][2]][i_cov]*u[j][i];
					}
				}
			}

			svd->n_ind = n_geno_pheno;

			free_dmatrix(u,1,n_geno_pheno,1,n_geno_pheno);
			free_dmatrix(phi_rr,1,n_geno_pheno,1,n_geno_pheno);
		} else {
			svd->n_ind = 0;
			svd->trait_u = 0;
			svd->cov_u = 0;
			svd->lambda = 0;
		}
	} /* end subset = 1*/

	if (subset == 2) { /* subset 2 is the subset of phenotyped individuals, who are genotyped or have genotyped relative */

		/* pull out variables from the struct to use internally in the function */
		int n_pheno = family.n_pheno;
		int *pheno_typed = family.pheno_typed;
        int *geno_typed = family.geno_typed;
		double **cov = family.cov;
		double *trait = family.trait;
		double **phi = family.phi;
		int n_cov = data_struct.n_cov;
        int *pheno_subset = family.pheno_subset;

 		/* count individuals to include */
        int n_pheno_genorel = 0;
      	for(i=1; i<=n_pheno; i++) {
            if(pheno_subset[i] == 1 || pheno_subset[i] == 2) {
                n_pheno_genorel++;
            }
        }

        if (n_pheno_genorel != 0) {
            /* allocate memory for the vectors/matrices of svd */
            svd->trait_u = dvector(1,n_pheno_genorel);
            svd->cov_u = dmatrix(1,n_pheno_genorel,1,n_cov);
            svd->lambda = dvector(1,n_pheno_genorel);

            /* set up Phi_rr (subset of phi for those pheno-/cov-typed and genotyped or with at least one genotyped relative) */
            int i0=0, j0=0;
            double **phi_rr;
            phi_rr = dmatrix(1,n_pheno_genorel,1,n_pheno_genorel);
            for(i=1; i<=n_pheno; i++) {
                j0=0;
                if(pheno_subset[i] == 1 || pheno_subset[i] == 2) {
                    i0++;
                    for(j=1; j<=n_pheno; j++) {
                        if(pheno_subset[j] == 1 || pheno_subset[j] == 2) {
                            j0++;
                            phi_rr[i0][j0] = phi[pheno_typed[i]][pheno_typed[j]];
                        }
                    }
                }
            }

            /* ADDED TO temporarily to TEST NOV 5*/
            /*
            double **u,**v;
            u = dmatrix(1,n_pheno_genorel,1,n_pheno_genorel);
            v= dmatrix(1,n_pheno_genorel,1,n_pheno_genorel);

            svdcomp_u_and_v(phi_rr, n_pheno_genorel, n_pheno_genorel,  svd->lambda, u, v);
            i0=0, j0=0;
            for(i=1; i<=n_pheno; i++) {
                j0=0;
                if(pheno_subset[i] == 1 || pheno_subset[i] == 2) {
                    i0++;
                    svd->trait_u[i0] = 0.0;
                    for(i_cov=1; i_cov<=n_cov; i_cov++) {
                        svd->cov_u[i0][i_cov] = 0.0;
                    }
                    for(j=1; j<=n_pheno; j++) {
                        if(pheno_subset[j] == 1 || pheno_subset[j] == 2) {
                            j0++;
                            svd->trait_u[i0] = svd->trait_u[i0] + trait[j]*u[j0][i0];
                            for(i_cov=1; i_cov<=n_cov; i_cov++) {
                                svd->cov_u[i0][i_cov] = svd->cov_u[i0][i_cov] + cov[j][i_cov]*u[j0][i0];
                            }
                        }
                    }
                }
            }

            double eucleadDist2=0;
            for(i=1; i<=n_pheno_genorel;i++) {
                for(j=1;j<=n_pheno_genorel;j++) {
                    eucleadDist2 = eucleadDist2 + (u[i][j]-v[i][j])*(u[i][j]-v[i][j]);
                }
            }

            printf("The squared Eucleadian distance between the left and right %f\n",eucleadDist2);fflush(stdout);
            */

            double **u;
            u = dmatrix(1,n_pheno_genorel,1,n_pheno_genorel);

            svdcomp(phi_rr, n_pheno_genorel, n_pheno_genorel, svd->lambda, u);

            /* fill svd structure for those pheno-/cov-type and genotyped or with at leaset one genotyped relative */

            i0=0, j0=0;
            for(i=1; i<=n_pheno; i++) {
                j0=0;
                if(pheno_subset[i] == 1 || pheno_subset[i] == 2) {
                    i0++;
                    svd->trait_u[i0] = 0.0;
                    for(i_cov=1; i_cov<=n_cov; i_cov++) {
                        svd->cov_u[i0][i_cov] = 0.0;
                    }
                    for(j=1; j<=n_pheno; j++) {
                        if(pheno_subset[j] == 1 || pheno_subset[j] == 2) {
                            j0++;
                            svd->trait_u[i0] = svd->trait_u[i0] + trait[j]*u[j0][i0];
                            for(i_cov=1; i_cov<=n_cov; i_cov++) {
                                svd->cov_u[i0][i_cov] = svd->cov_u[i0][i_cov] + cov[j][i_cov]*u[j0][i0];
                            }
                        }
                    }
                }
            }

            /* This was a bug and is wrong. Discovered by Sheng Zhong */
            /*	for(i=1; i<=n_pheno_genorel; i++) {
             svd->trait_u[i] = 0.0;
             for(i_cov=1; i_cov<=n_cov; i_cov++) {
             svd->cov_u[i][i_cov] = 0.0;
             }
             for(j=1; j<=n_pheno_genorel; j++) {
             svd->trait_u[i] = svd->trait_u[i] + trait[j]*u[j][i];
             for(i_cov=1; i_cov<=n_cov; i_cov++) {
             svd->cov_u[i][i_cov] = svd->cov_u[i][i_cov] + cov[j][i_cov]*u[j][i];
             }
             }
             }*/

            svd->n_ind = n_pheno_genorel;

            free_dmatrix(u,1,n_pheno_genorel,1,n_pheno_genorel);
            free_dmatrix(phi_rr,1,n_pheno_genorel,1,n_pheno_genorel);

        } else {
			svd->n_ind = 0;
			svd->trait_u = 0;
			svd->cov_u = 0;
			svd->lambda = 0;
		}

	} /* end subset = 2*/

    if (subset == 3) {

        /* pull out variables from the struct to use internally in the function */
		int n_pheno = family.n_pheno;
		int *pheno_typed = family.pheno_typed;
        int *geno_typed = family.geno_typed;
		double **cov = family.cov;
		double *trait = family.trait;
		double **phi = family.phi;
		int n_cov = data_struct.n_cov;
        int *pheno_subset = family.pheno_subset;

 		/* count individuals to include */
        int n_pheno_phenorel_w_genorel = 0;
      	for(i=1; i<=n_pheno; i++) {
            if(pheno_subset[i] == 1 || pheno_subset[i] == 2 || pheno_subset[i] == 3) {
                n_pheno_phenorel_w_genorel++;
            }
        }

        if (n_pheno_phenorel_w_genorel != 0) {
            /* allocate memory for the vectors/matrices of svd */
            svd->trait_u = dvector(1,n_pheno_phenorel_w_genorel);
            svd->cov_u = dmatrix(1,n_pheno_phenorel_w_genorel,1,n_cov);
            svd->lambda = dvector(1,n_pheno_phenorel_w_genorel);

            /* set up Phi_rr (subset of phi for those pheno-/cov-typed and genotyped, with a genotyped relative, or phenotyped relative with genotyped relative) */
            int i0=0, j0=0;
            double **phi_rr;
            phi_rr = dmatrix(1,n_pheno_phenorel_w_genorel,1,n_pheno_phenorel_w_genorel);
            for(i=1; i<=n_pheno; i++) {
                j0=0;
                if(pheno_subset[i] == 1 || pheno_subset[i] == 2 || pheno_subset[i] == 3) {
                    i0++;
                    for(j=1; j<=n_pheno; j++) {
                        if(pheno_subset[j] == 1 || pheno_subset[j] == 2 || pheno_subset[j] == 3) {
                            j0++;
                            phi_rr[i0][j0] = phi[pheno_typed[i]][pheno_typed[j]];
                        }
                    }
                }
            }

            double **u;
            u = dmatrix(1,n_pheno_phenorel_w_genorel,1,n_pheno_phenorel_w_genorel);

            svdcomp(phi_rr, n_pheno_phenorel_w_genorel, n_pheno_phenorel_w_genorel, svd->lambda, u);

            /* fill svd structure for those pheno-/cov-type and genotyped or with at leaset one genotyped relative */
            i0=0, j0=0;
            for(i=1; i<=n_pheno; i++) {
                j0=0;
                if(pheno_subset[i] == 1 || pheno_subset[i] == 2 || pheno_subset[i] == 3) {
                    i0++;
                    svd->trait_u[i0] = 0.0;
                    for(i_cov=1; i_cov<=n_cov; i_cov++) {
                        svd->cov_u[i0][i_cov] = 0.0;
                    }
                    for(j=1; j<=n_pheno; j++) {
                        if(pheno_subset[j] == 1 || pheno_subset[j] == 2 || pheno_subset[j] == 3) {
                            j0++;
                            svd->trait_u[i0] = svd->trait_u[i0] + trait[j]*u[j0][i0];
                            for(i_cov=1; i_cov<=n_cov; i_cov++) {
                                svd->cov_u[i0][i_cov] = svd->cov_u[i0][i_cov] + cov[j][i_cov]*u[j0][i0];
                            }
                        }
                    }
                }
            }

            svd->n_ind = n_pheno_phenorel_w_genorel;

            free_dmatrix(u,1,n_pheno_phenorel_w_genorel,1,n_pheno_phenorel_w_genorel);
            free_dmatrix(phi_rr,1,n_pheno_phenorel_w_genorel,1,n_pheno_phenorel_w_genorel);

        } else {
			svd->n_ind = 0;
			svd->trait_u = 0;
			svd->cov_u = 0;
			svd->lambda = 0;
		}

    } /* end subset = 3 */

} /* end function fill_struct_svd */

/* estimate_alpha_hat, assumes that the inital values bracket the mininum (i.e. mnbrak,
 initial_brak, initial_ax called already)

 alpha_hat is returned
 */
double estimate_alpha_hat (struct MLE_PARAM param) {

	double ax = param.ax;
	double bx = param.bx;
	double cx = param.cx;
	double tol = param.tol;

	double neg_log_like_mle_value;
	double alpha_hat;

	neg_log_like_mle_value = brent(ax, bx, cx, &f_neg_log_like_mle, tol, &alpha_hat);

	return alpha_hat;
}

/* estimate_mle (alpha, beta's, sigma_a^2, sigma_e^2) assumes alpha_hat has been found
 and structure mle_res updated accordingly
 */
void estimate_mle (struct MLE_RESULTS *mle_res, struct FAMILY_SVD *svd, struct DATA_STRUCT data_struct) {

	int n = 0, n_f=0, i, j, n_fam, n_cov;
	n_fam = data_struct.n_fam;
	n_cov = data_struct.n_cov;
	double alpha = mle_res->alpha_hat;
	double denominator_alpha;
    double add, err, var_add, var_err, cov_add_err;


    int i_cov, j_cov;
	double trait_u_2=0, i_22=0, i_23=0, i_33=0;
	double *trait_cov_u = dvector(1,n_cov);
	double **cov_u_2 = dmatrix(1,n_cov,1,n_cov);
	for(i_cov=1; i_cov<=n_cov; i_cov++) {
		trait_cov_u[i_cov] = 0;
		for(j_cov=1;j_cov<=n_cov;j_cov++) {
			cov_u_2[i_cov][j_cov] = 0;
		}
	}

	for(i=1; i<= n_fam; i++) {
		n_f = svd[i].n_ind;
		n = n + n_f;
		for(j=1; j<=n_f; j++) {
			trait_u_2 = trait_u_2 + (svd[i].trait_u[j]*svd[i].trait_u[j])/(alpha*svd[i].lambda[j]+1.0);
			denominator_alpha = (alpha*svd[i].lambda[j]+1.0)*(alpha*svd[i].lambda[j]+1.0);
			i_22 = i_22 + (svd[i].lambda[j]*svd[i].lambda[j])/(2.0*denominator_alpha);
			i_23 = i_23 + (svd[i].lambda[j])/(2.0*denominator_alpha);
			i_33 = i_33 + 1.0/(2.0*denominator_alpha);
		}
	}

	for	(i_cov=1; i_cov<=n_cov; i_cov++) {
		for(i=1; i<= n_fam; i++) {
			n_f = svd[i].n_ind;
			for(j=1; j<=n_f; j++) {
				trait_cov_u[i_cov] = trait_cov_u[i_cov] + (svd[i].trait_u[j]*svd[i].cov_u[j][i_cov])/(alpha*svd[i].lambda[j]+1.0);
			}
		}
		for(j_cov=1; j_cov<=n_cov; j_cov++) {
			for(i=1; i<= n_fam; i++) {
				n_f = svd[i].n_ind;
				for(j=1; j<=n_f; j++) {
					cov_u_2[i_cov][j_cov] = cov_u_2[i_cov][j_cov] + (svd[i].cov_u[j][i_cov]*svd[i].cov_u[j][j_cov])/(alpha*svd[i].lambda[j]+1.0);
				}
			}
		}
	}


	double **chol=dmatrix(1,n_cov,1,n_cov);
	double **aug=dmatrix(1,n_cov,1,n_cov+1);
	double **cholaug=dmatrix(1,n_cov,1,n_cov+1);

	for(i_cov=1;i_cov<=n_cov;i_cov++) {
		for(j_cov=1;j_cov<=n_cov;j_cov++) {
			chol[i_cov][j_cov]=0.0;
			aug[i_cov][j_cov]=0.0;
			cholaug[i_cov][j_cov]=0.0;
		}
		aug[i_cov][i_cov] = 1.0;
		aug[i_cov][n_cov+1] = trait_cov_u[i_cov];
	}

    int posdef=1;
    posdef = cholesky(cov_u_2,n_cov, aug, n_cov+1, chol, cholaug,1);
    if(posdef == 0) {
        printf("\nERROR: Exiting program!\n"
               "The matrix W^T Omega^-1 W is not positive semi-definite\n"
               "(W is the matrix of covariates)\n"
               "Likely this is due to a covarite having the same value for everyone\n"
               "used in the estimation of variance component and fixed effects\n");
        exit(1);
    }

	double err_tmp = 0.0;
	for(i_cov=1; i_cov<=n_cov;i_cov++) {
		err_tmp = err_tmp+cholaug[i_cov][n_cov+1]*cholaug[i_cov][n_cov+1];
	}

	mle_res->err_var_hat = (trait_u_2-err_tmp)/n;
	mle_res->add_var_hat = alpha*mle_res->err_var_hat;

	for(i_cov=1;i_cov<=n_cov;i_cov++) {
		mle_res->beta_hat[i_cov] = 0.0;
		for(j_cov=1;j_cov<=n_cov;j_cov++) {
			mle_res->beta_hat[i_cov] = mle_res->beta_hat[i_cov] + cholaug[j_cov][i_cov]*cholaug[j_cov][n_cov+1];
		}
	}

	int cov;
	double **beta_var_mat = dmatrix(1,n_cov,1,n_cov);
	for(i_cov=1;i_cov<=n_cov;i_cov++) {
		for(j_cov=1;j_cov<=n_cov;j_cov++) {
			beta_var_mat[i_cov][j_cov] = 0.0;
			for(cov=1;cov<=n_cov;cov++) {
				beta_var_mat[i_cov][j_cov] = beta_var_mat[i_cov][j_cov] + cholaug[cov][i_cov]*cholaug[cov][j_cov];
			}
		}
	}

    add = mle_res->add_var_hat;
    err = mle_res->err_var_hat;
    mle_res->herit = add/(add+err);

	i_22 = i_22/(mle_res->err_var_hat*mle_res->err_var_hat);
	i_23 = i_23/(mle_res->err_var_hat*mle_res->err_var_hat);
	i_33 = i_33/(mle_res->err_var_hat*mle_res->err_var_hat);

	mle_res->var_err_var_hat = i_22/(i_22*i_33-i_23*i_23);
	mle_res->var_add_var_hat = i_33/(i_22*i_33-i_23*i_23);

	for(i_cov=1;i_cov<=n_cov;i_cov++) {
		mle_res->var_beta_hat[i_cov] = mle_res->err_var_hat*beta_var_mat[i_cov][i_cov];
	}

    var_add = mle_res->var_add_var_hat;
    var_err = mle_res->var_err_var_hat;

    cov_add_err = -i_23/(i_22*i_33-i_23*i_23);

    mle_res->var_herit = (err*err*var_add+add*add*var_err-2*add*err*cov_add_err)/((add+err)*(add+err)*(add+err)*(add+err));


    free_dmatrix(beta_var_mat,1,n_cov,1,n_cov);
    free_dmatrix(cov_u_2,1,n_cov,1,n_cov);
    free_dvector(trait_cov_u,1,n_cov);
    free_dmatrix(chol,1,n_cov,1,n_cov);
    free_dmatrix(aug,1,n_cov,1,n_cov+1);
    free_dmatrix(cholaug,1,n_cov,1,n_cov+1);

}

void estimate_mle_ols (struct MLE_RESULTS *mle_res, struct FAMILY *family, struct DATA_STRUCT data_struct, int subset, int alpha_neg) {

	int n = 0, n_f=0, i, j, i_cov, j_cov, n_fam, n_cov;
	n_fam = data_struct.n_fam;
	n_cov = data_struct.n_cov;

	/* estimate the mle using everyone who is phenotyped */
	if(subset == 0) {


        double err_var_numerator=0.0;
        double **ww = dmatrix(1,n_cov,1,n_cov);
        double *wd = dvector(1,n_cov);

        for (i_cov=1;i_cov<=n_cov;i_cov++) {
            wd[i_cov] = 0.0;
            for (j_cov=1; j_cov<=n_cov;j_cov++) {
                ww[i_cov][j_cov] = 0.0;
            }
        }

        for(i=1; i<= n_fam; i++) {
            n_f = family[i].n_pheno;
            for(j=1; j<=n_f; j++) {
                for(i_cov=1;i_cov<=n_cov;i_cov++) {
                    wd[i_cov] = wd[i_cov] + family[i].cov[j][i_cov]*family[i].trait[j];
                    for(j_cov=1;j_cov<=n_cov;j_cov++) {
                        ww[i_cov][j_cov] = ww[i_cov][j_cov] + family[i].cov[j][i_cov]*family[i].cov[j][j_cov];
                    }
                }
            }
        }

        double **chol=dmatrix(1,n_cov,1,n_cov);
        double **aug=dmatrix(1,n_cov,1,n_cov+1);
        double **cholaug=dmatrix(1,n_cov,1,n_cov+1);

        for(i_cov=1;i_cov<=n_cov;i_cov++) {
            for(j_cov=1;j_cov<=n_cov;j_cov++) {
                chol[i_cov][j_cov]=0.0;
                aug[i_cov][j_cov]=0.0;
                cholaug[i_cov][j_cov]=0.0;
            }
            aug[i_cov][i_cov] = 1.0;
            aug[i_cov][n_cov+1] = wd[i_cov];
        }

        int posdef=1;
        posdef = cholesky(ww,n_cov, aug, n_cov+1, chol, cholaug,1);
        if(posdef == 0) {
            printf("\nERROR: Exiting program!\n"
                    "The matrix W^T W is not positive semi-definite\n"
                    "(W is the matrix of covariates)\n"
                    "Likely this is due to a covarite having the same value for everyone\n"
                    "used in the estimation of variance component and fixed effects\n");
            exit(1);
        }


        for(i_cov=1;i_cov<=n_cov;i_cov++) {
            mle_res->beta_ols[i_cov] = 0.0;
            for(j_cov=1;j_cov<=n_cov;j_cov++) {
                mle_res->beta_ols[i_cov] = mle_res->beta_ols[i_cov] + cholaug[j_cov][i_cov]*cholaug[j_cov][n_cov+1];
            }
        }

        if (alpha_neg==1) {

            mle_res->beta_hat = mle_res->beta_ols;

            int cov_adj;
            for(i=1; i<=n_fam; i++) {
                n_f = family[i].n_pheno;
                n = n + n_f;
                for(j=1; j<=n_f; j++) {
                    cov_adj = 0.0;
                    for (i_cov=1; i_cov<= n_cov; i_cov++) {
                        cov_adj = cov_adj + family[i].cov[j][i_cov]*mle_res->beta_hat[i_cov];
                    }
                    err_var_numerator = err_var_numerator + (family[i].trait[j]-cov_adj)*(family[i].trait[j]-cov_adj);
                }
            }

            mle_res->add_var_hat = 0.0;
            mle_res->err_var_hat = err_var_numerator/n;
            mle_res->herit = 0.0;

            for(i_cov=1;i_cov<=n_cov;i_cov++) {
                mle_res->var_beta_hat[i_cov] = 0.0;
            }
            mle_res->var_add_var_hat = 0.0;
            mle_res->var_err_var_hat = 0.0;
            mle_res->var_herit = 0.0;
        }

        free_dmatrix(ww,1,n_cov,1,n_cov);
        free_dvector(wd,1,n_cov);
        free_dmatrix(chol,1,n_cov,1,n_cov);
        free_dmatrix(aug,1,n_cov,1,n_cov+1);
        free_dmatrix(cholaug,1,n_cov,1,n_cov+1);

	} /* end if subset==0 */

    /* estimate using those phenotyped and genotyped */
	if(subset == 1) {


		double err_var_numerator=0.0;
        double **ww = dmatrix(1,n_cov,1,n_cov);
		double *wd = dvector(1,n_cov);

        for (i_cov=1;i_cov<=n_cov;i_cov++) {
			wd[i_cov] = 0.0;
			for (j_cov=1; j_cov<=n_cov;j_cov++) {
				ww[i_cov][j_cov] = 0.0;
			}
        }

		for(i=1; i<= n_fam; i++) {
			n_f = family[i].n_geno_pheno;
			for(j=1; j<=n_f; j++) {
				for(i_cov=1;i_cov<=n_cov;i_cov++) {
					wd[i_cov] = wd[i_cov] + family[i].cov[family[i].geno_pheno_typed[j][2]][i_cov]*family[i].trait[family[i].geno_pheno_typed[j][2]];
					for(j_cov=1;j_cov<=n_cov;j_cov++) {
						ww[i_cov][j_cov] = ww[i_cov][j_cov] + family[i].cov[family[i].geno_pheno_typed[j][2]][i_cov]*family[i].cov[family[i].geno_pheno_typed[j][2]][j_cov];
					}
				}
			}
		}

		double **chol=dmatrix(1,n_cov,1,n_cov);
		double **aug=dmatrix(1,n_cov,1,n_cov+1);
		double **cholaug=dmatrix(1,n_cov,1,n_cov+1);

		for(i_cov=1;i_cov<=n_cov;i_cov++) {
			for(j_cov=1;j_cov<=n_cov;j_cov++) {
				chol[i_cov][j_cov]=0.0;
				aug[i_cov][j_cov]=0.0;
				cholaug[i_cov][j_cov]=0.0;
			}
			aug[i_cov][i_cov] = 1.0;
			aug[i_cov][n_cov+1] = wd[i_cov];
		}


        int posdef=1;
        posdef = cholesky(ww,n_cov, aug, n_cov+1, chol, cholaug,1);
        if(posdef == 0) {
            printf("\nERROR: Exiting program!\n"
                   "The matrix W^T W is not positive semi-definite\n"
                   "(W is the matrix of covariates)\n"
                   "Likely this is due to a covarite having the same value for everyone\n"
                   "used in the estimation of variance component and fixed effects\n");
            exit(1);
        }


        for(i_cov=1;i_cov<=n_cov;i_cov++) {
			mle_res->beta_ols[i_cov] = 0.0;
			for(j_cov=1;j_cov<=n_cov;j_cov++) {
				mle_res->beta_ols[i_cov] = mle_res->beta_ols[i_cov] + cholaug[j_cov][i_cov]*cholaug[j_cov][n_cov+1];
			}
		}

		if (alpha_neg==1) {

            mle_res->beta_hat = mle_res->beta_ols;

            n = 0;

			int cov_adj;
			for(i=1; i<=n_fam; i++) {
				n_f = family[i].n_pheno;
				n = n + n_f;
				for(j=1; j<=n_f; j++) {
					cov_adj = 0.0;
					for (i_cov=1; i_cov<= n_cov; i_cov++) {
						cov_adj = cov_adj + family[i].cov[family[i].geno_pheno_typed[j][2]][i_cov]*mle_res->beta_hat[i_cov];
					}
					err_var_numerator = err_var_numerator + (family[i].trait[family[i].geno_pheno_typed[j][2]]-cov_adj)*(family[i].trait[family[i].geno_pheno_typed[j][2]]-cov_adj);
                }
            }

			mle_res->add_var_hat = 0.0;
			mle_res->err_var_hat = err_var_numerator/n;
            mle_res->herit = 0.0;


			for(i_cov=1;i_cov<=n_cov;i_cov++) {
				mle_res->var_beta_hat[i_cov] = 0.0;
			}
			mle_res->var_add_var_hat = 0.0;
			mle_res->var_err_var_hat = 0.0;
            mle_res->var_herit = 0.0;

        }

		free_dmatrix(ww,1,n_cov,1,n_cov);
		free_dvector(wd,1,n_cov);
		free_dmatrix(chol,1,n_cov,1,n_cov);
		free_dmatrix(aug,1,n_cov,1,n_cov+1);
		free_dmatrix(cholaug,1,n_cov,1,n_cov+1);

	} /* end if subset==1 */

	/* estimate the mle using everyone who is phenotyped and either genotyped or with at least one genotyped relative */
	if(subset == 2) {

        double err_var_numerator=0.0;
        double **ww = dmatrix(1,n_cov,1,n_cov);
        double *wd = dvector(1,n_cov);

        /* initialize */
        for (i_cov=1;i_cov<=n_cov;i_cov++) {
            wd[i_cov] = 0.0;
            for (j_cov=1; j_cov<=n_cov;j_cov++) {
                ww[i_cov][j_cov] = 0.0;
            }
        }

        for(i=1; i<= n_fam; i++) {
            n_f = family[i].n_pheno;
            for(j=1; j<=n_f; j++) {
                if(family[i].pheno_subset[j] == 1 || family[i].pheno_subset[j] == 2) {
                    for(i_cov=1;i_cov<=n_cov;i_cov++) {
                        wd[i_cov] = wd[i_cov] + family[i].cov[j][i_cov]*family[i].trait[j];
                        for(j_cov=1;j_cov<=n_cov;j_cov++) {
                            ww[i_cov][j_cov] = ww[i_cov][j_cov] + family[i].cov[j][i_cov]*family[i].cov[j][j_cov];
                        }
                    }
                }
            }
        }


        double **chol=dmatrix(1,n_cov,1,n_cov);
        double **aug=dmatrix(1,n_cov,1,n_cov+1);
        double **cholaug=dmatrix(1,n_cov,1,n_cov+1);

        for(i_cov=1;i_cov<=n_cov;i_cov++) {
            for(j_cov=1;j_cov<=n_cov;j_cov++) {
                chol[i_cov][j_cov]=0.0;
                aug[i_cov][j_cov]=0.0;
                cholaug[i_cov][j_cov]=0.0;
            }
            aug[i_cov][i_cov] = 1.0;
            aug[i_cov][n_cov+1] = wd[i_cov];
        }

        int posdef=1;
        posdef = cholesky(ww,n_cov, aug, n_cov+1, chol, cholaug,1);
        if(posdef == 0) {
            printf("\nERROR: Exiting program!\n"
                   "The matrix W^T W is not positive semi-definite\n"
                   "(W is the matrix of covariates)\n"
                   "Likely this is due to a covarite having the same value for everyone\n"
                   "used in the estimation of variance component and fixed effects\n");
            exit(1);
        }


        for(i_cov=1;i_cov<=n_cov;i_cov++) {
            mle_res->beta_ols[i_cov] = 0.0;
            for(j_cov=1;j_cov<=n_cov;j_cov++) {
                mle_res->beta_ols[i_cov] = mle_res->beta_ols[i_cov] + cholaug[j_cov][i_cov]*cholaug[j_cov][n_cov+1];
            }
        }

        if (alpha_neg==1) {

            mle_res->beta_hat = mle_res->beta_ols;

            n = 0;

            int cov_adj;
            for(i=1; i<=n_fam; i++) {
                n_f = family[i].n_pheno;
                n = n + n_f;
                for(j=1; j<=n_f; j++) {
                    cov_adj = 0.0;
                    if(family[i].pheno_subset[j]==1 || family[i].pheno_subset[j]==2) {
                        for (i_cov=1; i_cov<= n_cov; i_cov++) {
                            cov_adj = cov_adj + family[i].cov[j][i_cov]*mle_res->beta_hat[i_cov];
                        }
                        err_var_numerator = err_var_numerator + (family[i].trait[j]-cov_adj)*(family[i].trait[j]-cov_adj);
                    }
                }
            }

            mle_res->add_var_hat = 0.0;
            mle_res->err_var_hat = err_var_numerator/n;
            mle_res->herit = 0.0;


            for(i_cov=1;i_cov<=n_cov;i_cov++) {
                mle_res->var_beta_hat[i_cov] = 0.0;
            }
            mle_res->var_add_var_hat = 0.0;
            mle_res->var_err_var_hat = 0.0;
            mle_res->var_herit = 0.0;
        }

        free_dmatrix(ww,1,n_cov,1,n_cov);
        free_dvector(wd,1,n_cov);
        free_dmatrix(chol,1,n_cov,1,n_cov);
        free_dmatrix(aug,1,n_cov,1,n_cov+1);
        free_dmatrix(cholaug,1,n_cov,1,n_cov+1);


	} /* end if subset==2 */

    /* estimate the mle using everyone who is phenotyped and either genotyped, with genotyped relative, or phenotyped relative with genotyped relative */
	if(subset == 3) {

        double err_var_numerator=0.0;
        double **ww = dmatrix(1,n_cov,1,n_cov);
        double *wd = dvector(1,n_cov);

        /* initialize */
        for (i_cov=1;i_cov<=n_cov;i_cov++) {
            wd[i_cov] = 0.0;
            for (j_cov=1; j_cov<=n_cov;j_cov++) {
                ww[i_cov][j_cov] = 0.0;
            }
        }

        for(i=1; i<= n_fam; i++) {
            n_f = family[i].n_pheno;
            for(j=1; j<=n_f; j++) {
                if(family[i].pheno_subset[j] == 1 || family[i].pheno_subset[j] == 2 || family[i].pheno_subset[j] == 3) {
                    for(i_cov=1;i_cov<=n_cov;i_cov++) {
                        wd[i_cov] = wd[i_cov] + family[i].cov[j][i_cov]*family[i].trait[j];
                        for(j_cov=1;j_cov<=n_cov;j_cov++) {
                            ww[i_cov][j_cov] = ww[i_cov][j_cov] + family[i].cov[j][i_cov]*family[i].cov[j][j_cov];
                        }
                    }
                }
            }
        }


        double **chol=dmatrix(1,n_cov,1,n_cov);
        double **aug=dmatrix(1,n_cov,1,n_cov+1);
        double **cholaug=dmatrix(1,n_cov,1,n_cov+1);

        for(i_cov=1;i_cov<=n_cov;i_cov++) {
            for(j_cov=1;j_cov<=n_cov;j_cov++) {
                chol[i_cov][j_cov]=0.0;
                aug[i_cov][j_cov]=0.0;
                cholaug[i_cov][j_cov]=0.0;
            }
            aug[i_cov][i_cov] = 1.0;
            aug[i_cov][n_cov+1] = wd[i_cov];
        }

        int posdef=1;
        posdef = cholesky(ww,n_cov, aug, n_cov+1, chol, cholaug,1);
        if(posdef == 0) {
            printf("\nERROR: Exiting program!\n"
                   "The matrix W^T W is not positive semi-definite\n"
                   "(W is the matrix of covariates)\n"
                   "Likely this is due to a covarite having the same value for everyone\n"
                   "used in the estimation of variance component and fixed effects\n");
            exit(1);
        }


        for(i_cov=1;i_cov<=n_cov;i_cov++) {
            mle_res->beta_ols[i_cov] = 0.0;
            for(j_cov=1;j_cov<=n_cov;j_cov++) {
                mle_res->beta_ols[i_cov] = mle_res->beta_ols[i_cov] + cholaug[j_cov][i_cov]*cholaug[j_cov][n_cov+1];
            }
        }

        if (alpha_neg==1) {

            mle_res->beta_hat = mle_res->beta_ols;

            n = 0;

            int cov_adj;
            for(i=1; i<=n_fam; i++) {
                n_f = family[i].n_pheno;
                n = n + n_f;
                for(j=1; j<=n_f; j++) {
                    cov_adj = 0.0;
                    if(family[i].pheno_subset[j]==1 || family[i].pheno_subset[j]==2 || family[i].pheno_subset[j]==3) {
                        for (i_cov=1; i_cov<= n_cov; i_cov++) {
                            cov_adj = cov_adj + family[i].cov[j][i_cov]*mle_res->beta_hat[i_cov];
                        }
                        err_var_numerator = err_var_numerator + (family[i].trait[j]-cov_adj)*(family[i].trait[j]-cov_adj);
                    }
                }
            }

            mle_res->add_var_hat = 0.0;
            mle_res->err_var_hat = err_var_numerator/n;
            mle_res->herit = 0.0;


            for(i_cov=1;i_cov<=n_cov;i_cov++) {
                mle_res->var_beta_hat[i_cov] = 0.0;
            }
            mle_res->var_add_var_hat = 0.0;
            mle_res->var_err_var_hat = 0.0;
            mle_res->var_herit = 0.0;
        }

        free_dmatrix(ww,1,n_cov,1,n_cov);
        free_dvector(wd,1,n_cov);
        free_dmatrix(chol,1,n_cov,1,n_cov);
        free_dmatrix(aug,1,n_cov,1,n_cov+1);
        free_dmatrix(cholaug,1,n_cov,1,n_cov+1);


	} /* end if subset==3 */


}
