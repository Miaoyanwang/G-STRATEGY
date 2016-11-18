#include "calc.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nrutil.h"
#include "svdcomp.h"
#include "cholesky.h"

/* function to get the a_all vector, a_all = Omega^-1(D-W*beta) if person is in a phenotyped subset */
void a_vector(struct FAMILY *family, struct DATA_STRUCT data_struct, struct MLE_RESULTS mle_res, int subset) {

	/* pull out variables from the struct family to use internally in the function */
	int n_total = family->n_total;
	int n_pheno = family->n_pheno;
	int *pheno_typed = family->pheno_typed;
	int *geno_typed = family->geno_typed;
    int *pheno_subset = family->pheno_subset;
	double **cov = family->cov;
	double *trait = family->trait;
	double **phi = family->phi;
	int n_cov = data_struct.n_cov;
	double *beta = mle_res.beta_hat;

	double err_var, add_var;
	err_var = mle_res.err_var_hat;
	add_var = mle_res.add_var_hat;

	int i, j, k;

	if(subset == 0) {
		/* get trait_adj = D-W*beta */
		double tmp_adj;
		double *trait_adj;
		trait_adj = dvector(1,n_pheno);
		for(i=1; i<=n_pheno; i++) {
			tmp_adj = 0.0;
			for(j=1; j<= n_cov; j++) {
				tmp_adj = tmp_adj + cov[i][j]*beta[j];
			}
			trait_adj[i] = trait[i] - tmp_adj;
		}

		/* set up Omega = err_var*I_rr + add_var*Phi_rr */
		double **omega;
		omega = dmatrix(1,n_pheno,1,n_pheno);
		for(i=1; i<=n_pheno; i++) {
			for(j=1; j<=n_pheno; j++) {
				if(i==j) {
					omega[i][j] = err_var + add_var*phi[pheno_typed[i]][pheno_typed[j]];
				} else {
					omega[i][j] = add_var*phi[pheno_typed[i]][pheno_typed[j]];
				}
			}
		}

        /* singular value decomposition Omega = U Gamma V^T, for symmertic matrices U = V.
         Omega^-1 = U Gamma^-1 V^T */
        double *eigen; eigen = dvector(1,n_pheno);
        double **u; u = dmatrix(1,n_pheno,1,n_pheno);
        svdcomp(omega, n_pheno, n_pheno, eigen, u);

        double omega_inv_tmp=0;
		for (i=1; i<=n_pheno; i++) {
			family->a_store[1][i] = 0.0;
            for (j=1; j<=n_pheno; j++) {
                omega_inv_tmp = 0;
                for(k=1; k<=n_pheno; k++) {
                    omega_inv_tmp = omega_inv_tmp + u[i][k]*u[j][k]/eigen[k];
                }
                family->a_store[1][i] = family->a_store[1][i] + omega_inv_tmp*trait_adj[j];
            }
        }


		/* free matrices */
		free_dmatrix(omega,1,n_pheno,1,n_pheno);
		free_dvector(trait_adj,1,n_pheno);
        free_dmatrix(u,1,n_pheno,1,n_pheno);
        free_dvector(eigen,1,n_pheno);


	} /* end subset = 0 */

	if(subset == 1) {

        int i0=0, j0=0, k0=0;

		/* count individuals to include */
        int n_pheno_geno = 0;
      	for(i=1; i<=n_pheno; i++) {
            if(pheno_subset[i] == 1) {
                n_pheno_geno++;
            }
        }

  		/* get trait_adj = D-W*beta (for subset of phi for those pheno-/cov-typed and genotyped)*/
		double tmp_adj;
		double *trait_adj;
		trait_adj = dvector(1,n_pheno_geno);
		for(i=1; i<=n_pheno; i++) {
			tmp_adj = 0.0;
			if(pheno_subset[i] == 1) {
				i0++;
				for(j=1; j<= n_cov; j++) {
					tmp_adj = tmp_adj + cov[i][j]*beta[j];
				}
				trait_adj[i0] = trait[i] - tmp_adj;
			}
		}

		/* set up Phi_rr (for subset of phi for those pheno-/cov-typed and genotyped) */
    	double **phi_rr;
		i0=0, j0=0;
  		phi_rr = dmatrix(1,n_pheno_geno,1,n_pheno_geno);
		for(i=1; i<=n_pheno; i++) {
            j0=0;
            if(pheno_subset[i] == 1) {
                i0++;
                for(j=1; j<=n_pheno; j++) {
					if(pheno_subset[j] == 1) {
						j0++;
						phi_rr[i0][j0] = phi[pheno_typed[i]][pheno_typed[j]];
					}
				}
            }
		}

		/* set up omega = err_var*I_rr + add_var*phi_rr for the subset */
		double **omega;
		omega = dmatrix(1,n_pheno_geno, 1, n_pheno_geno);
		for(i=1; i<= n_pheno_geno; i++) {
			for (j=1; j<= n_pheno_geno; j++) {
				if (i==j) {
					omega[i][j] =  err_var + add_var * phi_rr[i][j];
				} else {
					omega[i][j] = add_var*phi_rr[i][j];
				}
			}
		}


        /* singular value decomposition Omega = U Gamma V^T, for symmertic matrices U = V.
         Omega^-1 = U Gamma^-1 V^T */
        double *eigen; eigen = dvector(1,n_pheno_geno);
        double **u; u = dmatrix(1,n_pheno_geno,1,n_pheno_geno);
        svdcomp(omega, n_pheno_geno, n_pheno_geno, eigen, u);

        i0=0;j0=0;k0=0;
        double omega_inv_tmp=0;
		for (i=1; i<=n_pheno; i++) {
			family->a_store[1][i] = 0.0;
			if (pheno_subset[i] == 1) {
				i0++;
				j0=0;
				for (j=1; j<= n_pheno; j++) {
					if (pheno_subset[j] == 1) {
						j0++;
						k0=0;
                        omega_inv_tmp = 0;
                        for(k=1; k<=n_pheno; k++) {
                            if (pheno_subset[k] == 1) {
                                k0++;
                                omega_inv_tmp = omega_inv_tmp + u[i0][k0]*u[j0][k0]/eigen[k0];
                            }
                        }
                        family->a_store[1][i] = family->a_store[1][i] + omega_inv_tmp*trait_adj[j0];
                    }
				}
			}
		}

		/* free matrices */
		free_dmatrix(omega,1,n_pheno_geno,1,n_pheno_geno);
		free_dmatrix(phi_rr,1,n_pheno_geno,1,n_pheno_geno);
		free_dvector(trait_adj,1,n_pheno_geno);
        free_dmatrix(u,1,n_pheno_geno,1,n_pheno_geno);
        free_dvector(eigen,1,n_pheno_geno);


	} /* end subset = 1 */


	if(subset == 2) {

		int i0=0, j0=0, k0=0;

		/* count individuals to include */
        int n_pheno_genorel = 0;
      	for(i=1; i<=n_pheno; i++) {
            if(pheno_subset[i] == 1 || pheno_subset[i] == 2) {
                n_pheno_genorel++;
            }
        }

  		/* get trait_adj = D-W*beta (for subset of phi for those pheno-/cov-typed and genotyped or with at least one genotyped relative)*/
		double tmp_adj;
		double *trait_adj;
		trait_adj = dvector(1,n_pheno_genorel);
		for(i=1; i<=n_pheno; i++) {
			tmp_adj = 0.0;
			if(pheno_subset[i] == 1 || pheno_subset[i] == 2) {
				i0++;
				for(j=1; j<= n_cov; j++) {
					tmp_adj = tmp_adj + cov[i][j]*beta[j];
				}
				trait_adj[i0] = trait[i] - tmp_adj;
			}
		}

		/* set up Phi_rr (for subset of phi for those pheno-/cov-typed and genotyped or with at least one genotyped relative) */
    	double **phi_rr;
		i0=0, j0=0;
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

		/* set up omega = err_var*I_rr + add_var*phi_rr for the subset */
		double **omega;
		omega = dmatrix(1,n_pheno_genorel, 1, n_pheno_genorel);
		for(i=1; i<= n_pheno_genorel; i++) {
			for (j=1; j<= n_pheno_genorel; j++) {
				if (i==j) {
					omega[i][j] =  err_var + add_var * phi_rr[i][j];
				} else {
					omega[i][j] = add_var*phi_rr[i][j];
				}
			}
		}


        /* singular value decomposition Omega = U Gamma V^T, for symmertic matrices U = V.
           Omega^-1 = U Gamma^-1 V^T */
        double *eigen; eigen = dvector(1,n_pheno_genorel);
        double **u; u = dmatrix(1,n_pheno_genorel,1,n_pheno_genorel);
        svdcomp(omega, n_pheno_genorel, n_pheno_genorel, eigen, u);

        i0=0;j0=0;k0=0;
        double omega_inv_tmp=0;
		for (i=1; i<=n_pheno; i++) {
			family->a_store[1][i] = 0.0;
			if (pheno_subset[i] == 1 || pheno_subset[i] == 2) {
				i0++;
				j0=0;
				for (j=1; j<= n_pheno; j++) {
					if (pheno_subset[j] == 1 || pheno_subset[j] == 2) {
						j0++;
						k0=0;
                        omega_inv_tmp = 0;
                        for(k=1; k<=n_pheno; k++) {
                            if (pheno_subset[k] == 1 || pheno_subset[k] == 2) {
                                k0++;
                                omega_inv_tmp = omega_inv_tmp + u[i0][k0]*u[j0][k0]/eigen[k0];
                            }
                        }
                        family->a_store[1][i] = family->a_store[1][i] + omega_inv_tmp*trait_adj[j0];
                    }
				}
			}
		}

		/* free matrices */
		free_dmatrix(omega,1,n_pheno_genorel,1,n_pheno_genorel);
		free_dmatrix(phi_rr,1,n_pheno_genorel,1,n_pheno_genorel);
		free_dvector(trait_adj,1,n_pheno_genorel);
        free_dmatrix(u,1,n_pheno_genorel,1,n_pheno_genorel);
        free_dvector(eigen,1,n_pheno_genorel);

	} /* end subset = 2 */

	if(subset == 3) {

        int i0=0, j0=0, k0=0;

        /* count individuals to include */
        int n_pheno_phenorel_w_genorel = 0;
      	for(i=1; i<=n_pheno; i++) {
            if(pheno_subset[i] == 1 || pheno_subset[i] == 2 || pheno_subset[i] == 3) {
                n_pheno_phenorel_w_genorel++;
            }
        }

  		/* get trait_adj = D-W*beta (for subset of phi for those pheno-/cov-typed and genotyped, with at least one genotyped relative, or phenotyped relative with genotyped relative) */
		double tmp_adj;
		double *trait_adj;
		trait_adj = dvector(1,n_pheno_phenorel_w_genorel);
		for(i=1; i<=n_pheno; i++) {
			tmp_adj = 0.0;
			if(pheno_subset[i] == 1 || pheno_subset[i] == 2 || pheno_subset[i] == 3) {
				i0++;
				for(j=1; j<= n_cov; j++) {
					tmp_adj = tmp_adj + cov[i][j]*beta[j];
				}
				trait_adj[i0] = trait[i] - tmp_adj;
			}
		}

		/* set up Phi_rr (for subset of phi for those pheno-/cov-typed and genotyped, with at least one genotyped relative, or phenotyped relative with genotyped relative) */
    	double **phi_rr;
		i0=0, j0=0;
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

		/* set up omega = err_var*I_rr + add_var*phi_rr for the subset */
		double **omega;
		omega = dmatrix(1,n_pheno_phenorel_w_genorel, 1, n_pheno_phenorel_w_genorel);
		for(i=1; i<= n_pheno_phenorel_w_genorel; i++) {
			for (j=1; j<= n_pheno_phenorel_w_genorel; j++) {
				if (i==j) {
					omega[i][j] =  err_var + add_var * phi_rr[i][j];
				} else {
					omega[i][j] = add_var*phi_rr[i][j];
				}
			}
		}

        /* singular value decomposition Omega = U Gamma V^T, for symmertic matrices U = V.
         Omega^-1 = U Gamma^-1 V^T */
        double *eigen; eigen = dvector(1,n_pheno_phenorel_w_genorel);
        double **u; u = dmatrix(1,n_pheno_phenorel_w_genorel,1,n_pheno_phenorel_w_genorel);
        svdcomp(omega, n_pheno_phenorel_w_genorel, n_pheno_phenorel_w_genorel, eigen, u);

        i0=0;j0=0;k0=0;
        double omega_inv_tmp=0;
		for (i=1; i<=n_pheno; i++) {
			family->a_store[1][i] = 0.0;
			if (pheno_subset[i] == 1 || pheno_subset[i] == 2 || pheno_subset[i] == 3) {
				i0++;
				j0=0;
				for (j=1; j<= n_pheno; j++) {
					if (pheno_subset[j] == 1 || pheno_subset[j] == 2 || pheno_subset[j] == 3) {
						j0++;
						k0=0;
                        omega_inv_tmp = 0;
                        for(k=1; k<=n_pheno; k++) {
                            if (pheno_subset[k] == 1 || pheno_subset[k] == 2 || pheno_subset[k] == 3) {
                                k0++;
                                omega_inv_tmp = omega_inv_tmp + u[i0][k0]*u[j0][k0]/eigen[k0];
                            }
                        }
                        family->a_store[1][i] = family->a_store[1][i] + omega_inv_tmp*trait_adj[j0];
                    }
				}
			}
		}


		/* free matrices */
		free_dmatrix(omega,1,n_pheno_phenorel_w_genorel,1,n_pheno_phenorel_w_genorel);
		free_dmatrix(phi_rr,1,n_pheno_phenorel_w_genorel,1,n_pheno_phenorel_w_genorel);
		free_dvector(trait_adj,1,n_pheno_phenorel_w_genorel);
        free_dmatrix(u,1,n_pheno_phenorel_w_genorel,1,n_pheno_phenorel_w_genorel);
        free_dvector(eigen,1,n_pheno_phenorel_w_genorel);

    } /* end subset = 3 */

	if (subset != 0 && subset != 1 && subset != 2 && subset !=3) {
		printf("\nERROR\nInvalid subset value, should be 0, 1, 2, or 3\n");
		exit(1);
	}

} /* end function a_vector */

/* function to get the a_all vector, a_all = Omega^-1(D-W*beta) if person is phenotyped
using Omega = I (the idenity matrix)
 */
void a_vector_i(struct FAMILY *family, struct DATA_STRUCT data_struct, struct MLE_RESULTS mle_res, int subset) {

	/* pull out variables from the struct family to use internally in the function */
	int n_total = family->n_total;
	int n_pheno = family->n_pheno;
	int *pheno_typed = family->pheno_typed;
	double **cov = family->cov;
	double *trait = family->trait;
	double **phi = family->phi;
	int n_cov = data_struct.n_cov;
	double *beta = mle_res.beta_ols;
	int *pheno_subset = family->pheno_subset;
	int i, j;

	/* create a_all = I*trait_adj = D-W*beta */

    /*	double tmp_adj;
	for(i=1; i<=n_pheno; i++) {
		tmp_adj = 0.0;
		for(j=1; j<= n_cov; j++) {
			tmp_adj = tmp_adj + cov[i][j]*beta[j];
		}
		family->a_store[3][i] = trait[i] - tmp_adj;
	}*/



    double tmp_adj;

    if (subset == 0) {
        for (i=1; i<=n_pheno; i++) {
            family->a_store[3][i] = 0.0;
            tmp_adj = 0.0;
            for(j=1; j<= n_cov; j++) {
                tmp_adj = tmp_adj + cov[i][j]*beta[j];
            }
            family->a_store[3][i] = trait[i] - tmp_adj;
        }
	} /* end subset = 0 */

    if (subset == 1) {
        for (i=1; i<=n_pheno; i++) {
            family->a_store[3][i] = 0.0;
            tmp_adj = 0.0;
            if (pheno_subset[i] == 1) {
                for(j=1; j<= n_cov; j++) {
                    tmp_adj = tmp_adj + cov[i][j]*beta[j];
                }
                family->a_store[3][i] = trait[i] - tmp_adj;
            }
        }
	} /* end subset = 2 */

    if (subset == 2) {
        for (i=1; i<=n_pheno; i++) {
            family->a_store[3][i] = 0.0;
            tmp_adj = 0.0;
            if (pheno_subset[i] == 1 || pheno_subset[i] == 2) {
                for(j=1; j<= n_cov; j++) {
                    tmp_adj = tmp_adj + cov[i][j]*beta[j];
                }
                family->a_store[3][i] = trait[i] - tmp_adj;
            }
        }
	} /* end subset = 2 */

    if (subset == 3) {
        for (i=1; i<=n_pheno; i++) {
            family->a_store[3][i] = 0.0;
            tmp_adj = 0.0;
            if (pheno_subset[i] == 1 || pheno_subset[i] == 2 || pheno_subset[i] == 3) {
                for(j=1; j<= n_cov; j++) {
                    tmp_adj = tmp_adj + cov[i][j]*beta[j];
                }
                family->a_store[3][i] = trait[i] - tmp_adj;
            }
        }
	} /* end subset = 3 */

} /* end function a_vector_i */

/* aug_phi_nn reads in one family at a time and fills in the values necessary to calculate
 the statistic and allele frequency for each marker, before the first family is processed,
 it is assumed that the components of MARKER_STATS stat have been initialized and set to zero, then
 for each family the components are updated
 */
void aug_phi_nn(int i_marker, int i_fam, struct FAMILY family, struct MARKER_STATS *stat) {

	/* pull out variables from the structures family to use internally in the function */
	int n_pheno = family.n_pheno;
	int n_geno = family.n_geno;
	int n_geno_pheno = family.n_geno_pheno;
	int n_geno_notpheno = family.n_geno_notpheno;
	int n_notgeno_pheno = family.n_notgeno_pheno;
	int *pheno_typed = family.pheno_typed;
	int **geno_pheno_typed = family.geno_pheno_typed;
	int **notgeno_pheno_typed = family.notgeno_pheno_typed;
	int *geno_notpheno_typed = family.geno_notpheno_typed;
	double **phi = family.phi;
	double *y = family.y;
	double *a_all = family.a_all;

	double **phi_nn;
	phi_nn = dmatrix(1,n_geno,1,n_geno);

	int i,j;

	/* fill in the entries for phi_nn such that entries for those geno&pheno come first and then
	 the entries for those geno&notpheno */
	for(i=1; i<=n_geno_pheno; i++) {
		for(j=1; j<=n_geno_pheno; j++) {
			phi_nn[i][j] = phi[geno_pheno_typed[i][1]][geno_pheno_typed[j][1]];
		}
		for(j=1; j<=n_geno_notpheno; j++) {
			phi_nn[i][j+n_geno_pheno] = phi[geno_pheno_typed[i][1]][geno_notpheno_typed[j]];
		}
	}

	for(i=1; i<=n_geno_notpheno; i++) {
		for(j=1; j<=n_geno_notpheno; j++) {
			phi_nn[i+n_geno_pheno][j+n_geno_pheno] = phi[geno_notpheno_typed[i]][geno_notpheno_typed[j]];
		}
		for(j=1; j<=n_geno_pheno; j++) {
			phi_nn[i+n_geno_pheno][j] = phi[geno_notpheno_typed[i]][geno_pheno_typed[j][1]];
		}
	}


	/* if there are individuals with observed pheno/cov but not genotyped then create v_m = phi_nm*a_m, where
	 the entries ordered the same way as in phi_nn above */
	double *v_m;
	v_m = dvector(1,n_geno);
	if(n_notgeno_pheno != 0) {
		for(i=1; i<=n_geno_pheno; i++) {
			v_m[i] = 0.0;
			for(j=1; j<=n_notgeno_pheno; j++) {
				v_m[i] = v_m[i] + phi[geno_pheno_typed[i][1]][notgeno_pheno_typed[j][1]]*a_all[notgeno_pheno_typed[j][2]];
			}
		}
		for(i=1; i<=n_geno_notpheno; i++) {
			v_m[i+n_geno_pheno] = 0.0;
			for(j=1; j<=n_notgeno_pheno; j++) {
				v_m[i+n_geno_pheno] = v_m[i+n_geno_pheno] + phi[geno_notpheno_typed[i]][notgeno_pheno_typed[j][1]]*a_all[notgeno_pheno_typed[j][2]];
			}
		}
	}

	/* create aug by merging 1_n, y, v_m, run cholesky
	 (if n_notgeno_pheno = 0 then do not try to add v_m, o.w. add v_m) */
	double **aug;
	double **chol, **cholaug;

	if(n_notgeno_pheno != 0) {
		aug = dmatrix(1,n_geno,1,3);
		for(i=1; i<=n_geno; i++) {
			aug[i][1] = 1.0;
			aug[i][2] = y[i];
			aug[i][3] = v_m[i];
		}

		/* run the cholesky */
		chol = dmatrix(1,n_geno,1,n_geno);
		cholaug = dmatrix(1,n_geno,1,3);

		/* initalize values with zero */
		for (i=1; i<=n_geno; i++) {
			for (j=1; j<=n_geno; j++) {
				chol[i][j] = 0.0;
			}
			cholaug[i][1] = 0.0;
			cholaug[i][2] = 0.0;
			cholaug[i][3] = 0.0;
		}

        int posdef=1;
        posdef = cholesky(phi_nn, n_geno, aug, 3, chol, cholaug, 1);
        if(posdef == 0) {
            printf("\nERROR: Exiting program!\n"
                   "The matrix Phi_NN is not positive semi-definite\n"
                   "Likely this is due to both individuals from a MZ twin pair"
                   "being genotyped\n");
            exit(1);
        }


	} else {
		aug = dmatrix(1,n_geno,1,2);
		for(i=1; i<=n_geno; i++) {
			aug[i][1] = 1.0;
			aug[i][2] = y[i];
		}

		/* run the cholesky */
		chol = dmatrix(1,n_geno,1,n_geno);
		cholaug = dmatrix(1,n_geno,1,2);

		/* initalize values with zero */
		for (i=1; i<=n_geno; i++) {
			for(j=1; j<= n_geno; j++) {
				chol[i][j] = 0.0;
			}
			cholaug[i][1] = 0.0;
			cholaug[i][2] = 0.0;
		}

        int posdef=1;
        posdef = cholesky(phi_nn, n_geno, aug, 2, chol, cholaug, 1);
        if(posdef == 0) {
            printf("\nERROR: Exiting program!\n"
                   "The matrix Phi_NN is not positive semi-definite\n"
                   "Likely this is due to both individuals from a MZ twin pair"
                   "being genotyped\n");
            exit(1);
        }


	}

	/* create
	 dot_1 = gamma_1^T*gamma_1
	 dot_1_y = gamma_1^T*gamma_y
	 dot_m = gamma_m^T*gamma_m
	 dot_1_m = gamma_1^T*gamma_m
	 dot_y_m = gamma_y^T*gamma_m
	 dot_y = gamma_y^T*gamma_y^T

	 dot_a_delta = a_n^T*delta_r
	 a_phi_a = a_m^T*phi_nm^T*a_n
	 dot_1_a = 1^T*a_n
	 dot_y_a = y^T*a_n
	 */

	double dot_1, dot_1_y, dot_m, dot_1_m, dot_y_m, dot_y, dot_a_delta, a_phi_a, dot_1_a, dot_y_a;
	dot_1 = 0.0;
	dot_1_y = 0.0;
	dot_m = 0.0;
	dot_1_m = 0.0;
	dot_y_m = 0.0;
	dot_y = 0.0;
	dot_a_delta = 0.0;
	a_phi_a = 0.0;
	dot_1_a = 0.0;
	dot_y_a = 0.0;


	if(n_notgeno_pheno != 0) {
		/* dot_1, dot_1_y dot_m, dot_1_m, dot_y_m */
		for(i=1; i<=n_geno; i++) {
			dot_1 = dot_1 + cholaug[i][1]*cholaug[i][1];
			dot_1_y = dot_1_y + cholaug[i][1]*cholaug[i][2];
			dot_m = dot_m + cholaug[i][3]*cholaug[i][3];
			dot_1_m = dot_1_m + cholaug[i][1]*cholaug[i][3];
			dot_y_m = dot_y_m + cholaug[i][2]*cholaug[i][3];
			dot_y = dot_y + cholaug[i][2]*cholaug[i][2];
		}

		/* dot_a_delta */
		double delta_r_tmp;
		for (i=1; i<= n_geno_pheno; i++) {
			delta_r_tmp = 0.0;
			for (j=1; j<=n_pheno; j++) {
				delta_r_tmp = delta_r_tmp + phi[geno_pheno_typed[i][1]][pheno_typed[j]]*a_all[j];
			}
			dot_a_delta = dot_a_delta + delta_r_tmp*a_all[geno_pheno_typed[i][2]];
		}

		/* a_phi_a */
		double a_phi_a_tmp;
		for (i=1; i<= n_geno_pheno; i++) {
			a_phi_a_tmp = 0.0;
			for(j=1; j<= n_notgeno_pheno; j++) {
				a_phi_a_tmp = a_phi_a_tmp + phi[geno_pheno_typed[i][1]][notgeno_pheno_typed[j][1]]*a_all[notgeno_pheno_typed[j][2]];
			}
			a_phi_a = a_phi_a + a_phi_a_tmp*a_all[geno_pheno_typed[i][2]];
		}

		/* dot_1_a, dot_y_a */
		for (i=1; i<=n_geno_pheno; i++) {
			dot_1_a = dot_1_a + a_all[geno_pheno_typed[i][2]];
			dot_y_a = dot_y_a + y[i]*a_all[geno_pheno_typed[i][2]];
		}

        printf("dot_1 %lf\n",dot_1);fflush(stdout);
        printf("dot_1_y %lf\n",dot_1_y);fflush(stdout);
        printf("dot_1_a %lf\n",dot_1_a);fflush(stdout);
        printf("dot_1_m %lf\n",dot_1_m);fflush(stdout);
        printf("dot_y %lf\n",dot_y);fflush(stdout);
        printf("dot_y_a %lf\n",dot_y_a);fflush(stdout);
        printf("dot_y_m %lf\n",dot_y_m);fflush(stdout);
        printf("dot_a_delta %lf\n",dot_a_delta);fflush(stdout);
        printf("dot_m %lf\n",dot_m);fflush(stdout);
        printf("a_phi_a %lf\n",a_phi_a);fflush(stdout);
        printf("n_geno %d\n",n_geno);fflush(stdout);

		/* update the structure stat */
		stat->v_1 += dot_1;
		stat->v_2 += dot_1_y;
		stat->v_3 += dot_a_delta + a_phi_a + dot_m;
		stat->v_4 += dot_1_a + dot_1_m;
		stat->v_5 += dot_y_a + dot_y_m;
		stat->v_6 += dot_y;
		stat->n += n_geno;
	} else {
		/* dot_1, dot_1_y */
		for(i=1; i<=n_geno; i++) {
			dot_1 = dot_1 + cholaug[i][1]*cholaug[i][1];
			dot_1_y = dot_1_y + cholaug[i][1]*cholaug[i][2];
			dot_y = dot_y + cholaug[i][2]*cholaug[i][2];
		}

		/* dot_a_delta */
		double delta_r_tmp;
		for (i=1; i<= n_geno_pheno; i++) {
			delta_r_tmp = 0.0;
			for (j=1; j<=n_pheno; j++) {
				delta_r_tmp = delta_r_tmp + phi[geno_pheno_typed[i][1]][pheno_typed[j]]*a_all[j];
			}
			dot_a_delta = dot_a_delta + delta_r_tmp*a_all[geno_pheno_typed[i][2]];
		}

		/* dot_1_a, dot_y_a */
		for (i=1; i<=n_geno_pheno; i++) {
			dot_1_a = dot_1_a + a_all[geno_pheno_typed[i][2]];
			dot_y_a = dot_y_a + y[i]*a_all[geno_pheno_typed[i][2]];
		}

		/* update the structure stat */
		stat->v_1 += dot_1;
		stat->v_2 += dot_1_y;
		stat->v_3 += dot_a_delta;
		stat->v_4 += dot_1_a;
		stat->v_5 += dot_y_a;
		stat->v_6 += dot_y;
		stat->n += n_geno;
	}



	/* free internal/temporary matrices/vectors */
	free_dmatrix(phi_nn,1,n_geno,1,n_geno);
	free_dvector(v_m,1,n_geno);
	free_dmatrix(chol,1,n_geno,1,n_geno);

	if(n_notgeno_pheno != 0) {
		free_dmatrix(aug,1,n_geno,1,3);
		free_dmatrix(cholaug,1,n_geno,1,3);
	} else {
		free_dmatrix(aug,1,n_geno,1,2);
		free_dmatrix(cholaug,1,n_geno,1,2);
	}

} /* end function aug_phi_nn */


/* svd_phi_nn reads in one family at a time and fills in the values necessary to calculate
 the statistic and allele frequency for each marker, before the first family is processed,
 it is assumed that the components of MARKER_STATS stat have been initialized and set to zero, then
 for each family the components are updated

 we can use this instead of aug_phi_nn in case the Cholesky does not work, i.e. when Phi_NN is not p.s.d.
 */
void svd_phi_nn(int i_marker, int i_fam, struct FAMILY family, struct MARKER_STATS *stat) {

	/* pull out variables from the structures family to use internally in the function */
	int n_pheno = family.n_pheno;
	int n_geno = family.n_geno;
	int n_geno_pheno = family.n_geno_pheno;
	int n_geno_notpheno = family.n_geno_notpheno;
	int n_notgeno_pheno = family.n_notgeno_pheno;
	int *pheno_typed = family.pheno_typed;
	int **geno_pheno_typed = family.geno_pheno_typed;
	int **notgeno_pheno_typed = family.notgeno_pheno_typed;
	int *geno_notpheno_typed = family.geno_notpheno_typed;
	double **phi = family.phi;
	double *y = family.y;
	double *a_all = family.a_all;

	double **phi_nn;
	phi_nn = dmatrix(1,n_geno,1,n_geno);

	int i,j;

	/* fill in the entries for phi_nn such that entries for those geno&pheno come first and then
	 the entries for those geno&notpheno */
	for(i=1; i<=n_geno_pheno; i++) {
		for(j=1; j<=n_geno_pheno; j++) {
			phi_nn[i][j] = phi[geno_pheno_typed[i][1]][geno_pheno_typed[j][1]];
		}
		for(j=1; j<=n_geno_notpheno; j++) {
			phi_nn[i][j+n_geno_pheno] = phi[geno_pheno_typed[i][1]][geno_notpheno_typed[j]];
		}
	}

	for(i=1; i<=n_geno_notpheno; i++) {
		for(j=1; j<=n_geno_notpheno; j++) {
			phi_nn[i+n_geno_pheno][j+n_geno_pheno] = phi[geno_notpheno_typed[i]][geno_notpheno_typed[j]];
		}
		for(j=1; j<=n_geno_pheno; j++) {
			phi_nn[i+n_geno_pheno][j] = phi[geno_notpheno_typed[i]][geno_pheno_typed[j][1]];
		}
	}


	/* if there are individuals with observed pheno/cov but not genotyped then create v_m = phi_nm*a_m, where
	 the entries ordered the same way as in phi_nn above */
	double *v_m;
	v_m = dvector(1,n_geno);
	if(n_notgeno_pheno != 0) {
		for(i=1; i<=n_geno_pheno; i++) {
			v_m[i] = 0.0;
			for(j=1; j<=n_notgeno_pheno; j++) {
				v_m[i] = v_m[i] + phi[geno_pheno_typed[i][1]][notgeno_pheno_typed[j][1]]*a_all[notgeno_pheno_typed[j][2]];
			}
		}
		for(i=1; i<=n_geno_notpheno; i++) {
			v_m[i+n_geno_pheno] = 0.0;
			for(j=1; j<=n_notgeno_pheno; j++) {
				v_m[i+n_geno_pheno] = v_m[i+n_geno_pheno] + phi[geno_notpheno_typed[i]][notgeno_pheno_typed[j][1]]*a_all[notgeno_pheno_typed[j][2]];
			}
		}
	}

    /* singular value decompostion of phi_nn */
    double *eigen; eigen = dvector(1,n_geno);
    double **u; u = dmatrix(1,n_geno,1,n_geno);

    svdcomp(phi_nn, n_geno, n_geno, eigen, u);

    /* create:
        u_1 = 1^T x u
        u_y = y^T x u
        u_m = v_m^T x u (phi_nm x a_m)^T x u
     */

    double *u_1, *u_y, *u_m;
    u_1=dvector(1,n_geno);
    u_y=dvector(1,n_geno);
    u_m=dvector(1,n_geno);

    for(i=1; i<=n_geno; i++) {
        u_1[i] = 0;
        u_y[i] = 0;
        u_m[i] = 0;

        for(j=1; j<=n_geno; j++) {
            u_1[i] = u_1[i] + u[j][i];   // u_1[i] is the sum of the ith column u (dot product with vector of ones)
            u_y[i] = u_y[i] + y[j]*u[j][i]; // u_y[i] is the dot product of y and the ith column of u, y's order geno_pheno then geno_notpheno
            u_m[i] = u_m[i] + v_m[j]*u[j][i]; // u_m[i] is the dot prodcut of v_m and the ith column of u
        }
    }

    /* create:
        dot_1   = 1 x Phi_nn^-1 x 1^T       = u_1 x eigen^-1 x u_1^T
        dot_1_y = 1 x Phi_nn^-1 x y^T       = u_1 x eigen^-1 x u_y^T
        dot_m   = v_m x Phi_nn^-1 x v_m^T   = u_m x eigen^-1 x u_m^T
        dot_1_m = 1 x Phi_nn^-1 x v_m^T     = u_1 x eigen^-1 x u_m^T
        dot_y_m = y x Phi_nn^-1 x v_m^T     = u_y x eigen^-1 x u_m^T
        dot_y   = y x Phi_nn^-1 x y^T       = u_y x eigen^-1 x u_y^T

        dot_a_delta = a_n^T*delta_r
        a_phi_a = a_m^T*phi_nm^T*a_n
        dot_1_a = 1^T*a_n
        dot_y_a = y^T*a_n

     */

    double dot_1, dot_1_y, dot_m, dot_1_m, dot_y_m, dot_y, dot_a_delta, a_phi_a, dot_1_a, dot_y_a;
	dot_1 = 0.0;
	dot_1_y = 0.0;
	dot_m = 0.0;
	dot_1_m = 0.0;
	dot_y_m = 0.0;
	dot_y = 0.0;
	dot_a_delta = 0.0;
	a_phi_a = 0.0;
	dot_1_a = 0.0;
	dot_y_a = 0.0;

    if(n_notgeno_pheno != 0) {
		/* dot_1, dot_1_y dot_m, dot_1_m, dot_y_m */
		for(i=1; i<=n_geno; i++) {
			dot_1 = dot_1 + u_1[i]*u_1[i]/eigen[i];
			dot_1_y = dot_1_y + u_1[i]*u_y[i]/eigen[i];
			dot_m = dot_m + u_m[i]*u_m[i]/eigen[i];
			dot_1_m = dot_1_m + u_1[i]*u_m[i]/eigen[i];
			dot_y_m = dot_y_m + u_y[i]*u_m[i]/eigen[i];
			dot_y = dot_y + u_y[i]*u_y[i]/eigen[i];
		}

		/* dot_a_delta */
		double delta_r_tmp;
		for (i=1; i<= n_geno_pheno; i++) {
			delta_r_tmp = 0.0;
			for (j=1; j<=n_pheno; j++) {
				delta_r_tmp = delta_r_tmp + phi[geno_pheno_typed[i][1]][pheno_typed[j]]*a_all[j];
			}
			dot_a_delta = dot_a_delta + delta_r_tmp*a_all[geno_pheno_typed[i][2]];
		}

		/* a_phi_a */
		double a_phi_a_tmp;
		for (i=1; i<= n_geno_pheno; i++) {
			a_phi_a_tmp = 0.0;
			for(j=1; j<= n_notgeno_pheno; j++) {
				a_phi_a_tmp = a_phi_a_tmp + phi[geno_pheno_typed[i][1]][notgeno_pheno_typed[j][1]]*a_all[notgeno_pheno_typed[j][2]];
			}
			a_phi_a = a_phi_a + a_phi_a_tmp*a_all[geno_pheno_typed[i][2]];
		}

		/* dot_1_a, dot_y_a */
		for (i=1; i<=n_geno_pheno; i++) {
			dot_1_a = dot_1_a + a_all[geno_pheno_typed[i][2]];
			dot_y_a = dot_y_a + y[i]*a_all[geno_pheno_typed[i][2]];
		}

		/* update the structure stat */
		stat->v_1 += dot_1;
		stat->v_2 += dot_1_y;
		stat->v_3 += dot_a_delta + a_phi_a + dot_m;
		stat->v_4 += dot_1_a + dot_1_m;
		stat->v_5 += dot_y_a + dot_y_m;
		stat->v_6 += dot_y;
		stat->n += n_geno;

	} else {
		/* dot_1, dot_1_y */
		for(i=1; i<=n_geno; i++) {
			dot_1 = dot_1 + u_1[i]*u_1[i]/eigen[i];
			dot_1_y = dot_1_y + u_1[i]*u_y[i]/eigen[i];
			dot_y = dot_y + u_y[i]*u_y[i]/eigen[i];
		}

		/* dot_a_delta */
		double delta_r_tmp;
		for (i=1; i<= n_geno_pheno; i++) {
			delta_r_tmp = 0.0;
			for (j=1; j<=n_pheno; j++) {
				delta_r_tmp = delta_r_tmp + phi[geno_pheno_typed[i][1]][pheno_typed[j]]*a_all[j];
			}
			dot_a_delta = dot_a_delta + delta_r_tmp*a_all[geno_pheno_typed[i][2]];
		}

		/* dot_1_a, dot_y_a */
		for (i=1; i<=n_geno_pheno; i++) {
			dot_1_a = dot_1_a + a_all[geno_pheno_typed[i][2]];
			dot_y_a = dot_y_a + y[i]*a_all[geno_pheno_typed[i][2]];
		}

		/* update the structure stat */
		stat->v_1 += dot_1;
		stat->v_2 += dot_1_y;
		stat->v_3 += dot_a_delta;
		stat->v_4 += dot_1_a;
		stat->v_5 += dot_y_a;
		stat->v_6 += dot_y;
		stat->n += n_geno;
	}

	/* free internal/temporary matrices/vectors */
	free_dmatrix(phi_nn,1,n_geno,1,n_geno);
	free_dvector(v_m,1,n_geno);
    free_dvector(eigen,1,n_geno);
    free_dmatrix(u,1,n_geno,1,n_geno);
    free_dvector(u_1,1,n_geno);
    free_dvector(u_y,1,n_geno);
    free_dvector(u_m,1,n_geno);


} /* end function svd_phi_nn */
