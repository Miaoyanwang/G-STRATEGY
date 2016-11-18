#ifndef MASTOR_H
#define MASTOR_H

#define MAXLEN 1024
#define MISSVAL -9.0
#define NUMTOL 1e-6
#define ANALYSIS_STATUS 50000

struct DATA_STRUCT {
  int n_fam;    /* total number of families */
  int n_marker;    /* total number of markers */
  int n_total;    /* total number of individuals (across all families) */
  int n_typed;    /* total number of geno/phenotyped individuals (across all families) */
  int n_cov;    /* total number of covariates including intercept */
  double missing_val;
  double tol;
};

struct FAMILY {
  int n_total;            /* total number individuals in family */
  int n_pheno;            /* number of phenotyped individuals in family */
  double **phi;           /* phi[ind i][ind j] is the kinship matrix: i=j then 1+h_i, o.w. 2*phi_ij */
  int *pheno_typed;       /* pheno_typed[index among pheno/cov typed] = index among all */
  int *pheno_typed_inv;   /* pheno_typed_inv[index among all] = index among pheno/cov typed, 0 o.w. */
  double **cov;           /* cov[index among phenotyped][x] = cov value for ind for cov x */
  double *trait;          /* trait[index among phenotyped] */
  double *a_all;          /* alpha vector, a_all[index among phenotyped] */
  double **a_store;       /* alpha matrix for all types: a[1][pheno index] = estimed var, a[2][] = true var, a[3][] = omega=I */
  int n_geno;             /* number of individuals genotyped at current marker */
  int n_geno_pheno;       /* number of individuals genotyped and pheno/cov typed at current marker */
  int n_geno_notpheno;    /* number of individuals genotyped but not pheno/cov typed at current marker */
  int n_notgeno_pheno;    /* number of individuals not genotyped but pheno/cov typed at current marker */
  int **geno_pheno_typed; /* geno_pheno_typed[index among geno/pheno typed][1] = index among all
                             geno_pheno_typed[index among geno/pheno typed][2] = index among pheno */
  int **notgeno_pheno_typed;  /* notgeno_pheno_typed[index among notgeno/pheno typed][1] = index among all
                                 notgeno_pheno_typed[index among notgeno/pheno typed][2] = index among pheno */
  int *geno_notpheno_typed;   /* geno_notpheno[index among geno/notpheno] = index among all */
  double *y;              /* y[index] = 0.5x(no. of alleles 1);
                              index = 1...n_geno_pheno, 1+n_geno_pheno...n_geno_notpheno+n_geno_pheno */
  int *geno_typed;        /* geno_typed[index among all] = 0 (not genotyped) or 1 (genotyped) NOTE only for the readmarker_header */
  int *pheno_subset;      /* pheno_subset[index among pheno] = 1 (genotyped), 2 (not geno, genotyped relative),
                                 3 (not geno, phenotyped relative with genotyped relative) NOTE only for the fill_pheno_subset */
};


struct FAMILY_SVD {
  double *trait_u;
  double **cov_u;
  double *lambda;
  int n_ind;
};

struct OPT_STRUCT {
  struct FAMILY *family;
  struct FAMILY_SVD *svd;
  struct DATA_STRUCT data_struct;
};

struct MARKER_STATS {
  double v_1;
  double v_2;
  double v_3;
  double v_4;
  double v_5;
  double v_6;
  int n;
};

struct ANNOTATION {
  char *chr;
  char *markername;
  int centimorgan;
  int basepair;
};

struct RESULTS {
  double p_0_blue;
  double mqlsq0, mqlsq1;
  double p_value0, p_value1;
};

struct MLE_PARAM {
  double ax,bx,cx,tol;
};


struct MLE_RESULTS {
  double alpha_hat;
  double err_var_hat,add_var_hat,herit,*beta_hat, *beta_ols;
  double var_err_var_hat,var_add_var_hat,var_herit,*var_beta_hat;
};

struct GTAM {
  double dd_est;
  double *d_est;
  double **m_est;
  double t_gtam_est, p_value_est;
  int n_typed_all;
};

#endif
