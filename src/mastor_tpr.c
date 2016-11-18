#include "mastor.h"

#include "nrutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>

#include "brak.h"
#include "brent.h"
#include "calc.h"
#include "cholesky.h"
#include "mle.h"
#include "svdcomp.h"
#include "read.h"
#include "hashes.h"
#include "datacheck.h"
#include "tpr.h"



struct OPT_STRUCT f_opt_struct;

void init_opt_struct(struct FAMILY *family, struct FAMILY_SVD *svd, struct DATA_STRUCT data_struct) {
  f_opt_struct.family = family;
  f_opt_struct.svd = svd;
  f_opt_struct.data_struct = data_struct;
}

int main (int argc, char **argv) {

  struct DATA_STRUCT data_struct;
  struct FAMILY *family;
  struct MARKER_STATS stat_marker;
  struct RESULTS results;
  struct FAMILY_SVD *svd;
  struct MLE_RESULTS mle_res;
  struct MLE_PARAM mle_initial;
  struct ANNOTATION anno;

  int i,j,k;
  double alpha;
  int n_pheno;
  int arg;

  double like_from, like_to, like_by;

  struct HASH hash;
  /*hash = (struct HASH) malloc(sizeof(struct HASH));*/
  hash.fam2fam = kh_init(str);
  hash.ind2fam = kh_init(str);
  hash.ind2ind = kh_init(str);

  char *pedfile, *genofile, *kinfile;
  pedfile = (char*) malloc(MAXLEN);
  genofile = (char*) malloc(MAXLEN);
  kinfile = (char*) malloc(MAXLEN);

  char *prefix;
  prefix = (char*) malloc(MAXLEN);
  strcpy(prefix,"");

    char *mle_outfile;
    mle_outfile = (char*) malloc(MAXLEN);
    strcpy(mle_outfile, "");
    
  strcpy(pedfile,"ped.txt");
  strcpy(genofile,"geno.txt");
  strcpy(kinfile,"kin.txt");

  int pfile=0, gfile=0, kfile=0;
  int likeplot = 0;
  int nullmle = 0;
  int option = 0;
  int subset = -9;
  int tpr = 0;

  const char* opt_string = "p:g:k:a:m:s:";
  static struct option long_options[] =
  {
    {"pedfile",required_argument,0,'p'},
    {"genofile", required_argument,0,'g'},
    {"kinfile",required_argument,0,'k'},
    {"tpr",no_argument,0,0},
    {0,0,0,0}
  };

  int option_index = 0, c;
  while (1) {

    c = getopt_long(argc,argv,opt_string,long_options, &option_index);
    if(c == -1) {
      break;
    }

    switch(c) {
      case 'p':
        strcpy(pedfile, optarg);
        pfile = 1;
        break;
      case 'g':
        strcpy(genofile, optarg);
        gfile = 1;
        break;
      case 'k':
        strcpy(kinfile, optarg);
        kfile = 1;
        break;
     case 0:
        if (strcmp("tpr",long_options[option_index].name) == 0) {
          tpr = 1;
          subset = 0;
          option = 1;
          nullmle = 0;
          break;
        }
        break;
      default:
        exit(1);
    } // close switch
  } // close while

    
    
  data_struct.tol = NUMTOL;
  data_struct.missing_val = MISSVAL;

    
  format_checks(pedfile,genofile,kinfile,&data_struct);
    


    
  int n_fam = data_struct.n_fam;
  int n_total = data_struct.n_total;
  int n_typed = data_struct.n_typed;

  char **id_geno2all = (char**)malloc((n_typed+1) * sizeof(char*)); /* id_geno2all[unique index id over all genotyped ind across all ped] = unique index id over all ind across all ped*/

  int n_cov = data_struct.n_cov;
  int unr;

  FILE *out_param;
  FILE *outfile;
  gzFile markfile;

  anno.chr = (char*) malloc(4);
  anno.markername = (char*) malloc(MAXLEN);

  /* ====================================================================================== */
  /* this part is for MASTOR analysis only (not GTAM analysis)                              */
  /* ====================================================================================== */

    /* option 1 = run MASTOR */
    if(option == 1) {


      family = (struct FAMILY *)malloc(sizeof(struct FAMILY)*(n_fam+1));
      svd = (struct FAMILY_SVD *)malloc(sizeof(struct FAMILY_SVD)*(n_fam+1));

      readped(pedfile,family, &hash, data_struct);
      readkin(kinfile,family, hash);
      markfile = gzopen(genofile,"r");
      marker_skrats_allocate(family, data_struct);

      if(markfile == NULL) {
        printf("Can't open file with marker genotypes\n");
        exit(1);
      }

      /* void readmarker_header(FILE *markfile, char **id_geno2all, struct FAMILY *family, struct HASH hash, struct DATA_STRUCT data_struct); */
      readmarker_header(markfile, id_geno2all, family, hash, data_struct);

      /* void fill_pheno_subset(struct FAMILY *family, struct DATA_STRUCT data_struct); */
      fill_pheno_subset(family, data_struct);

      /* int unrelated_check (struct FAMILY *family, struct DATA_STRUCT data_struct, int subset, int relpair); */
      unr = unrelated_check(family, data_struct, subset, 1);

      mle_initial.tol = NUMTOL;
      initial_ax(family, &mle_initial, data_struct);

      for(i=1; i<= n_fam;i++) {
        /* void fill_struct_svd(struct FAMILY family, struct DATA_STRUCT data_struct, struct FAMILY_SVD *svd, int subset); */
        fill_struct_svd(family[i], data_struct, &svd[i],subset);
      }
      init_opt_struct(family, svd, data_struct);

      /* if sample contains enough phenotyped relative pairs */
      if(unr==0) {
        initial_brak(&mle_initial, &f_neg_log_like_mle);

        alpha = estimate_alpha_hat(mle_initial);
        mle_res.alpha_hat = alpha;
        mle_res.var_beta_hat = dvector(1,data_struct.n_cov);
        mle_res.beta_hat = dvector(1,data_struct.n_cov);
        mle_res.beta_ols = dvector(1,data_struct.n_cov);

        /* if additive variance < 0 try to include more individuals in the estimation */
        if (alpha < 0.0) {
          for(i=1; i<= n_fam;i++) {
            /* void fill_struct_svd(struct FAMILY family, struct DATA_STRUCT data_struct, struct FAMILY_SVD *svd, int subset);
             subset = 0 includes all phenotyped individuals in the estimation */
            fill_struct_svd(family[i], data_struct, &svd[i],0);
          }
          init_opt_struct(family, svd, data_struct);
          initial_brak(&mle_initial, &f_neg_log_like_mle);

          alpha = estimate_alpha_hat(mle_initial);
          mle_res.alpha_hat = alpha;
          mle_res.var_beta_hat = dvector(1,data_struct.n_cov);
          mle_res.beta_hat = dvector(1,data_struct.n_cov);
          mle_res.beta_ols = dvector(1,data_struct.n_cov);

        } /* done including more individuals */

        /* if additive variance is still < 0 after including additional individuals in the estimation.
            Then use OLS to replace the GLS. Use the original smaller subset
                else proceed with using the GLS and calcualte the OLS but not replace GLS with OLS */
        /* void estimate_mle_ols (struct MLE_RESULTS *mle_res, struct FAMILY *family, struct DATA_STRUCT data_struct, int subset, int alpha_neg); */

        if (alpha < 0.0) {
            mle_res.alpha_hat = 0.0;
            estimate_mle_ols (&mle_res, family, data_struct, subset, 1);
        } else {
          estimate_mle (&mle_res, svd, data_struct);
          estimate_mle_ols (&mle_res, family, data_struct, subset, 0);
        }
      }


      /* if sample does not contain enough (or none) phenotyped relative pairs */
      if(unr==1) {
        mle_res.alpha_hat = 0.0;
        mle_res.var_beta_hat = dvector(1,data_struct.n_cov);
        mle_res.beta_hat = dvector(1,data_struct.n_cov);
        mle_res.beta_ols = dvector(1,data_struct.n_cov);
        estimate_mle_ols (&mle_res, family, data_struct, 2, 1);

        printf("WARNING: Not enough relative pairs for VC estimation, using OLS estimators for betas\n");

        if(likeplot==1) {
          printf("WARNING: User selected option to calculate various values of likelihood\n"
                 "Calculations were not performed because sample did not contain any\n"
                 "relative pairs to estimate the additive variance, thus only the error\n"
                 "variance is estimated using ordinary least squares\n");
        }
      } // end if(unr==1)

        
      /* output TPR (a_vector) = transformed phenotypic residuals */
      if(tpr == 1) {
        char *tpr_outfile;
        tpr_outfile = (char*) malloc(MAXLEN);
        strcpy(tpr_outfile, "");
        strcpy(tpr_outfile,prefix);
        if(strcmp("",prefix) == 0) {
          strcat(tpr_outfile,"TPR.txt");
        } else {
          char postfix[20];
          strcpy(postfix,"_TPR.txt");
          strcat(tpr_outfile,postfix);
        }
        for(i=1; i<= n_fam;i++) {
          /* void a_vector(struct FAMILY *family, struct DATA_STRUCT data_struct, struct MLE_RESULTS mle_res, int subset); */
          a_vector(&family[i], data_struct, mle_res,subset);
        }
        print_tpr(tpr_outfile,family, hash, data_struct);
        free(tpr_outfile);
     } // end if(tpr==1)

    }
    
    for(i=1; i<= n_fam;i++) {
      n_pheno = family[i].n_pheno;
      if (svd[i].n_ind != 0) {
        free_dvector(svd[i].trait_u,1,n_pheno);
        free_dmatrix(svd[i].cov_u,1,n_pheno,1,n_cov);
        free_dvector(svd[i].lambda,1,n_pheno);
      }

      free_dmatrix(family[i].phi,1,n_total,1,n_total);
      free_ivector(family[i].geno_typed,1,n_total);
      free_ivector(family[i].pheno_typed,1,n_pheno);
      free_ivector(family[i].pheno_typed_inv,1,n_total);
      free_dvector(family[i].trait,1,n_pheno);
      free_dmatrix(family[i].cov,1,n_pheno,1,n_cov);
      free_dmatrix(family[i].a_store,1,3,1,n_pheno);

    }

    free_dvector(mle_res.var_beta_hat,1,data_struct.n_cov);
    free_dvector(mle_res.beta_hat,1,data_struct.n_cov);
    if(unr==0){ free_dvector(mle_res.beta_ols,1,data_struct.n_cov); }

    free(family);
    free(svd);

   
  return 0;
}
