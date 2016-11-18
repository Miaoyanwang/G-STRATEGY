#ifndef MLE_H
#define MLE_H

#include <math.h>
#include "nrutil.h"
#include "svdcomp.h"
#include "brak.h"
#include "brent.h"
#include "mastor.h"

void fill_struct_svd(struct FAMILY family, struct DATA_STRUCT data_struct, struct FAMILY_SVD *svd, int subset);
double calc_err_var_hat(double alpha, struct FAMILY_SVD *svd, struct DATA_STRUCT data_struct);
double neg_log_like_mle(double alpha, struct FAMILY_SVD *svd, struct DATA_STRUCT data_struct);
double f_neg_log_like_mle(double alpha);
double estimate_alpha_hat (struct MLE_PARAM param);
void estimate_mle (struct MLE_RESULTS *mle_res, struct FAMILY_SVD *svd, struct DATA_STRUCT data_struct);
void estimate_mle_ols (struct MLE_RESULTS *mle_res, struct FAMILY *family, struct DATA_STRUCT data_struct, int subset, int alpha_neg);

double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin);

#endif
