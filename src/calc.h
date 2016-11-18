#ifndef CALC_H
#define CALC_H

#include "mastor.h"

void a_vector(struct FAMILY *family, struct DATA_STRUCT data_struct, struct MLE_RESULTS mle_res, int subset);
void a_vector_i(struct FAMILY *family, struct DATA_STRUCT data_struct, struct MLE_RESULTS mle_res, int subset);
void aug_phi_nn(int i_marker, int i_fam, struct FAMILY family, struct MARKER_STATS *stat);
void svd_phi_nn(int i_marker, int i_fam, struct FAMILY family, struct MARKER_STATS *stat);
#endif
