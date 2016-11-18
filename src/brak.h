#ifndef BRAK_H
#define BRAK_H

#include "mastor.h"

void brak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double));
void initial_ax(struct FAMILY *family, struct MLE_PARAM *mle_initial, struct DATA_STRUCT data_struct);
void initial_brak(struct MLE_PARAM *mle_initial, double (*func)(double));

#endif
