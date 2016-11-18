#ifndef DATACHECK_H
#define DATACHECK_H

#include <math.h>
#include "nrutil.h"
#include "mastor.h"
#include "hashes.h"
#include <stdio.h>
#include <zlib.h>

int unrelated_check (struct FAMILY *family, struct DATA_STRUCT data_struct, int subset, int relpair);
void format_checks(char *pedfilename, char *markerfilename, char *kinfilename, struct DATA_STRUCT *data_struct);
int ped_famorder(khash_t(str) *fam_list, const char *fam);
int ped_uniqueind(khash_t(str) *ind_list, const char *ind );
int marker_uniqueind(khash_t(str) *ind_to_typed, const char *ind);

#endif
