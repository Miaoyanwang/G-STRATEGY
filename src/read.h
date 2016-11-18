#ifndef READ_H
#define READ_H

#include <stdio.h>
#include "mastor.h"
#include "hashes.h"
#include <zlib.h>

void readcov(char *filename, struct DATA_STRUCT *data_struct);
void readstruct(char *filename, struct DATA_STRUCT *data_struct);
void readped(char *filename, struct FAMILY *family, struct HASH *hash, struct DATA_STRUCT data_struct);
void addgeno(int *geno, int allele, double *y);
void readkin(char *filename, struct FAMILY *family, struct HASH hash);
void marker_skrats_allocate(struct FAMILY *family, struct DATA_STRUCT data_struct);
void marker_skrats_initialize(struct FAMILY *family, struct DATA_STRUCT data_struct);
void marker_skrats_free(struct FAMILY *family, struct DATA_STRUCT data_struct);

void readmarker_header(gzFile markfile, char **id_geno2all, struct FAMILY *family, struct HASH hash, struct DATA_STRUCT data_struct);
void readmarker(gzFile markfile, struct FAMILY *family, struct HASH hash, char **id_geno2all, struct DATA_STRUCT data_struct, struct ANNOTATION *anno);

void fill_pheno_subset(struct FAMILY *family, struct DATA_STRUCT data_struct);

#endif
