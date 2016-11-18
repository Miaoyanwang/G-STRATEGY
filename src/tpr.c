#include "tpr.h"

void print_tpr(char *filename, struct FAMILY *family, struct HASH hash, struct DATA_STRUCT data_struct) {
  int n_fam = data_struct.n_fam;
  int n_typed = data_struct.n_typed;
  int n_total = data_struct.n_total;
  int fam_val, fam_val2,ind_val, f_k, i_k, i_k_2;
  const char *ind_id;
  const char *fam_id;
  FILE *outfile;
  outfile = fopen(filename,"w");

  for(f_k=kh_begin(hash.fam2fam);f_k!=kh_end(hash.fam2fam);f_k++) {
    if(kh_exist(hash.fam2fam,f_k)) {
      fam_id = kh_key(hash.fam2fam,f_k);
      fam_val = kh_value(hash.fam2fam,f_k);
      for(i_k=kh_begin(hash.ind2fam);i_k!=kh_end(hash.ind2fam);i_k++) {
        if(kh_exist(hash.ind2fam,i_k)) {
          ind_id = kh_key(hash.ind2fam,i_k);
          fam_val2 = kh_value(hash.ind2fam,i_k);
          if(fam_val == fam_val2) {
            i_k_2=kh_get(str,hash.ind2ind,ind_id);
            ind_val = kh_value(hash.ind2ind,i_k_2);
            if(family[fam_val].pheno_typed_inv[ind_val] !=0) {
              //fprintf(outfile,"%s %s %lf\n",fam_id,ind_id,family[fam_val].a_store[1][family[fam_val].pheno_typed_inv[ind_val]]);
              fprintf(outfile,"%s %s %lf %lf\n",fam_id,ind_id,family[fam_val].a_store[1][family[fam_val].pheno_typed_inv[ind_val]],family[fam_val].trait[family[fam_val].pheno_typed_inv[ind_val]]);
            } else {
              fprintf(outfile,"%s %s %d\n",fam_id,ind_id,0);
            }
          }
        }
      }
    }
  }

  fclose(outfile);

}
