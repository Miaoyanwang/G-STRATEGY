#include "read.h"
#include "hashes.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nrutil.h"


/* function to read in the pedigree information, the trait, and covariates, assumes number of covarites
 is known,  that is the function format_checks function has been called */
void readped(char *filename, struct FAMILY *family, struct HASH *hash, struct DATA_STRUCT data_struct) {

	if (hash->fam2fam != NULL) {
		kh_destroy(str,hash->fam2fam);
		kh_destroy(str,hash->ind2fam);
		kh_destroy(str,hash->ind2ind);
	}

	hash->fam2fam = kh_init(str);
	hash->ind2fam = kh_init(str);
	hash->ind2ind = kh_init(str);


	FILE *pedfile;
	int n_cov = data_struct.n_cov;
	int ind, fam, fam_total = 0, n_total, n_pheno;
	int ret;
	khiter_t k;
	fam_total = 0;
	char *fam_id, *ind_id, *fa, *mo;
	fam_id = (char*) malloc(MAXLEN);
	ind_id = (char*) malloc(MAXLEN);
	fa = (char*) malloc(MAXLEN);
	mo = (char*) malloc(MAXLEN);

	int sex;
	double missing_val = data_struct.missing_val;
	double trait = missing_val;
	double cov_tmp = missing_val;
	double tol = data_struct.tol;
	int f_read=1, i_cov;

	double *cov;
	cov = dvector(1,n_cov);

	/* read through the file once to count n_total (total number of people in the pedigree)
	 and n_pheno (total number of individuals phenotyped, where phenotyped means that trait
	 and all covariates are non-missing).
	 It is assumed that no person has pariatial trait/covariate information, either all is
	 missing or nothing is missing. */
	pedfile = fopen(filename, "r");

	if(pedfile == NULL) {
		printf("Can't open pedigree and trait file\n");
		exit(1);
	}

	while(!feof(pedfile)) {

		f_read=1;
        if(n_cov>1){
            fscanf(pedfile,"%s %s %s %s %d %lf ", fam_id, ind_id, fa, mo, &sex, &trait);
        }
        if (n_cov==1) {
            fscanf(pedfile,"%s %s %s %s %d %lf \n", fam_id, ind_id, fa, mo, &sex, &trait);
        }

		if(fabs(trait-missing_val) < tol) {
			f_read = 0;
		}

        if (n_cov>1) {
            for(i_cov=2; i_cov<n_cov; i_cov++) {
                fscanf(pedfile,"%lf ", &cov_tmp);
                if(fabs(cov_tmp-missing_val) < tol) {
                    f_read = 0;
                }
            }

            fscanf(pedfile, "%lf \n", &cov_tmp);
            if(fabs(cov_tmp-missing_val) < tol) {
                f_read = 0;
            }
        }

		/* update family id */
		k = kh_get(str,hash->fam2fam,fam_id);
		if(k == kh_end(hash->fam2fam)) {
			fam_total = fam_total + 1;
			fam = fam_total;

			k = kh_put(str,hash->fam2fam, strdup(fam_id), &ret);
			kh_value(hash->fam2fam,k) = fam;

			k = kh_put(str,hash->ind2fam, strdup(ind_id), &ret);
			kh_value(hash->ind2fam,k) = fam;

			family[fam].n_total = 1;
			if (f_read != 0) {
				family[fam].n_pheno = 1;
			} else {
				family[fam].n_pheno = 0;
			}

			k = kh_put(str,hash->ind2ind, strdup(ind_id), &ret);
			kh_value(hash->ind2ind,k) = family[fam].n_total;

		} else {
			fam = kh_value(hash->fam2fam,k);

			k = kh_put(str,hash->ind2fam,strdup(ind_id), &ret);
			kh_value(hash->ind2fam,k) = fam;

			family[fam].n_total++;
			if (f_read != 0) {
				family[fam].n_pheno++;
			}

			k = kh_put(str,hash->ind2ind,strdup(ind_id),&ret);
			kh_value(hash->ind2ind,k) = family[fam].n_total;

		}

	}

	fclose(pedfile);


	/* for each family allocate memory for necessary vectors and matrices */
	for(k = kh_begin(hash->fam2fam); k != kh_end(hash->fam2fam); k++) {
		if (kh_exist(hash->fam2fam,k)) {
			fam = kh_value(hash->fam2fam,k);
			n_total = family[fam].n_total;
			n_pheno = family[fam].n_pheno;
			family[fam].phi = dmatrix(1,n_total,1,n_total);
            family[fam].geno_typed = ivector(1,n_total);
			family[fam].pheno_typed = ivector(1,n_pheno);
			family[fam].pheno_typed_inv = ivector(1,n_total);
			family[fam].trait = dvector(1,n_pheno);
			family[fam].cov = dmatrix(1,n_pheno,1,n_cov);
			family[fam].a_store = dmatrix(1,3,1,n_pheno);
            family[fam].pheno_subset = ivector(1,n_pheno);
		}
	}

	/* read through the file again to index individuals and families, and to save trait and covariates */
	pedfile = fopen(filename, "r");

	if(pedfile == NULL) {
		printf("Can't open pedigree and trait file\n");
		exit(1);
	}

	n_total = 0; /* temp variable for the n_total at the family currently read in */
	int i_total = 0; /* temp variable that counts how many have been read in of the current family */
	int i_pheno = 0; /* temp variable that coutns how many phenotyped have been read in of the current family */

	while(!(feof(pedfile))) {

		f_read=1;

        if(n_cov>1){
            fscanf(pedfile,"%s %s %s %s %d %lf ", fam_id, ind_id, fa, mo, &sex, &trait);
        }
        if (n_cov==1) {
            fscanf(pedfile,"%s %s %s %s %d %lf \n", fam_id, ind_id, fa, mo, &sex, &trait);
        }

		k = kh_get(str,hash->ind2fam,ind_id);
		fam = kh_value(hash->ind2fam,k);

		k = kh_get(str,hash->ind2ind,ind_id);
		ind = kh_value(hash->ind2ind,k);

		if(fabs(trait-missing_val) < tol) {
			f_read = 0;
		}

		cov[1] = 1.0; /* force the user to use an intercept, which is assumed to be absent in input file */

        if (n_cov>1) {
            for(i_cov=2; i_cov<n_cov; i_cov++) {
                fscanf(pedfile,"%lf ", &cov_tmp);
                if(fabs(cov_tmp-missing_val) < tol) {
                    f_read = 0;
                }
                cov[i_cov] = cov_tmp;
            }

            fscanf(pedfile, "%lf \n", &cov_tmp);
            if(fabs(cov_tmp-missing_val) < tol) {
                f_read = 0;
            }
            cov[n_cov] = cov_tmp;
        }

		/* add current individual to the count */
		i_total++;
		n_total = family[fam].n_total;



		if(f_read != 0) {

			i_pheno++; /* add current phenotyped individual to the count */

			/* link indexes of all people and phenotype people in a pedigree */
			family[fam].pheno_typed[i_pheno]=ind;
			family[fam].pheno_typed_inv[ind]=i_pheno;

			/* save the trait value */
			family[fam].trait[i_pheno]=trait;

			/* save the covariate value */
			for(i_cov=1; i_cov<n_cov; i_cov++) {
				family[fam].cov[i_pheno][i_cov] = cov[i_cov];
			}
			family[fam].cov[i_pheno][n_cov] = cov[n_cov];
		} else {
			family[fam].pheno_typed_inv[ind]=0; /*pheno_typed_inv[ind] = 0 if ind hasn't pheno/cov */
		}

		/* if everyone in the family has been read in then initilize the count variables */
		if(i_total == n_total) {
			i_total = 0;
			i_pheno = 0;
		}
	}

	fclose(pedfile);

	free(fam_id);
	free(ind_id);
	free(fa);
	free(mo);
	free_dvector(cov,1,n_cov);
}

/* allocate memory for the vectors of the sturcture FAMILY that link rank IDs (index) of everyone and those
 (not)phenotyped to those (not)genotyped, these vectors are reused for each marker and the memory
 allocated is more than needed */
void marker_skrats_allocate(struct FAMILY *family, struct DATA_STRUCT data_struct) {

	int n_fam = data_struct.n_fam;

	int i;
	for (i=1; i<=n_fam; i++) {
		family[i].geno_pheno_typed = imatrix(1,family[i].n_pheno,1,2);
		family[i].notgeno_pheno_typed = imatrix(1,family[i].n_pheno,1,2);
		if (family[i].n_total>family[i].n_pheno) {
			family[i].geno_notpheno_typed = ivector(1,family[i].n_total-family[i].n_pheno);
		}
		family[i].y = dvector(1,family[i].n_total);
	}
}

/* initialize vectors/matrices that link rank IDs */
void marker_skrats_initialize(struct FAMILY *family, struct DATA_STRUCT data_struct) {

	int n_fam = data_struct.n_fam;
	int n_total, n_pheno, n_notpheno;

	int i,j;
	for (i=1; i<=n_fam; i++) {
		n_total = family[i].n_total;
		n_pheno = family[i].n_pheno;
		n_notpheno = n_total-n_pheno;
		for(j=1; j<= n_pheno; j++) {
			family[i].geno_pheno_typed[j][1] = 0;
			family[i].notgeno_pheno_typed[j][1] = 0;
			family[i].geno_pheno_typed[j][2] = 0;
			family[i].notgeno_pheno_typed[j][2] = 0;
		}
		for(j=1; j<= n_notpheno; j++) {
			family[i].geno_notpheno_typed[j] = 0;
		}
		for(j=1; j<= n_total; j++) {
			family[i].y[j] = 9.0;
		}
	}

}

/* free the memory for the vectors/matrices that rank IDs*/
void marker_skrats_free(struct FAMILY *family, struct DATA_STRUCT data_struct) {

	int n_fam = data_struct.n_fam;

	int i;
	for (i=1; i<=n_fam; i++) {
		free_imatrix(family[i].geno_pheno_typed,1,family[i].n_pheno,1,2);
		free_imatrix(family[i].notgeno_pheno_typed,1,family[i].n_pheno,1,2);
		if (family[i].n_total>family[i].n_pheno) {
			free_ivector(family[i].geno_notpheno_typed,1,family[i].n_total-family[i].n_pheno);
		}
		free_dvector(family[i].y,1,family[i].n_total);
	}
}

/* function to read in the header line of the marker genotype file, it assumed that the file has
 already been opened and that the total number of genotyped individuals is known */

//void readmarker_header(FILE *markfile, char **id_geno2all, struct FAMILY *family, struct HASH hash, struct DATA_STRUCT data_struct) {
void readmarker_header(gzFile markfile, char **id_geno2all, struct FAMILY *family, struct HASH hash, struct DATA_STRUCT data_struct) {

	int marker=0, i, j, k, ind, fam;
	int n_fam = data_struct.n_fam;
	int n_typed = data_struct.n_typed;
	int n_total = data_struct.n_total;

	char *ind_id;
	ind_id = (char *) malloc(MAXLEN);

	int max_char = 10*n_total;
	char *line = (char*) malloc(max_char);
	//fgets(line,max_char,markfile);
	gzgets(markfile,line,max_char);

	int ret;

	int n_read=0; /* to keep track of how much of the line has been read with sscanf */
	int n_read_tmp=0; /* to store how much was read in the last call to sscanf */

	int i_geno_pheno=0, i_notgeno_pheno=0, i_geno_notpheno=0;

	/* initialize counts for each family */
	for(i=1; i<= n_fam; i++) {
		family[i].n_geno = 0;
		family[i].n_geno_pheno = 0;
		family[i].n_notgeno_pheno = 0;
		family[i].n_geno_notpheno = 0;
	}

	/* === scan the line === */

    char *chr, *markerid, *cm, *bp;
    chr = (char*) malloc(4);
    markerid = (char*) malloc(MAXLEN);
    cm = (char*) malloc(MAXLEN);
    bp = (char*) malloc(MAXLEN);
    sscanf(line, "%s %s %s %s %n", chr, markerid, cm, bp, &n_read_tmp);
    free(chr);
    free(markerid);
    free(cm);
    free(bp);


	//sscanf(line, "%d %n", &marker, &n_read_tmp); /* read in the first column, the marker number */
	n_read += n_read_tmp; /* update how much of the string has been read by sscanf */

	khash_t(str) *ind2typed = kh_init(str);

	for(i=1;i<=n_typed;i++) {
		sscanf(line+n_read, "%s %n", ind_id, &n_read_tmp);
		n_read += n_read_tmp;

		k = kh_get(str,hash.ind2ind,ind_id);
		ind = kh_value(hash.ind2ind,k);

		k = kh_get(str,hash.ind2fam,ind_id);
		fam = kh_value(hash.ind2fam,k);

		family[fam].n_geno++;
		if(family[fam].pheno_typed_inv[ind] != 0) {
			family[fam].n_geno_pheno++;
			i_geno_pheno = family[fam].n_geno_pheno;
			family[fam].geno_pheno_typed[i_geno_pheno][1]=ind;
			family[fam].geno_pheno_typed[i_geno_pheno][2]=family[fam].pheno_typed_inv[ind];
		} else {
			family[fam].n_geno_notpheno++;
			i_geno_notpheno = family[fam].n_geno_notpheno;
			family[fam].geno_notpheno_typed[i_geno_notpheno]=ind;
		}

		id_geno2all[i] = strdup(ind_id);

		k = kh_put(str,ind2typed,strdup(ind_id),&ret);
		kh_value(ind2typed,k) = 1;
	}

	/* find everyone who is not genotyped, loop through everyone in the data and check
	 if they are among those genotyped */
	const char *tmp_ind_id;
	for (k = kh_begin(hash.ind2ind); k != kh_end(hash.ind2ind); k++) {
		if (kh_exist(hash.ind2ind,k)) {
			tmp_ind_id = kh_key(hash.ind2ind,k);
			if (kh_get(str,ind2typed,tmp_ind_id) == kh_end(ind2typed)) {
				j = kh_get(str,hash.ind2ind,tmp_ind_id);
				ind = kh_value(hash.ind2ind,j);
				j = kh_get(str,hash.ind2fam,tmp_ind_id);
				fam = kh_value(hash.ind2fam,j);

				if(family[fam].pheno_typed_inv[ind] != 0) {
					family[fam].n_notgeno_pheno++;
					i_notgeno_pheno = family[fam].n_notgeno_pheno;
					family[fam].notgeno_pheno_typed[i_notgeno_pheno][1]=ind;
					family[fam].notgeno_pheno_typed[i_notgeno_pheno][2]=family[fam].pheno_typed_inv[ind];
				}

			}
		}
	}

    /* find who is genotyped (=1) and who is not genotyped (=0) */
	for (k = kh_begin(hash.ind2ind); k != kh_end(hash.ind2ind); k++) {
		if (kh_exist(hash.ind2ind,k)) {
			tmp_ind_id = kh_key(hash.ind2ind,k);
			j = kh_get(str,hash.ind2ind,tmp_ind_id);
			ind = kh_value(hash.ind2ind,j);
			j = kh_get(str,hash.ind2fam,tmp_ind_id);
			fam = kh_value(hash.ind2fam,j);
			if (kh_get(str,ind2typed,tmp_ind_id) == kh_end(ind2typed)) {
				family[fam].geno_typed[ind] = 0;
			} else {
				family[fam].geno_typed[ind] = 1;
			}
		}
	}

	free(ind_id);
	free(line);
}

/* function to fill in the info on which subset phenotyped individuals belongs to
 called after readmarker_header and not after readmarker, for the purpose of finding
 the subset of phenotyped individuals to use in VC estimation */

void fill_pheno_subset(struct FAMILY *family, struct DATA_STRUCT data_struct) {
 	int n_fam, n_total, n_pheno, i_fam, j_total, k_total, i_pheno, j_pheno;
    int geno=0, genorel=0, phenorel_w_genorel;
    double phi_ij=0, phi_jk=0;

    int *geno_typed;
    int *pheno_typed;
    double **phi;

	n_fam = data_struct.n_fam;

    for(i_fam = 1; i_fam <= n_fam; i_fam++) {
        n_total = family[i_fam].n_total;
        n_pheno = family[i_fam].n_pheno;
        geno_typed = family[i_fam].geno_typed;
        pheno_typed = family[i_fam].pheno_typed;
        phi = family[i_fam].phi;

        if(n_pheno > 0) {
            for(i_pheno = 1; i_pheno<=n_pheno; i_pheno++) {

                family[i_fam].pheno_subset[i_pheno] = 0; /* initialize */

                geno = 0;
                geno = geno_typed[pheno_typed[i_pheno]];

               if(geno == 1) {
                    /*check if phenotyped person is genotyped */
                   family[i_fam].pheno_subset[i_pheno] = 1;

                } else {
                    /* else if not genotyped check if any relatives genotyped */
                    j_total = 1;
                    genorel = 0;
                    while(j_total<=n_total && genorel==0) {
                        phi_ij = phi[pheno_typed[i_pheno]][j_total];
                        if(phi_ij>0 && geno_typed[j_total]==1) {
                            family[i_fam].pheno_subset[i_pheno] = 2;
                            genorel = 1;
                        }
                        j_total++;
                    }
                    if(genorel == 0) {
                        /* if gone through all relatives and no one genotyped then
                         check if a phenotyped relative has a genotyped relative */
                        for(j_pheno = i_pheno+1; j_pheno <= n_pheno; j_pheno++){
                            phi_ij = phi[pheno_typed[i_pheno]][pheno_typed[j_pheno]];
                            k_total = 1;
                            phenorel_w_genorel = 0;
                            while(phi_ij>0 && k_total<=n_total && phenorel_w_genorel==0) {
                                phi_jk = phi[pheno_typed[j_pheno]][k_total];
                                if(phi_jk>0 && geno_typed[k_total]==1) {
                                    family[i_fam].pheno_subset[i_pheno] = 3;
                                    phenorel_w_genorel = 1;
                                }
                                k_total++;
                            }
                        }
                    } /*end: if genorel == 0 */
                } /* end: if person genotyped else not genotyped*/
            } /*end: for loop through all phenotyped */
        } /* end: if n_pheno > 1 */
    } /* end: for(i_fam = 1; i_fam <= n_fam; i_fam++) */
 } /* end function: void fill_pheno_subset(struct FAMILY *family, struct DATA_STRUCT data_struct)*/

/* function to read in genotypes for one marker at a time, it is assumed that the file has
 already been opened for the first time with the function readmarker_header and the first
 line read in, so that the first time readmarker is called the second line of the file is read
 in (the genotypes for the first marker are stored in the second line) */
//void readmarker(FILE *markfile, struct FAMILY *family, struct HASH hash, char **id_geno2all, struct DATA_STRUCT data_struct, struct ANNOTATION *anno) {
void readmarker(gzFile markfile, struct FAMILY *family, struct HASH hash, char **id_geno2all, struct DATA_STRUCT data_struct, struct ANNOTATION *anno) {

	int n_fam = data_struct.n_fam;
	int n_typed = data_struct.n_typed;
	int n_total = data_struct.n_total;

	/* allocate enough memory so that the whole line is read in with the fgets call */
	int max_char = 4*n_total;
	char line[max_char];
	//fgets(line,max_char,markfile);
	gzgets(markfile,line,max_char);


	int marker=0;
	int n_read=0; /* to keep track of how much of the line has been read with sscanf */
	int n_read_tmp=0; /* to store how much was read in the last call to sscanf */

	char *ind_id;
	int ind, fam;

	int ret;
	int i_geno_pheno, i_notgeno_pheno, i_geno_notpheno;


	int i, j, k;
	int *geno;
	geno = ivector(1,2);

	khash_t(str) *ind2typed = kh_init(str);

	/* initialize counts for each family */
	for(i=1; i<= n_fam; i++) {
		family[i].n_geno = 0;
		family[i].n_geno_pheno = 0;
		family[i].n_notgeno_pheno = 0;
		family[i].n_geno_notpheno = 0;
	}

	/* === scan the line to count n_geno_pheno in each family === */

	int *fam_n_geno_pheno = ivector(1,n_fam);
	for(i=1; i<=n_fam; i++) {
		fam_n_geno_pheno[i] = 0;
	}

    char *fam_id, *ind1_id, *ind2_id;
	fam_id = (char*) malloc(MAXLEN);

    char *chr, *markerid;
    int cm, bp;
    chr = (char*) malloc(4);
    markerid = (char*) malloc(MAXLEN);
    sscanf(line, "%s %s %d %d %n", chr, markerid, &cm, &bp, &n_read_tmp);
    anno->chr = chr;
    anno->markername = markerid;
    anno->centimorgan = cm;
    anno->basepair =  bp;

	//sscanf(line, "%d %n", &marker, &n_read_tmp); /* read in the first column, the marker number */
	n_read += n_read_tmp; /* update how much of the string has been read by sscanf */

	for(i=1;i<=n_typed;i++) {

		ind_id = id_geno2all[i];

		k = kh_get(str,hash.ind2ind,ind_id);
		ind = kh_value(hash.ind2ind,k);

		k = kh_get(str,hash.ind2fam,ind_id);
		fam = kh_value(hash.ind2fam,k);

		sscanf(line+n_read, "%1d %1d %n", &geno[1], &geno[2], &n_read_tmp);
		n_read += n_read_tmp;



		if(geno[1]!=0 && geno[2]!=0) {
			family[fam].n_geno++;
			if(family[fam].pheno_typed_inv[ind]!=0) {
				fam_n_geno_pheno[fam]++;
			}
		}
	}

	/* === scan the line again to create links and fill in genotypes === */

	n_read = 0;
	n_read_tmp = 0;
    sscanf(line, "%s %s %d %d %n", chr, markerid, &cm, &bp, &n_read_tmp);
	//sscanf(line, "%d %n", &marker, &n_read_tmp); /* read in the first column, the marker number */
	n_read += n_read_tmp; /* update how much of the string has been read by sscanf */

	for(i=1;i<=n_typed;i++) {

		ind_id = id_geno2all[i];

		k = kh_get(str,hash.ind2ind,ind_id);
		ind = kh_value(hash.ind2ind,k);

		k = kh_get(str,hash.ind2fam,ind_id);
		fam = kh_value(hash.ind2fam,k);

		/* read in each marker and save link between indexes, and additive genotype code (calls
		 the function addgeno) */
		sscanf(line+n_read, "%1d %1d %n", &geno[1], &geno[2], &n_read_tmp);
		n_read += n_read_tmp;

		if(geno[1]!=0 && geno[2]!=0) {
			if(family[fam].pheno_typed_inv[ind]!=0) {
				family[fam].n_geno_pheno++;
				i_geno_pheno = family[fam].n_geno_pheno;
				family[fam].geno_pheno_typed[i_geno_pheno][1] = ind;
				family[fam].geno_pheno_typed[i_geno_pheno][2] = family[fam].pheno_typed_inv[ind];
				addgeno(geno,1,&(family[fam].y[i_geno_pheno]));
			} else {
				family[fam].n_geno_notpheno++;
				i_geno_notpheno = family[fam].n_geno_notpheno;
				family[fam].geno_notpheno_typed[i_geno_notpheno] = ind;
				addgeno(geno,1,&(family[fam].y[i_geno_notpheno+fam_n_geno_pheno[fam]]));
			}

			k = kh_put(str,ind2typed,strdup(ind_id),&ret);
			kh_value(ind2typed,k) = 1;
		}
	}

	/* find everyone who is not genotyped, loop through everyone in the data and check
	   if they are among those genotyped */
	const char *tmp_ind_id;
	for (k = kh_begin(hash.ind2ind); k != kh_end(hash.ind2ind); k++) {
		if (kh_exist(hash.ind2ind,k)) {
			tmp_ind_id = kh_key(hash.ind2ind,k);
			if (kh_get(str,ind2typed,tmp_ind_id) == kh_end(ind2typed)) {
				j = kh_get(str,hash.ind2ind,tmp_ind_id);
				ind = kh_value(hash.ind2ind,j);
				j = kh_get(str,hash.ind2fam,tmp_ind_id);
				fam = kh_value(hash.ind2fam,j);

				if(family[fam].pheno_typed_inv[ind] != 0) {
					family[fam].n_notgeno_pheno++;
					i_notgeno_pheno = family[fam].n_notgeno_pheno;
					family[fam].notgeno_pheno_typed[i_notgeno_pheno][1]=ind;
					family[fam].notgeno_pheno_typed[i_notgeno_pheno][2]=family[fam].pheno_typed_inv[ind];
				}

			}
		}
	}

	/* free all keys in ind2typed and the hashtable itself */
	for (k = kh_begin(ind2typed); k != kh_end(ind2typed); k++){
		if (kh_exist(ind2typed,k)) {
			free((char*) kh_key(ind2typed,k)); /* cast away constness */
		}
	}
	kh_destroy(str,ind2typed);


	free_ivector(geno,1,2);
	free_ivector(fam_n_geno_pheno,1,n_fam);
}


/* function to recode the genotype into an addtive model, allele = non-risk allele
 assumes markers have only 2 alleles (e.g. SNPs) */
void addgeno(int *geno, int allele, double *y) {

	if(geno[1] == allele && geno[2] == allele) {
		*y = 0.0;
	}

	if(geno[1] == allele && geno[2] != allele) {
		*y = 0.5;
	}

	if(geno[2] == allele && geno[1] != allele) {
		*y = 0.5;
	}

	if(geno[1] != allele && geno[2] != allele) {
		*y = 1.0;
	}

	if(geno[1] == 0 && geno[2] == 0) {
		*y = 9.0;
	}


}

/* function to read int the kinship coefficient matrix and save in the phi matrix
 phi_ii = 1+coef_ii and phi_ij = 2*coef_ij */
void readkin(char *filename, struct FAMILY *family, struct HASH hash) {

	char *fam_id, *ind1_id, *ind2_id;
	fam_id = (char*) malloc(MAXLEN);
	ind1_id = (char*) malloc(MAXLEN);
	ind2_id = (char*) malloc(MAXLEN);

	double coef = 0;
	int fam = 0, ind1 = 0, ind2 = 0, k;


	FILE *kinfile;
	kinfile = fopen(filename, "r");

	if(kinfile == NULL) {
		printf("Can't open file kinship coefficent file\n");
		exit(1);
	}

	while (!feof(kinfile)) {
		if (fscanf(kinfile,"%s %s %s %lf ", fam_id, ind1_id, ind2_id, &coef) != 4) {
			break;
		}

		k = kh_get(str,hash.fam2fam,fam_id);
		fam = kh_value(hash.fam2fam,k);

		k = kh_get(str,hash.ind2ind,ind1_id);
		ind1 = kh_value(hash.ind2ind,k);

		k = kh_get(str,hash.ind2ind,ind2_id);
		ind2 = kh_value(hash.ind2ind,k);

		if(ind1==ind2) {
			coef = coef + 1.0;
			family[fam].phi[ind1][ind2] = coef;
		}
		else {
			coef = 2.0*coef;
			family[fam].phi[ind1][ind2] = coef;
			family[fam].phi[ind2][ind1] = coef;
		}
	}

	fclose(kinfile);
	free(fam_id);
	free(ind1_id);
	free(ind2_id);

}
