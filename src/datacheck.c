#include "datacheck.h"
#include <string.h>

//KHASH_MAP_INIT_STR(str,cons char*);

/* function to count a minimum (relpair) number of relative pairs who are
   phenotyped (subset=0),
   both phenotyped and genotyped (subset=1)
   phenotyped and 1) genotyped or 2) with genotyped relative (subset=2)
   phenotyped and 1) genotyed or 2) with genotyped relative or 3) with phenotyped relative who has genotyped relative (subset=3)
   Returns integer = 0 if enough pairs are found o.w. 1 */
int unrelated_check (struct FAMILY *family, struct DATA_STRUCT data_struct, int subset, int relpair) {
	int n_fam, n_pheno, n_geno_pheno, i, j, k, unr=1, rel = 0;
	n_fam = data_struct.n_fam;
	double phi;
    int *pheno_subset;

    /* all pheno */
	if(subset == 0) {
		for(i = 1; i<=n_fam; i++) {
			n_pheno = family[i].n_pheno;
			if(n_pheno > 1) {
				for(j = 1; j<n_pheno; j++) {
					for(k=j+1; k<=n_pheno; k++) {
						phi = family[i].phi[family[i].pheno_typed[j]][family[i].pheno_typed[k]];
						if(phi > 0) {rel = rel + 1;}
					} /* end for k */
				} /* end for j*/
			} /* end if n_pheno > 1*/
		} /* end for i*/

		if(rel >= relpair) {
			unr = 0;
		}
	}

    /* all pheno_geno */
    if(subset == 1) {
		for(i = 1; i<=n_fam; i++) {
			n_geno_pheno = family[i].n_geno_pheno;
			if(n_geno_pheno > 1) {
				for(j = 1; j<n_geno_pheno; j++) {
					for(k=j+1; k<=n_geno_pheno; k++) {
						phi = family[i].phi[family[i].geno_pheno_typed[j][1]][family[i].geno_pheno_typed[k][1]];
						if(phi > 0) {rel = rel + 1;}
					}
				}
			}
		}
		if(rel >= relpair) {
			unr = 0;
		}
	}

    /* all pheno_geno or pheno_genorel */
	if(subset == 2) {
		for(i=1; i<=n_fam; i++) {
            pheno_subset = family[i].pheno_subset;
			n_pheno = family[i].n_pheno;
			if(n_pheno > 1) {
				for(j = 1; j<n_pheno;j++) {
					for(k=j+1;k<=n_pheno;k++){
                        if((pheno_subset[j] == 1 || pheno_subset[j] == 2) && (pheno_subset[k] == 1 || pheno_subset[k] == 2)) {
                            phi = family[i].phi[family[i].pheno_typed[j]][family[i].pheno_typed[k]];
                            if(phi >0) { rel = rel+1;}
                        }
					}
				}
			}
		}
		if(rel >= relpair) {
			unr = 0;
		}
	}

    /* all pheno_geno, pheno_genorel, pheno_phenorel_w_genorel*/
	if(subset == 3) {
		for(i=1; i<=n_fam; i++) {
            pheno_subset = family[i].pheno_subset;
			n_pheno = family[i].n_pheno;
			if(n_pheno > 1) {
				for(j = 1; j<n_pheno;j++) {
					for(k=j+1;k<=n_pheno;k++){
                        if((pheno_subset[j] == 1 || pheno_subset[j] == 2 || pheno_subset[j] == 3) && (pheno_subset[k] == 1 || pheno_subset[k] == 2 || pheno_subset[k] == 3)) {
                            phi = family[i].phi[family[i].pheno_typed[j]][family[i].pheno_typed[k]];
                            if(phi >0) { rel = rel+1;}
                        }
					}
				}
			}
		}
		if(rel >= relpair) {
			unr = 0;
		}
	}


	return unr;
}

/* ================================================================================================= */
/* ================================================================================================= */
/* ================================================================================================= */

void format_checks(char *pedfilename, char *markerfilename, char *kinfilename, struct DATA_STRUCT *data_struct) {

    /* define, initialize, and allocate memory for variables */
    int n_cov, total, typed, n_fam;
    int n_mark=-1; /* -1 to subtract the header line */
    int marker;
    int i;
    int i_cov;

    int ret,k;

    char *fam_id0,*fam_id, *ind_id, *fa, *mo, *ind_id1, *ind_id2;
    fam_id0 = (char*) malloc(MAXLEN);
    fam_id = (char*) malloc(MAXLEN);
    ind_id = (char*) malloc(MAXLEN);
    fa = (char*) malloc(MAXLEN);
    mo = (char*) malloc(MAXLEN);
    ind_id1 = (char*) malloc(MAXLEN);
    ind_id2 = (char*) malloc(MAXLEN);

    double kincoef;

    int sex;
    double missing_val = data_struct->missing_val;
    double trait = missing_val;
    double cov_tmp = missing_val;
    int n_read=0; /* to keep track of how much of the line has been read with sscanf */
    int n_read_tmp=0; /* to store how much was read in the last call to sscanf */

    int col=0;
    char car;

    
    /* create internal hash tables */

    khash_t(str) *indlist = kh_init(str);
    khash_t(str) *famlist = kh_init(str);
    khash_t(str) *falist = kh_init(str);
    khash_t(str) *molist = kh_init(str);
    khash_t(str) *ind2kin = kh_init(str);

    /* ================================================================================================= */
    /* check if files are present and can be opened input files */
	FILE *pedfile, *kinfile;
	gzFile markerfile;

    pedfile = fopen(pedfilename, "r");
    markerfile = gzopen( markerfilename, "r");
    kinfile = fopen(kinfilename, "r");

    if(pedfile == NULL) {
        printf("ERROR: Can't open pedigree and trait file\n");
        exit(1);
    }
    if(markerfile == NULL) {
        printf("ERROR: Can't open genotype data file\n");
        exit(1);
    }
    if(kinfile == NULL) {
        printf("ERROR: Can't open kinship information file\n");
        exit(1);
    }

    fclose(pedfile);
	gzclose(markerfile);
    fclose(kinfile);
    
    /* ================================================================================================= */


    /* ================================================================================================= */
    /* read through the pedfile to count the number of columns, to then find the number of covariates */
    pedfile = fopen(pedfilename, "r");

    car = 0;
    while(car != EOF && car != '\n') {
        car = getc(pedfile);
        if(car == ' ' || car == '\t') {
            col++;
        }
    }
    col++;
    n_cov = col-5; /* no. of covariates + 1 (+1 for the intercept which is assumed to be absent in pedfile) */
    fclose(pedfile);
    
   
    
    
    /* ================================================================================================= */


    /* ================================================================================================= */
    /* read through the pedfile to create lists of pedigree and individual ids,
     count number of families (n_fam) and number of individuals (total) */
    //fam_id0 = strdup("0");
    strcpy(fam_id0,"0");
	total = 0;
    n_fam = 0;
    pedfile = fopen(pedfilename, "r");

    
    while(!feof(pedfile)) {

        /* ------------------------------------------------ */

        if(n_cov>1){
            fscanf(pedfile,"%s %s %s %s %d %lf ", fam_id, ind_id, fa, mo, &sex, &trait);
        }
        if (n_cov==1) {
            fscanf(pedfile,"%s %s %s %s %d %lf \n", fam_id, ind_id, fa, mo, &sex, &trait);
        }

        if(n_cov>1) {
            for(i_cov=2; i_cov<n_cov; i_cov++) {
                fscanf(pedfile,"%lf ", &cov_tmp);
            }
            fscanf(pedfile, "%lf \n", &cov_tmp);
        }
        total++;
        /* ------------------------------------------------ */

        /* check that the current family ID is also the last family ID, if it is not the last
         family ID check if it appears before that, indicating that the family IDs appear in
         two distinct groups, then exit o.w. update n_fam and "last family" aka fam_id0 */
         if(strcmp(fam_id,fam_id0)!=0) {
            ped_famorder(famlist, fam_id);
            n_fam++;
        }
        
        
        strcpy(fam_id0,fam_id);

        /* ------------------------------------------------ */

        /* check that the current individuals ID has not been added to the hash table, if it has
         that indicates that the individual appears twice in the pedigree file */
        ped_uniqueind(indlist, ind_id);

        /* ------------------------------------------------ */

        /* create hash tables */
        k = kh_put(str,famlist,strdup(fam_id), &ret);
        kh_value(famlist,k)=1;

        k = kh_put(str,indlist,strdup(ind_id), &ret);
        kh_value(indlist,k)=1;

        k = kh_put(str,falist,strdup(fa), &ret);
        kh_value(falist,k)=1;

        k = kh_put(str,molist,strdup(mo), &ret);
        kh_value(molist,k)=1;

    } /* closes while(!feof(pedfile)), which finishes counting total and n_fam */
    
    
    fclose(pedfile);
    /* ================================================================================================= */


    /* ================================================================================================= */
    /* check if all father ids and mother ids present as ids in the pedfile
       does not check whether the father and mother have the same family id as the child */
    const char *tmp_ind_id;
    for(k=kh_begin(falist); k!=kh_end(falist);k++) {
        if(kh_exist(falist,k)) {
            tmp_ind_id = kh_key(falist,k);
            if(strcmp("0",tmp_ind_id)!=0) {
                if(kh_get(str,indlist, tmp_ind_id) == kh_end(indlist)) {
                    printf("ERROR: Father id %s is not present as id in the pedigree file\n",tmp_ind_id);
                    exit(1);
                }
            }
        }
    }
    for(k=kh_begin(molist); k!=kh_end(molist);k++) {
        if(kh_exist(molist,k)) {
            tmp_ind_id = kh_key(molist,k);
            if(strcmp("0",tmp_ind_id)!=0) {
                if(kh_get(str,indlist, tmp_ind_id) == kh_end(indlist)) {
                    printf("ERROR: Mother id %s is not present as id in the pedigree file\n",tmp_ind_id);
                    exit(1);
                }
            }
        }
    }
    /* ================================================================================================= */


    /* ================================================================================================= */
    /* read through the markerfile to count the number of markers */
	markerfile = gzopen(markerfilename, "r");

    car = 0;
    while (car != EOF) {
		car = gzgetc(markerfile);
        if(car == '\n') {
            n_mark++;
        }
    }
	gzclose(markerfile);

    
    /* read through the first line of the markerfile to count the number of columns, to then find the number genotyped individuals */
	markerfile = gzopen(markerfilename, "r");

    col = 1;
    car = 0;
    while(car != EOF && car != '\n') {
		car = gzgetc(markerfile);
        if(car == ' ' || car == '\t') {
            col++;
        }
    } /* closes while(car != EOF && car != '\n') */
    typed = col-4; /*no of genotyped inviduals, subtract the 4 annotation columns*/
	gzclose(markerfile);
    /* ================================================================================================= */


    /* ================================================================================================= */
    /* scan the line and create a list of all ids of genotyped individuals */
    int max_char = 10*total;
    char *line = (char*) malloc(max_char);
	markerfile = gzopen(markerfilename, "r");
	gzgets(markerfile,line,max_char);
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
    n_read += n_read_tmp; /* update how much of the string has been read by sscanf */

    khash_t(str) *ind2typed = kh_init(str);


    for(i=1;i<=typed;i++) {
        sscanf(line+n_read, "%s %n", ind_id, &n_read_tmp);
        n_read += n_read_tmp;

        /* check that the current individuals ID has not been added to the hash table, if it has
         that indicates that the individual appears twice in the marker file */
        marker_uniqueind(ind2typed, ind_id);

        k = kh_put(str,ind2typed,strdup(ind_id),&ret);
        kh_value(ind2typed,k) = 1;
    }

	gzclose(markerfile);
    
    
    
    /* ================================================================================================= */

     /* check if everyone in the markerfile is in the pedfile */
    for(k=kh_begin(ind2typed); k!=kh_end(ind2typed);k++) {
        if(kh_exist(ind2typed,k)) {
            tmp_ind_id = kh_key(ind2typed,k);
            if(kh_get(str,indlist, tmp_ind_id) == kh_end(indlist)) {
                printf("ERROR: Individual %s in the genotype data file is not in the pedigree file\n",tmp_ind_id);
                exit(1);
            }
        }
    }
    
    
    /* ================================================================================================= */
    /* read through kinship file to check if all individuals are in the pedigree file */

    kinfile = fopen(kinfilename, "r");
    
    while(!feof(kinfile)) {

        fscanf(kinfile,"%s %s %s %lf ", fam_id, ind_id1, ind_id2, &kincoef);
        
        k = kh_put(str,ind2kin,strdup(ind_id1),&ret);
        kh_value(ind2kin,k)=1;
        k = kh_put(str,ind2kin,strdup(ind_id2),&ret);
        kh_value(ind2kin,k)=1;
    }
    

    fclose(kinfile);

    for(k=kh_begin(ind2kin); k!=kh_end(ind2kin);k++) {
        if(kh_exist(ind2kin,k)) {
            tmp_ind_id = kh_key(ind2kin,k);
            if(kh_get(str,indlist, tmp_ind_id) == kh_end(indlist)) {
                printf("ERROR: Individual %s in the kinship file is not in the pedigree file\n",tmp_ind_id);
                exit(1);
            }
        }
    }
    
    

    /* ================================================================================================= */

     /* fill data_struct */
    data_struct->n_fam = n_fam;
    data_struct->n_marker = n_mark;
    data_struct->n_total = total;
    data_struct->n_typed = typed;
    data_struct->n_cov = n_cov;

    /* ================================================================================================= */

    /* free memory and destroy internal hash tables */

    free(ind_id);
    free(ind_id1);
    free(ind_id2);
    free(fam_id0);
    free(fam_id);
    free(fa);
    free(mo);
    free(line);

    kh_destroy(str,indlist);
    kh_destroy(str,famlist);
    kh_destroy(str,falist);
    kh_destroy(str,molist);
    kh_destroy(str,ind2kin);
    
    
    

} /* close format_checks */

/* ================================================================================================= */
/* ================================================================================================= */
/* ================================================================================================= */

int ped_famorder(khash_t(str) *fam_list, const char *fam) {
    int k;
    if(kh_get(str,fam_list,fam) == kh_end(fam_list)) {
        return 1;
    } else {
        printf("ERROR: Family IDs are not grouped together in the pedigree file. ID %s appears in at least two distinct groups\n",fam);
        exit(1);
    }
}

/* ================================================================================================= */
/* ================================================================================================= */
/* ================================================================================================= */

int ped_uniqueind(khash_t(str) *ind_list, const char *ind ) {
    int k;
    if(kh_get(str,ind_list,ind) == kh_end(ind_list)) {
        return 1;
    } else {
        printf("ERROR: Individual %s appears at least twice in the pedigree file \n",ind);
        exit(1);
    }
}

/* ================================================================================================= */
/* ================================================================================================= */
/* ================================================================================================= */

int marker_uniqueind(khash_t(str) *ind_to_typed, const char *ind){
    int k;
    if(kh_get(str,ind_to_typed,ind) == kh_end(ind_to_typed)) {
        return 1;
    } else {
        printf("ERROR: Individual %s appears at least twice in the marker file \n", ind);
        exit(1);
    }
}

/* ================================================================================================= */
/* ================================================================================================= */
/* ================================================================================================= */
