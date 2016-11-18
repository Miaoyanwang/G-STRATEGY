/* Miaoyan Wang, Aug 23, 2016 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "nrutil.c"
#include <unistd.h> /* for fork */
#include <time.h>
#include "clapack_routine/dposv.h"

#define MAXLEN 2048
#define TINY 1.0e-20
#define STOPTEM 1.0e-5

struct FAM {
    int famID;
    int *indID;
    int r;
    int *R;
    double *D;
    int *fix;
    int *G;
    int *subR; 
    double *subD;
    int famn;
    int c_famn;
    int famn_best;
    int c_famn_best;
    int *c_subR;
    double *c_subD;
    int *subR_best;
    double *subD_best;
    int *c_subR_best;
    double *c_subD_best;
    int nfix;
    int ngeno;
    double **phiN;
    double **invN;
    double *E;
    double *famcomp;
    double **matrix;
    //double **Cov;
    int *enrich_select;
}FAM;


struct sortarray{
    double value;
    int ind;
    int fam;
    int select;
};

int comparray(const void *i, const void *j);
void Sortarray(struct sortarray *array, int n);



char phenofile[MAXLEN]="phenotype.txt", sizefile[MAXLEN]="size.txt",kinship[MAXLEN]="kinship.txt",outfile[MAXLEN]="g_strategy.txt";
FILE *fpheno,*fsize,*fkinship,*fpheno_MASTOR,*fgeno_MASTOR, *ftpr,*fenrich,*seedfl,*fp_R,*fp_p,*fp_m,*fp_hp;
int n,nindiv, totalgeno, totalfix, f, fam,i,j,k,t;

void writevector(FILE *fp, struct FAM *Info, int f);
void readsize(FILE *fsize,int *nindiv, int *n, int *f, int *ncov, double *prevalance,int trait_type);
int checkformat(FILE *fpheno,  int ncov, int nindiv, int n, int f, int *familysize);
void readin(FILE *fpheno,FILE *fkinship, struct FAM *Info, int f,int *totalgeno, int *n,  int trait_type, double prevalance, int ncov, int nindiv,FILE *fpheno_MASTOR,FILE *fgeno_MASTOR,int *famIDvector);
int findInd(int n_indiv, int *indivlist, int indiv);
void MatchTPR(FILE *fpheno,FILE *ftpr, int *famIDvector, struct FAM *Info, int f);
double random2();
double tune(struct FAM *Info,int *tem_subR, int *tem_subR2,  int totalgeno, int f, int n,  int *temR, double *tem_m[], double *tem_x[], double *famcomp_sum, double *famcomp_sumnew, double **famcomp_pre, int nq, int tem1, double tem2, int ind_r, int ind_a, int fam_r, int fam_a,int totalfix);
void neigh(struct FAM *Info, int *ind_r, int *fam_r, int *ind_a, int *fam_a, int *tem1, double *tem2);
void unneigh(struct FAM *Info, double **famcomp_pre,int *fam_r,int *ind_r,int *fam_a,int *ind_a,int *tem1, double *tem2);
int findInd(int n_indiv, int *indivlist, int indiv);
void print_dvector(double *vector,int n);
int comp(const void *, const void *);
void update(struct FAM *Info, int f, int reverse);
void print_vector(struct FAM *Info, int n);
int metrop(double de, double tem);
void solve(double **y,double **a,int n,double *m,double *x);
void print_ivector(int *vector,int n, int sort);
void sample(int *tem_subR, int *tem_subR2, int r, int n,int *temR, int totalfix);
double tune(struct FAM *Info,int *tem_subR, int *tem_subR2,  int totalgeno, int f, int n, int *temR, double *tem_m[], double *tem_x[], double *famcomp_sum, double *famcomp_sumnew, double **famcomp_pre, int nq, int tem1, double tem2, int ind_r, int ind_a, int fam_r, int fam_a,int totalfix);
void neigh(struct FAM *Info, int *ind_r, int *fam_r, int *ind_a, int *fam_a, int *tem1, double *tem2);
void unneigh(struct FAM *Info, double **famcomp_pre,int *fam_r,int *ind_r,int *fam_a,int *ind_a,int *tem1, double *tem2);
void famsample(int *tem_subR,int *tem_subR2, int totalgeno,int n, struct FAM *Info);
double diffpower(struct FAM *Info, double *famcomp_sum,double *famcomp_sumnew,double **famcomp_pre,int *ind_r, int *fam_r, int *ind_a, int *fam_a, double *tem_m[],double *tem_x[]);
int uphill(struct FAM *Info, int f, double *powersubR, double *famcomp_sum, double *famcomp_sumnew, double **famcomp_pre,double *tem_m[], double *tem_x[], double *pbest, int *tem1, double *tem2);
int equal(int n);
void famcomponent(struct FAM *Info, double *tem_m[],double *tem_x[], int f,double *famcomp_sum);



int main(int argc, char *argv[]){
    
    printf("\n"
           "  |----------------------------------------------------------|\n"
           "  |                        G-STRATEGY:                       |\n"
           "  |      Optimal Selection of Individuals for Sequencing     |\n"
           "  |               in Genetic Association Studies             |\n"
           "  |                                                          |\n"
           "  |                Version 1.0 - Aug 22, 2016                |\n"
           "  |                                                          |\n"
           "  |                     Copyright(C) 2016                    |\n"
           "  | Miaoyan Wang, Johanna Jakobsdottir, and Mary Sara McPeek |\n"
           "  |                     License: GNU GPL v3                  |\n"
           "  |                                                          |\n"
           "  |  http://galton.uchicago.edu/∼mcpeek/software/GSTRATEGY   |\n"
           "  |----------------------------------------------------------|\n\n"
           );

    
    long seed;
    int trait_type=0;
    int pfile=0;
    int ncov=0;
    double prevalance;
    int *familysize;
    int strategy_option=1;
    double *famcomp_sum,*famcomp_sumnew,**famcomp_pre;
    int *temR,*tem_subR,*tem_subR2;
    clock_t begin,end;
    double time_spent;
    int nover,nlimit,para_opt=0,nsucc=1;
    int ans,q;
    int n,nindiv, totalgeno, totalfix, f, fam,i,j,k,t;
    int fam_r,ind_r,fam_a,ind_a,tem1;
    double tem2;
    double de,pbest,tem,powersubR_new, powersubR,alpha;
    

    if (argc>1){
        for (i=1;i<argc && argv[i][0]=='-';i++){
            switch(argv[i][1])
            {
                case 'b':
                    printf("user specified objective function: M_QLS \n");
                    trait_type=1;
                   break;
                case 'q':
                    printf("user specified objective function: MASTOR \n");
                    trait_type=2;
                    break;
                case 'p':
                    strncpy(phenofile,argv[++i],MAXLEN);
                    printf("user specified phenotype file: %s\n",phenofile);
                    break;
                case 'k':
                    strncpy(kinship,argv[++i],MAXLEN);
                    printf("user specified kinship file: %s\n",kinship);
                    break;
                case 's':
                    strncpy(sizefile,argv[++i],MAXLEN);
                    printf("user specified size file: %s\n",sizefile);
                    break;
                case 'e':
                    //strncpy(outfile2,argv[++i],MAXLEN);
                    printf("user requested strategy: extreme enrichment selection \n");
                    strategy_option=2;
                    break;
            }
        }
    }
    
    if(trait_type==0){
        printf("ERROR: you must specify either -b or -q flag! \n");
        exit(1);
    }
    
   /*phenotype file*/
    if((fpheno=fopen(phenofile,"r"))==NULL){
        printf("Cannot open phenotype file.\n");
        exit(1);
    }
   
    
    /* kinship coefficient file */
    if((fkinship=fopen(kinship,"r"))==NULL){
        printf("Cannot open kinship file.\n");
        exit(1);
    }
    
    
    /* size file */
    if((fsize=fopen(sizefile,"r"))==NULL){
        printf("Cannot open size file.\n");
        exit(1);
    }
    

    fpheno_MASTOR=fopen("phenotype_MASTOR.txt","w");
    fgeno_MASTOR=fopen("genotype_MASTOR.txt","w");
   
    
   readsize(fsize,&nindiv,&n,&f,&ncov,&prevalance,trait_type);
    
    
    
    familysize=ivector(1,f);
    for (fam=1;fam<=f;fam++){
        familysize[fam]=0;
    }
    
    int *famIDvector=ivector(1,f);
    
   totalfix=checkformat(fpheno,ncov,nindiv,n,f,familysize);
  
    
    fclose(fpheno);
    
     struct FAM Info[f+1];
    
    

    for (fam=1;fam<=f;fam++)
    {
        Info[fam].famID=0;
        Info[fam].indID=ivector(1,familysize[fam]);
        Info[fam].r=0; /* Number of phenotyped individuals in family "fam"*/
        Info[fam].R=ivector(1,familysize[fam]); /* Individual Id's */
        Info[fam].D=dvector(1,familysize[fam]); /* Transformed phenotypic residuals */
        Info[fam].fix=ivector(1,familysize[fam]); /* a 0-1 vector labelling whether the individual has been initially genotyped or not; 1 = yes, 0 = no */
        Info[fam].G=ivector(1,familysize[fam]);/* a 0-1 vector labeling whether the individual is available to be selected; 1 = yes, 0 = no */
        Info[fam].subR=ivector(1,n); /* Candidate subset of individuals selected from family "fam" */
        Info[fam].subD=dvector(1,n); /* Transformed phenotypic residuals for the individuals in subR */
        Info[fam].nfix=0; /* Total number of initially genotyped individuals in family "fam "*/
        Info[fam].ngeno=0; /* Number of individuals not available to be selected in family "fam" */
        Info[fam].famn=0; /* Number of individuals selected into candidate subset from family "fam" */
        Info[fam].c_subR=ivector(1,familysize[fam]); /* Complement set of subR w.r.t. R */
        Info[fam].c_subD=dvector(1,familysize[fam]); /* Complement set of subD w.r.t. D */
        Info[fam].c_famn=0; /* Number of individuals who are not selected into candidate subset from family "fam"*/
        Info[fam].subR_best=ivector(1,n); /* Best-so-far subset of individuals selected from family "fam" */
        Info[fam].subD_best=dvector(1,n); /* Transformed phenotypic residuals for the individuals in subR_best */
        Info[fam].c_subR_best=ivector(1,familysize[fam]); /* Complement set of subR_best w.r.t R */
        Info[fam].c_subD_best=dvector(1,familysize[fam]); /* Complement sset of subD_best w.r.t. D */
        Info[fam].famn_best=Info[fam].c_famn_best=0; /* Number of individuals selected into the best-so-far subset from family "fam" and its complement set */
        Info[fam].invN=dmatrix(1,n,1,n); /* Inverse of kinship matrix for individuals selected into candidate subset from family "fam" */
        Info[fam].phiN=dmatrix(1,n,1,n); /* Kinship matrix for individuals selected into candidate subset from family "fam" */
        Info[fam].E=dvector(1,familysize[fam]); /* Enrichment effect for phenotyped individuals in family "fam" */
        Info[fam].famcomp=dvector(1,3); /* Family component used for calculation of the noncentrality parameter */
        Info[fam].matrix=dmatrix(1,familysize[fam],1,familysize[fam]); /* Kinship matrix for phenotyped individuals in family "fam" */
        //Info[fam].Cov=dmatrix(1,familysize[f],1,ncov+1);
        Info[fam].enrich_select=ivector(1,familysize[fam]);
    }
   
    if((fpheno=fopen(phenofile,"r"))==NULL){
        printf("Cannot open phenotype file.\n");
        exit(1);
    }
    
    /* read in the phenotype and kinship file */
    readin(fpheno,fkinship,Info,f,&totalgeno,&n, trait_type,prevalance,ncov,nindiv,fpheno_MASTOR,fgeno_MASTOR,famIDvector);
 
    if(trait_type==2){
       // char cwd[1024];
      //  getcwd(cwd, sizeof(cwd));
        
    char *s1="./mastor_tpr -p phenotype_MASTOR.txt -g genotype_MASTOR.txt -k ";
    char *s3=" --tpr";
    char *result = malloc(strlen(s1)+strlen(kinship)+strlen(s3)+1);
    
    //strcpy(result,cwd);
    strcpy(result,s1);
    strcat(result,kinship);
    strcat(result,s3);
   
        
    printf("Calculating the transformed phenotypic residual using MASTOR_tpr program...\n");
    system(result);
        
    ftpr=fopen("TPR.txt","r");
    MatchTPR(fpheno,ftpr,famIDvector,Info,f);
        fclose(ftpr);
        fclose(fpheno);
    }
    
    
/**********************************************************************************/
 
 
    

/*prepare enrichment value for each feasible individual*/
 struct sortarray array[totalgeno];
 i=0;
 for(fam=1;fam<=f;fam++){
 for (j=1;j<=Info[fam].r;j++){
 array[i].select=0;
 Info[fam].E[j]=0;
 
 for (k=1;k<=Info[fam].r;k++){
 Info[fam].E[j]+=Info[fam].matrix[j][k]*Info[fam].D[k];
 }
 
 if(Info[fam].G[j]==1){
 array[i].fam=fam;
 array[i].ind=j;
 array[i].value=Info[fam].E[j];
 i++;
 }
 }
 }
 
if(strategy_option==2){
    fenrich=fopen("g_enrichment.txt","wb");
   
 Sortarray(array,totalgeno);
 
 for(fam=1;fam<=f;fam++){
 for (j=1;j<=Info[fam].r;j++){
 if(Info[fam].fix[j]==1){
 Info[fam].enrich_select[j]=1;
 }
 else
 Info[fam].enrich_select[j]=0;
 }
 }
 
 
 int ntotal=0;
 int count=0;
    
 while(ntotal<n-totalfix){
 if (count%2==0){
 if(Info[array[(count/2)].fam].fix[array[(count/2)].ind]==0){
Info[array[(count/2)].fam].enrich_select[array[(count/2)].ind]=1;
 if(Info[array[(count/2)].fam].G[array[(count/2)].ind]==1)
 ntotal++;
 }
 count++;
 }else if(count%2==1){
 if(Info[array[totalgeno-(count-1)/2-1].fam].fix[array[totalgeno-(count-1)/2-1].ind]==0){
     Info[array[totalgeno-(count-1)/2-1].fam].enrich_select[array[totalgeno-(count-1)/2-1].ind]=1;
 if(Info[array[totalgeno-(count-1)/2-1].fam].G[array[totalgeno-(count-1)/2-1].ind]==1)
 ntotal++;
 }
 count++;
 }
 }
 
     
fprintf(fenrich,"famID indID enrichment_value previously_genotyped_flag\n");
 
 for(fam=1;fam<=f;fam++){
 for (j=1;j<=Info[fam].r;j++){
 if( Info[fam].enrich_select[j]==1){
 
 fprintf(fenrich,"%d\t%d\t%lf\t%d\n",Info[fam].famID,Info[fam].indID[j],Info[fam].E[j],Info[fam].fix[j]);
 }
 }
 }
 
 
 
 fflush(fenrich);
 fclose(fenrich);
 printf("Extreme enrichment selection is done.\nOutput is in file: g_enrichment.txt\n");
}

    
 
 /************************************************/

else{/* Workspace variables; preparing for random sampling */
    double *tem_x[(n+1)];
    double *tem_m[(n+1)];
    temR=ivector(1,totalgeno-totalfix);
    tem_subR2=ivector(1,nindiv);
    tem_subR=ivector(1,n);
    
    
    /*
     temR - a vector storing new Id's for individuals available to be selected. The new Id is a set of consecutive integers starting from 0, which can be used as the address in an array to facilitate the need for random sampling.
     tem_subR - candidate subset of individuals to be selected based on the new Id's; tem_subR[1],...,tem_subR[totalfix] are ''initially genotyped indivdiuals''; tem_subR[totlafix+1],...,tem_subR[n] are ''additional genotyeped indivdiuals''.
     tem_subR2 - a 0-1 vector labeling whether or not the individual is in tem_subR; -9 means unavailable for selection.
     */
    
    k=1;q=0;j=1;
    
    for (fam=1;fam<=f;fam++){
        for (i=1;i<=Info[fam].r;i++){
            q++;
            if ((Info[fam].fix[i]==0) && (Info[fam].G[i]==1)){
                temR[k]=q;
                tem_subR2[q]=0;
                k++;
            } else if ((Info[fam].fix[i]==1) && (Info[fam].G[i]==1)){
                tem_subR2[q]=1;
                tem_subR[j]=q;
                j++;
            }
            else if ((Info[fam].fix[i]==1) && (Info[fam].G[i]==0))
            {
                
                printf("Invalid input! Individuals who are previously genotyped should also be considered available for selection!");
                exit(1);
            }
            else
                tem_subR2[q]=-9;
        }
    }

    
    /* output to files */
    fp_R=fopen(outfile,"wb");
    
    
    fp_p=fopen("objective_trajectory.txt","w");
    fp_m=fopen("accept_rates.txt","w");
    fp_hp=fopen("objective_up.txt","w");
    
    
    
    /* initial seed */
    if ((seedfl = fopen("seed.txt","r")) == NULL){
        seed = -1 * (long) time(NULL);
    }
    else {
        fscanf(seedfl,"%ld",&seed);
        fclose(seedfl);
    }

    
    srand((unsigned int)seed);
    
    
    /* initialize the family components in preparation for noncentrality parameter calculation */
    famcomp_sum=dvector(1,3);
    famcomp_sumnew=dvector(1,3);
    famcomp_pre=dmatrix(1,2,1,3);
    
    for (j=1;j<=3;j++){
        famcomp_sum[j]=0;
        famcomp_sumnew[j]=0;
    }
    
    for (i=1;i<=2;i++){
        for (j=1;j<=3;j++){
            famcomp_pre[i][j]=0;
        }
    }
    
    for (i=1;i<=n;i++){
        tem_m[i]=dvector(0,i*i-1);
        tem_x[i]=dvector(0,i*i-1);
    }
    
    
   // begin=clock();
    /* initialize the temperature by taking the maximum change in noncentrality parameters among, say, 500 random pairs of neighboring candidate subsets */
   
                                             
    q=500;
    
    tem=tune(Info,tem_subR, tem_subR2,totalgeno, f, n,  temR, tem_m, tem_x, famcomp_sum, famcomp_sumnew,famcomp_pre, q, tem1, tem2, ind_r, ind_a,fam_r, fam_a,totalfix);
    
   // printf("initial temperature is %f\n",tem);
    
    
    /* initialize the candidate subset by random sampling */
    sample(tem_subR,tem_subR2, totalgeno-totalfix,n-totalfix,temR,totalfix);
    famsample(tem_subR,tem_subR2,totalgeno,n,Info);
    
    
    
    /* calculate the corresponding noncentrality parameter of the initial candidate subset */
    famcomponent(Info, tem_m,tem_x,f,famcomp_sum);
    pbest=powersubR=famcomp_sum[1]-(famcomp_sum[2]*famcomp_sum[2])/famcomp_sum[3];
    
    
    
    printf("initial objective function value %lf\n",powersubR);
    
    
    /* update the best-so-far subset by coping from the current candidate subset */
    update(Info,f,0);
    
    
    /* Search the optimal candidate subset by simulated annealing */
    
    
    /* Controlling parameters in annealing schedule */
    alpha=0.95; /* Rate of geometric cooling of temperature */
    nover=5000;
    nlimit=0.5*5000;
    
    
    /* Stopping Criteria: until temperature decreases below some threshold, say, 1.0e-5, or no successful moves are found at some temperature stage. */
    while((tem>STOPTEM)&&(nsucc>para_opt)){
        
        
        nsucc=0; /* Number of successful moves at current temperature stage */
        
        tem=tem*alpha;  /* Temperature stage based on geometric cooling schedule */
        
        
        /* Repeat the following at each temperature stage */
        for (j=1;j<=nover;j++){
            
            
            fprintf(fp_p,"%lf\n",powersubR);
            fflush(fp_p);
            
            
            
            /* propose a move from current candidate subset to its neighboring set, by taking one person out and replace by a different person in feasible set*/
            i=1,k=1;
            ind_r=equal(n-totalfix);
            ind_a=equal(totalgeno-n);
            
            while((ind_r=ind_r-(Info[i].famn-Info[i].nfix))>0)
            {i++;}
            ind_r=ind_r+Info[i].famn;
            
            while((ind_a=ind_a-Info[k].c_famn)>0)
            {k++;}
            ind_a=ind_a+Info[k].c_famn;
            
            fam_r=i;
            fam_a=k;
            
            neigh(Info, &ind_r, &fam_r, &ind_a, &fam_a, &tem1, &tem2);
            
            /* calculate the change in noncentrality parameter between the two candidate subsets */
            de=diffpower(Info, famcomp_sum,famcomp_sumnew,famcomp_pre,&ind_r, &fam_r, &ind_a, &fam_a, tem_m,tem_x);
            
            
            /* Probabilistically decides whether or not to accept the move based on Metropolis rule */
            ans=metrop(de,tem);
            
            
            /* if accepted, call this step a successful move; add 1 to nsucc. */
            if (ans){
                
                nsucc++;
                
                /* update workspace variables; prepare for the next step */
                if (fam_r==fam_a){
                    Info[fam_a].c_subR[ind_a]=tem1;
                    Info[fam_a].c_subD[ind_a]=tem2;
                }
                else {
                    Info[fam_r].c_subR[Info[fam_r].c_famn]=tem1;
                    Info[fam_r].c_subD[Info[fam_r].c_famn]=tem2;
                    if (ind_a<Info[fam_a].c_famn+1){
                        Info[fam_a].c_subR[ind_a]=Info[fam_a].c_subR[Info[fam_a].c_famn+1];
                        Info[fam_a].c_subD[ind_a]=Info[fam_a].c_subD[Info[fam_a].c_famn+1];
                    }
                }
                
                for (i=1;i<=3;i++){
                    famcomp_sum[i]=famcomp_sumnew[i];
                }
                
                powersubR=famcomp_sum[1]-(famcomp_sum[2]*famcomp_sum[2])/famcomp_sum[3];
                
                
                /* update the best-so-far subset*/
                if(powersubR>pbest){
                    pbest=powersubR;
                    update(Info,f,0);
                }
                
                
            }
            
            /* if rejected, undo the changes and prepare for the next step */
            else{
                unneigh(Info,famcomp_pre,&fam_r,&ind_r,&fam_a,&ind_a,&tem1, &tem2);
            }
            
            
            if (nsucc>=nlimit) break;
            
        }
        
        
        
        fprintf(fp_m,"%lf\n",(double)nsucc/j);
        fflush(fp_m);
        
        
        
    }

    
    
    /* restore the best-so-far subset from simulated annealing */
    powersubR=pbest;
    update(Info,f,1);
    famcomponent(Info, tem_m,tem_x,f,famcomp_sum);
    pbest=powersubR=famcomp_sum[1]-(famcomp_sum[2]*famcomp_sum[2])/famcomp_sum[3];
    printf("objective function value after simulated annealing %f\n",pbest);
    fprintf(fp_hp,"%lf\n",pbest);
    fflush(fp_hp);
    
    
    /* post-proceed the best-so-far subset by running a locally up-climbing optimization; ensure the selected subset is a local maximum */
    while(uphill(Info, f, &powersubR, famcomp_sum, famcomp_sumnew, famcomp_pre,tem_m, tem_x,  &pbest, &tem1, &tem2)==1){
        update(Info,f,1);
        famcomponent(Info, tem_m,tem_x,f,famcomp_sum);
        pbest=powersubR=famcomp_sum[1]-(famcomp_sum[2]*famcomp_sum[2])/famcomp_sum[3];
        fprintf(fp_hp,"%lf\n",pbest);
        fflush(fp_hp);
    }

    fprintf(fp_R,"famID indID enrichment_value previously_genotyped_flag\n");
    writevector(fp_R,Info,f);
    fflush(fp_R);
    
    
    
    printf("objective function value after optimization %f\n",pbest);
    
    printf("G-STRATEGY selection is done.\nOutput is in file: g_strategy.txt\n");
    
    //end = clock();
    //time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("time spent in total %lf secs\n", time_spent);
  

    
    free_ivector(temR,1, totalgeno-n);
    free_ivector(tem_subR,1,n);
    free_ivector(tem_subR2,1,nindiv);
    free_dvector(famcomp_sum,1,3);
    free_dvector(famcomp_sumnew,1,3);
    free_dmatrix(famcomp_pre,1,2,1,3);
    for (i=1;i<=n;i++){
        free_dvector(tem_m[i],0,i*i-1);
        free_dvector(tem_x[i], 0, i*i-1);
    }
    
    seedfl=fopen("seed.txt","w");
    fprintf(seedfl,"%ld\n",-1*seed);
    fclose(seedfl);
    fclose(fp_R);
    fclose(fp_m);
    fclose(fp_p);
    fclose(fp_hp);
    
    
}
    
    
    for (fam=1;fam<=f;fam++){
        free_ivector(Info[fam].R,1,familysize[fam]);
        free_dvector(Info[fam].D,1,familysize[fam]);
        free_ivector(Info[fam].subR,1,n);
        free_dvector(Info[fam].subD,1,n);
        free_ivector(Info[fam].c_subR,1,familysize[fam]);
        free_dvector(Info[fam].c_subD,1,familysize[fam]);
        free_ivector(Info[fam].fix,1,familysize[fam]);
        free_ivector(Info[fam].subR_best,1,n);
        free_dvector(Info[fam].subD_best,1,n);
        free_ivector(Info[fam].c_subR_best,1,familysize[fam]);
        free_dvector(Info[fam].c_subD_best,1,familysize[fam]);
        free_dmatrix(Info[fam].invN,1,n,1,n);
        free_dmatrix(Info[fam].phiN,1,n,1,n);
        free_dvector(Info[fam].E,1,familysize[fam]);
        free_dvector(Info[fam].famcomp,1,3);
        free_dmatrix(Info[fam].matrix,1,familysize[fam],1,familysize[fam]);
        free_ivector(Info[fam].enrich_select,1,familysize[fam]);
    }
    
    
    
    
}


int findInd(int n_indiv, int *indivlist, int indiv){
    int index = 1;
    while (index <= n_indiv){
        if(indivlist[index] == indiv) return (index);
        index++;
    }
    printf("ERROR: Individual %d is present in the kinship file, but not in the phenotype file!\n", indiv);
    exit(1);
}


void readin(FILE *fpheno,FILE *fkinship, struct FAM *Info, int f,int *totalgeno, int *n, int trait_type, double prevalance, int ncov, int nindiv,FILE* fpheno_MASTOR,FILE *fgeno_MASTOR,int *famIDvector){
    
    int famID,indID,indID1,indID2,Geno,Fix,famIDold=0,totalfix=0,fatherID,motherID,Sex;
    double Phe;
    double kin;
    *totalgeno=0;
    int fam=0;
    double Cov;
    
    
    
    int totalcount=0;
    
    fprintf(fgeno_MASTOR,"Chr Marker cm bp");
//readin phenotype file
    while(fscanf(fpheno,"%d %d %d %d %d %d %d %lf",&famID, &indID, &fatherID, &motherID, &Sex, &Fix, &Geno,&Phe)==8){
        
        fprintf(fpheno_MASTOR,"%d %d %d %d %d %lf",famID, indID, fatherID, motherID, Sex, Phe);
        fprintf(fgeno_MASTOR," %d",indID);
       
        if (famID!=famIDold)
        {
            fam++;
            Info[fam].famID=famID;
            famIDold=famID;
            famIDvector[fam]=famID;
            
        }
       
         Info[fam].r++;
        
        if(trait_type==1){
        if(Phe==0)Phe=0;
        else if(Phe==1)Phe=-prevalance;
        else if(Phe==2)Phe=1-prevalance;
            else {
                 printf("ERROR: The column (8) in the phenotype file should be the case-control status (0=unknown, 1=unaffected, 2=affected)!\n");
                exit(1);
            }
        }
        
 
        if(trait_type==2){
          if(ncov>=1){
          for(i=1;i<=ncov;i++){
            fscanf(fpheno,"%lf",&Cov);
            fprintf(fpheno_MASTOR," %lf",Cov);
            //Info[fam].Cov[Info[fam].r][i+1]=Cov;
          }
               fprintf(fpheno_MASTOR,"\n");
          }
            else
                 fprintf(fpheno_MASTOR,"\n");
        }
        
         Info[fam].indID[Info[fam].r]=indID;
         Info[fam].R[Info[fam].r]=Info[fam].r;
         Info[fam].D[Info[fam].r]=Phe;
         Info[fam].fix[Info[fam].r]=Fix;
         Info[fam].G[Info[fam].r]=Geno;
                
                if (Fix==1){
                    Info[fam].nfix++;
                    Info[fam].subR[Info[fam].nfix]=Info[fam].r;
                    Info[fam].subD[Info[fam].nfix]=Phe;
                    
                }
        
        
                if (Geno==1){
                    Info[fam].ngeno++;
                    (*totalgeno)++;
                }
    }
    
    
      fprintf(fgeno_MASTOR,"\n");
    fprintf(fpheno_MASTOR,"\n");
    fclose(fgeno_MASTOR);
fclose(fpheno_MASTOR);
    
    //readin kinship file
    famIDold=0;
    int index1,index2,index1old=0,index2old=0,indID1old,indID2old;
    
    while(fscanf(fkinship, "%d %d %d %lf",&famID, &indID1, &indID2,&kin)==4){
        
        if (famID!=famIDold){
            
            
            fam=findInd(f,famIDvector,famID);
            
            
            index1=findInd(Info[fam].r,Info[fam].indID,indID1);
            
          
            
            index2=findInd(Info[fam].r,Info[fam].indID,indID2);
            
        }
       else{
           
            if(indID1==indID1old){
                index1=index1old;
                
            }
            else if (indID1==indID2old){
                index1=index2old;
                
            }
            else {
                index1=findInd(Info[fam].r,Info[fam].indID,indID1);
                
            }
            
            
            if(indID2==indID2old){
                index2=index2old;
            }
            else if (indID2==indID1old){
                index2=index1old;
            }
            else{
                index2=findInd(Info[fam].r,Info[fam].indID,indID2);
            }
           
        }
        
        if(index1!=index2){
            Info[fam].matrix[index1][index2] = Info[fam].matrix[index2][index1]=2*kin;
        }
        else{
            if(kin>0.1){
                printf("WARNING: Inbreeding coefficient for individual %d in family %d is %lf > 0.1. This is rare for most human samples. ", indID1, famID, kin);
                printf("User should make sure the kinship input file is correctly prepared. In particular, the file should contain inbreeding coefficients, not self-kinship coefficients.\n\n");
            }
            
            Info[fam].matrix[index1][index2]=(1+kin);
            
        }
        famIDold=famID;
        index1old=index1;
        index2old=index2;
        indID1old = indID1;
        indID2old = indID2;
        
    }
    
    
    if (fam!=f){
        printf("ERROR: Number of independent families in the kinship file does not match with that in the size file!\n");
        exit(1);
    }
    
    
}


/***********************************************
 This function checks the input format and reads the size of each family based on the input.
 As the same time, this function counts the number of previously genotyped individuals.
 ********************************************************/
int checkformat(FILE *fpheno, int ncov, int nindiv, int n, int f, int *familysize){
    

 int famID,indID,fatherID,motherID,Fix,Geno,famIDold=0,totalfix=0,totalfeasible=0,Sex;
    
    int totalcount=0;
    int fam=0;
    int i;
    double Cov;
    double Phe;
    
    
   while(fscanf(fpheno,"%d %d %d %d %d %d %d %lf",&famID, &indID, &fatherID, &motherID, &Sex, &Fix, &Geno,&Phe)==8){
       
       if(ncov>=1){
      for(i=1;i<=ncov;i++){
           fscanf(fpheno,"%lf",&Cov);
       }
       }
       
       
        totalcount++;
           
            
            if (famID!=famIDold)
            {
                fam++;
                famIDold=famID;
            }
       
       
       familysize[fam]++;
       
            if (fam>f){
                printf("ERROR: Number of families in the phenotype file does not match with that provided in the size file!\n");
                exit(1);
            }
            
             if (Fix==1){
                  totalfix++;
             }
            if (Geno==1){
                 totalfeasible++;
            }
   }
    
 
//check format
    if (fam!=f){
        printf("ERROR: Number of families in the phenotype file does not match with that provided in the size file!\n");
        exit(1);
    }
    
    if(totalcount!=nindiv){
        printf("Number of individuals in the phenotype file does not match with that provided in the size file!\n");
        exit(1);
    }
    
    if (totalfeasible<n){
        printf("ERROR: The target number of individuals to be genotyped should be smaller than the total number of feasible individuals.\n");
        exit(1);
    }

    if (n<totalfix){
        printf("ERROR: The target number of individuals to be genotyped should be larger than the total number of previously genotyped individuals. \n");
        exit(1);
    }

    if(totalfeasible==n){
        printf("Warning: No selection is performed! The target number of individuals to be genotyped exactly matches with the togal number of feasible indivdiuals.\n");
        exit(1);
    }
    
    if(n==totalfix){
        printf("Warning: No selection is performed! The target number of individuals to be genotyped exactly matches with the togal number of previously genotyped indivdiuals.\n");
        exit(1);

    }

    return(totalfix);
}


/***********************************************
 This function reads in the size file.
 
 nindiv - number of individuals in phenotype file
 n - target number of indivdiuals to be genotyped (#extended genotyped set).
 f - number of families in phenotype file
 ********************************************************/


void readsize(FILE *fsize,int *nindiv, int *n, int *f, int *ncov, double *prevalance,int trait_type){
    
    int count=0;
    double *vector;
    vector=dvector(1,4);
    float val;
    
    while(fscanf(fsize,"%f",&val)==1 && count<5){
        count++;
        vector[count]=val;
    }
    
    
    *nindiv=(int)vector[1];
    *f=(int)vector[2];
    *n=(int)vector[3];
    val=vector[4];
    
    if(trait_type==1){
        if((val>=1)||(val<=0)){
            printf("ERROR: The population prevalance should be an estimate between 0 and 1.\n");
            exit(1);
        }
        else
            *prevalance=val;
    }
    else if(trait_type==2){
        if(ceilf(val)-val>TINY){
            printf("ERROR:  Number of covariates should be a non-negative integer. \n");
            exit(1);
        }
        else
            *ncov=(int)(val);
    }
}

/***** this function writes tpr into the main program **/

 void MatchTPR(FILE *fpheno,FILE *ftpr, int *famIDvector, struct FAM *Info, int f){
 
 int famID,indID,fam,index;
 double tpr, res;

     
     while(fscanf(ftpr,"%d %d %lf %lf",&famID, &indID, &tpr, &res)==4){
     
     fam=findInd(f,famIDvector,famID);
     index=findInd(Info[fam].r,Info[fam].indID,indID);
         
    Info[fam].D[index]=tpr;
      
   // printf("%d\t%d\t%lf\n",fam,index, Info[fam].D[index]);
     }
 }

    



void Sortarray(struct sortarray *array, int n){
    qsort(array, n, sizeof(struct sortarray),comparray);
}

int comparray(const void *i, const void *j)
{

    struct sortarray *i1, *i2;
    i1=(struct sortarray*)i;
    i2=(struct sortarray*)j;
    
    if(i1->value==i2->value)
        return (random2()<0.5);
    
    else if (i1->value<i2->value)
        return -1;
    else
        return 1;
}



double random2()
{
    return (double)rand() / (double)RAND_MAX ;
}




/******************************************************
 copy the selected individuals between best-so-far subset and current candidate subset (reverse=0 means rewrite the best-so-far subset from current candidate subset) or vice verse (reverse=1)
 f - number of families in phenotype file
 *****************************************************/

void update(struct FAM *Info, int f, int reverse){
    int i,j;
    if (reverse==0){
        for (i=1;i<=f;i++){
            Info[i].famn_best=Info[i].famn;
            Info[i].c_famn_best=Info[i].c_famn;
            
            if (Info[i].famn!=0){
                for (j=1;j<=Info[i].famn;j++){
                    Info[i].subR_best[j]=Info[i].subR[j];
                    Info[i].subD_best[j]=Info[i].subD[j];
                }
            }
            if (Info[i].c_famn!=0){
                for (j=1;j<=Info[i].c_famn;j++){
                    Info[i].c_subR_best[j]=Info[i].c_subR[j];
                    Info[i].c_subD_best[j]=Info[i].c_subD[j];
                }
            }
        }
    }
    
    else if (reverse==1){
        for (i=1;i<=f;i++){
            Info[i].famn=Info[i].famn_best;
            Info[i].c_famn=Info[i].c_famn_best;
            
            if (Info[i].famn_best!=0){
                for (j=1;j<=Info[i].famn_best;j++){
                    Info[i].subR[j]=Info[i].subR_best[j];
                    Info[i].subD[j]=Info[i].subD_best[j];
                }
            }
            if (Info[i].c_famn_best!=0){
                for (j=1;j<=Info[i].c_famn_best;j++){
                    Info[i].c_subR[j]=Info[i].c_subR_best[j];
                    Info[i].c_subD[j]=Info[i].c_subD_best[j];
                }
            }
        }
    }
    
}




/******************************************************
 This function allocates the individuals in tem_subR to FAM structure
 
 totalgeno - number of individuals available for selection
 n - number of individuals to be genotyped (#extended genotyped set)
 tem_subR - the subset of individuals to be selected based on new ID (i.e. a set of positive and unique intergers for each individuals across families in the sample); tem_subR[1],...,tem_subR[totalfix] are "initially genotyped inviduals"; tem_subR[totlafix+1],...,total_subR[n] are "additionally genotyped individuasl".
 tem_subR2 - a 0-1 vector labeling whether or not the individual is in tem_subR. 1=yes, 0=no.
 *****************************************************/


void famsample(int *tem_subR,int *tem_subR2,int totalgeno,int n,struct FAM *Info){
    int i,k,j,t;
    
    
    for (i=1;i<=n;i++){
        j=1;
        t=tem_subR[i];
        while((t=t-Info[j].r)>0){
            j++;
        }
        Info[j].famn++;
        Info[j].subR[Info[j].famn]=Info[j].R[t+Info[j].r];
        Info[j].subD[Info[j].famn]=Info[j].D[t+Info[j].r];
    }
    
    
    i=1;j=1;k=0;
    
    while(k<totalgeno-n){
        
        if (tem_subR2[i]==0){
            t=i;
            j=1;
            while((t=t-Info[j].r)>0){
                j++;
            }
            Info[j].c_famn++;
            Info[j].c_subR[Info[j].c_famn]=Info[j].R[t+Info[j].r];
            Info[j].c_subD[Info[j].c_famn]=Info[j].D[t+Info[j].r];
            k++;
        }
        i++;
    }
    
    
    
}




/******************************************************
 This function prepares the calculation of noncentrality parameter when there are multiple families in the sample. Each family has a three-element workspace vector denoted as Info[fam].famcomp.
 f - number of families
 famcomp_sum - a vector summarizing all components across families and preparing for the calculation of noncentrality parameter
 tem_m, tem_x - workspace vectors used in matrix inversion
 *****************************************************/
void famcomponent(struct FAM *Info, double *tem_m[],double *tem_x[],int f, double *famcomp_sum){
    int i,j,k;
    
    famcomp_sum[1]=famcomp_sum[2]=famcomp_sum[3]=0;
    
    
    
    
    for (i=1;i<=f;i++){
        
        
        for (j=1;j<=Info[i].famn;j++){
            for (k=1;k<=Info[i].famn;k++){
                Info[i].phiN[j][k]=Info[i].matrix[Info[i].subR[j]][Info[i].subR[k]];
            }
        }
        
        
        Info[i].famcomp[1]=Info[i].famcomp[2]=Info[i].famcomp[3]=0;
        
        
        
        
        if (Info[i].famn!=0){
            
            
            solve(Info[i].invN,Info[i].phiN,Info[i].famn,tem_m[Info[i].famn],tem_x[Info[i].famn]);
            
            for (j=1;j<=Info[i].famn;j++){
                Info[i].E[j]=0;
                
                for (k=1;k<=Info[i].r;k++)
                    Info[i].E[j]+=Info[i].matrix[Info[i].subR[j]][k]*Info[i].D[k];
                
            }
            
            
            for (j=1;j<=Info[i].famn;j++){
                
                for (k=1;k<=Info[i].famn;k++){
                    Info[i].famcomp[1]+=Info[i].E[k]*Info[i].invN[j][k]*Info[i].E[j];
                    Info[i].famcomp[2]+=Info[i].invN[j][k]*Info[i].E[j];
                    Info[i].famcomp[3]+=Info[i].invN[j][k];
                    
                }
                
                
            }
            
            
        }
        
        famcomp_sum[1]+=Info[i].famcomp[1];
        famcomp_sum[2]+=Info[i].famcomp[2];
        famcomp_sum[3]+=Info[i].famcomp[3];
        
    }
    
}



/******************************************************
 This function calculates the change in noncentrality parameter between two neighboring candidate subsets.
 fam_r, ine_r - Family and Individual id removed in proposed set
 fam_a, ind_a - Family and Individual Id added in proposed set
 famcomp_sum, famcomp_sumnew, famcomp_sumpre - a three-component vector used to calculate noncentrality parameter
 tem_m, tem_x - workspace vector used in matrix inversion
 *****************************************************/

double diffpower(struct FAM *Info, double *famcomp_sum,double *famcomp_sumnew,double **famcomp_pre,int *ind_r, int *fam_r, int *ind_a, int *fam_a, double *tem_m[],double *tem_x[]){
    
    int i,j,k;
    double de;
    
    
    Info[*fam_r].E[*ind_r]=0;
    
    
    
    for (i=1;i<=3;i++){
        famcomp_pre[1][i]=Info[*fam_r].famcomp[i];
        famcomp_pre[2][i]=Info[*fam_a].famcomp[i];
        Info[*fam_r].famcomp[i]=0;
        Info[*fam_a].famcomp[i]=0;
    }
    
    
    
    if (*fam_r==*fam_a){
        
        
        
        for (k=1;k<=Info[*fam_r].famn;k++){
            Info[*fam_r].phiN[*ind_r][k]=Info[*fam_r].phiN[k][*ind_r]=Info[*fam_r].matrix[Info[*fam_r].subR[*ind_r]][Info[*fam_r].subR[k]];
        }
        
        
        for (k=1;k<=Info[*fam_r].r;k++){
            Info[*fam_r].E[*ind_r]+=Info[*fam_r].matrix[Info[*fam_r].subR[*ind_r]][k]*Info[*fam_r].D[k];
        }
        
        solve(Info[*fam_r].invN,Info[*fam_r].phiN,Info[*fam_r].famn,tem_m[Info[*fam_r].famn],tem_x[Info[*fam_r].famn]);
        
        for (j=1;j<=Info[*fam_r].famn;j++){
            for (k=1;k<=Info[*fam_r].famn;k++){
                Info[*fam_r].famcomp[1]+=Info[*fam_r].E[k]*Info[*fam_r].invN[j][k]*Info[*fam_r].E[j];
                Info[*fam_r].famcomp[2]+=Info[*fam_r].invN[j][k]*Info[*fam_r].E[j];
                Info[*fam_r].famcomp[3]+=Info[*fam_r].invN[j][k];
            }
        }
        
        
        
        
        
        famcomp_sumnew[1]=Info[*fam_r].famcomp[1]-famcomp_pre[1][1]+famcomp_sum[1];
        famcomp_sumnew[2]=Info[*fam_r].famcomp[2]-famcomp_pre[1][2]+famcomp_sum[2];
        famcomp_sumnew[3]=Info[*fam_r].famcomp[3]-famcomp_pre[1][3]+famcomp_sum[3];
        
        de=famcomp_pre[1][1]-Info[*fam_r].famcomp[1]-((famcomp_sum[2]*famcomp_sum[2])/famcomp_sum[3]-(famcomp_sumnew[2]*famcomp_sumnew[2])/famcomp_sumnew[3]);
        
        
        
        
        
        
    }
    else {
        
        Info[*fam_a].E[Info[*fam_a].famn]=0;
        
        for (k=1;k<=Info[*fam_a].famn;k++){
            Info[*fam_a].phiN[Info[*fam_a].famn][k]=Info[*fam_a].phiN[k][Info[*fam_a].famn]=Info[*fam_a].matrix[Info[*fam_a].subR[Info[*fam_a].famn]][Info[*fam_a].subR[k]];
        }
        
        
        for (k=1;k<=Info[*fam_a].r;k++){
            Info[*fam_a].E[Info[*fam_a].famn]+=Info[*fam_a].matrix[Info[*fam_a].subR[Info[*fam_a].famn]][k]*Info[*fam_a].D[k];
        }
        
        
        
        solve(Info[*fam_a].invN,Info[*fam_a].phiN,Info[*fam_a].famn,tem_m[Info[*fam_a].famn],tem_x[Info[*fam_a].famn]);
        
        for (j=1;j<=Info[*fam_a].famn;j++){
            for (k=1;k<=Info[*fam_a].famn;k++){
                Info[*fam_a].famcomp[1]+=Info[*fam_a].E[k]*Info[*fam_a].invN[j][k]*Info[*fam_a].E[j];
                Info[*fam_a].famcomp[2]+=Info[*fam_a].invN[j][k]*Info[*fam_a].E[j];
                Info[*fam_a].famcomp[3]+=Info[*fam_a].invN[j][k];
            }
            
        }
        
        
        
        if (Info[*fam_r].famn!=0){
            
            if (*ind_r<Info[*fam_r].famn+1){
                for (k=1;k<=Info[*fam_r].famn;k++){
                    Info[*fam_r].phiN[*ind_r][k]=Info[*fam_r].phiN[k][*ind_r]=Info[*fam_r].matrix[Info[*fam_r].subR[*ind_r]][Info[*fam_r].subR[k]];
                    
                }
                
                
                for (k=1;k<=Info[*fam_r].r;k++){
                    Info[*fam_r].E[*ind_r]+=Info[*fam_r].matrix[Info[*fam_r].subR[*ind_r]][k]*Info[*fam_r].D[k];
                }
                
            }
            
            solve(Info[*fam_r].invN,Info[*fam_r].phiN,Info[*fam_r].famn,tem_m[Info[*fam_r].famn],tem_x[Info[*fam_r].famn]);
            
            for (j=1;j<=Info[*fam_r].famn;j++){
                for (k=1;k<=Info[*fam_r].famn;k++){
                    Info[*fam_r].famcomp[1]+=Info[*fam_r].E[k]*Info[*fam_r].invN[j][k]*Info[*fam_r].E[j];
                    Info[*fam_r].famcomp[2]+=Info[*fam_r].invN[j][k]*Info[*fam_r].E[j];
                    Info[*fam_r].famcomp[3]+=Info[*fam_r].invN[j][k];
                }
            }
            
            
            
            
            
        }
        
        
        
        
        famcomp_sumnew[1]=famcomp_sum[1]+((Info[*fam_r].famcomp[1]-famcomp_pre[1][1])+(Info[*fam_a].famcomp[1]-famcomp_pre[2][1]));
        famcomp_sumnew[2]=famcomp_sum[2]+((Info[*fam_r].famcomp[2]-famcomp_pre[1][2])+(Info[*fam_a].famcomp[2]-famcomp_pre[2][2]));
        famcomp_sumnew[3]=famcomp_sum[3]+((Info[*fam_r].famcomp[3]-famcomp_pre[1][3])+(Info[*fam_a].famcomp[3]-famcomp_pre[2][3]));
        
        
        de=-(Info[*fam_r].famcomp[1]-famcomp_pre[1][1])-(Info[*fam_a].famcomp[1]-famcomp_pre[2][1])-((famcomp_sum[2]*famcomp_sum[2])/famcomp_sum[3]-(famcomp_sumnew[2]*famcomp_sumnew[2])/famcomp_sumnew[3]);
        
        
        
    }
    return de;
}




/******************************************************
 Print out the current candidate subset
 f - number of families
 Info[fam].subR - individual Id's in candidate subset; grouped by family index "fam".
 sort - a 0-1 command to sort the individuals in candidate set based on the id's within each family. 1 = increasing order; otherwise = unsorted
 *****************************************************/

void print_ivector(int *vector,int n,int sort){
    int i;
    if (sort==1){
        int *tem=ivector(0,n-1);
        
        for (i=1;i<=n;i++){
            tem[i-1]=vector[i];
        }
        qsort(tem,n,sizeof(int), comp);
        
        for (i=0;i<=n-1;i++)
            printf("%d\t", tem[i]);
        printf("\n");
        free_ivector(tem,0,n-1);
    }else{
        
        for (i=1;i<=n;i++)
            printf("%d\t", vector[i]);
        printf("\n");
        
        
    }
    
}




void print_dvector(double *vector,int n){
    int i;
    for (i=1;i<=n;i++)
        printf("%f\t", vector[i]);
    printf("\n");
}


int comp(const void *i, const void *j)
{
    return *(int *)i - *(int *)j;
}




/******************************************************
 print out the current candidate subset
 f - number of families
 Info[fam].subR -  the individual Id's in current candidate subset
 *****************************************************/
void print_vector(struct FAM *Info, int f){
    int i,j;
    for (i=1;i<=f;i++){
        
        if (Info[i].famn!=0){
            printf("fam %d\t ind\t",i);
            print_ivector(Info[i].subR,Info[i].famn,1);
        }
    }
    printf("number of genotyped individuals from each family\n");
    
    for (i=1;i<=f;i++){
        printf("%d\t", Info[i].famn);
    }
    
    printf("\n");
    
    
}




/******************************************************
 This function initializes the temperature by taking the max changes in noncentrality parameter among, say, 500 random pairs of neighboring candidate sets
 q - the number of of random pairs searched
 ******************************************************/
double tune(struct FAM *Info,int *tem_subR, int *tem_subR2, int totalgeno, int f, int n, int *temR, double *tem_m[], double *tem_x[], double *famcomp_sum, double *famcomp_sumnew, double **famcomp_pre, int nq, int tem1, double tem2, int ind_r, int ind_a, int fam_r, int fam_a,int totalfix){
    
    int j,i,k,fam;
    double de, debest=0;
    int q;
    
    for (j=1;j<=nq;j++){
        
        sample(tem_subR,tem_subR2, totalgeno-totalfix,n-totalfix,temR,totalfix);
        
        famsample(tem_subR,tem_subR2,totalgeno,n,Info);
        
        
        
        
        q=0;k=1;
        
        for (fam=1;fam<=f;fam++){
            for (i=1;i<=Info[fam].r;i++){
                q++;
                if((Info[fam].fix[i]==0)&&(Info[fam].G[i]==1)) {
                    temR[k]=q;
                    tem_subR2[q]=0;
                    k++;
                }
                
                else if (Info[fam].fix[i]==1){
                    tem_subR2[q]=1;
                }
                
            }
        }
        
        famcomponent(Info, tem_m,tem_x,f,famcomp_sum);
        
        
        i=1,k=1;
        ind_r=equal(n-totalfix);
        
        ind_a=equal(totalgeno-n);
        
        while((ind_r=ind_r-(Info[i].famn-Info[i].nfix))>0)
        {i++;}
        ind_r=ind_r+Info[i].famn;
        
        while((ind_a=ind_a-Info[k].c_famn)>0)
        {k++;}
        ind_a=ind_a+Info[k].c_famn;
        
        fam_r=i;
        fam_a=k;
        
        
        neigh(Info, &ind_r, &fam_r, &ind_a, &fam_a, &tem1, &tem2);
        
        
        de=diffpower(Info, famcomp_sum,famcomp_sumnew,famcomp_pre,&ind_r, &fam_r, &ind_a, &fam_a, tem_m,tem_x);
        
        if (fabs(de)>debest)debest=fabs(de);
        
        
        famcomp_sum[1]=famcomp_sum[2]=famcomp_sum[3]=0;
        for (i=1;i<=f;i++){
            Info[i].c_famn=Info[i].famn=0;
            for (k=1;k<=Info[i].famn+1;k++){
                Info[i].subR[k]=0;
            }
        }
        
    }
    
    
    return debest;
    
}





/******************************************************
 Metropolis rule: the acceptance probability depends on the change of noncentrality parameter between the pair of neighboring candidate subsets "de", and a time-varing temperature "tem".
 *****************************************************/

int metrop(double de, double tem){
    if(de<0)
        return 1;
    else if ((de>=0) && (random2() < exp(-de/tem)))
        return 1;
    else
        return 0;
}



/******************************************************
 This function post-preceeds the best-so-far candidate subset by running a local up-climbing until no better neighbors are found. It ensures the resulting subset is a local maximum.
 f - number of families
 powersubR - noncentrality parameter for the current candidate subset
 famcomp_sum, famcomp_sumnew, famcomp_sumpre - a three-component vector used to calculate noncentrality parameter
 tem_m, tem_x - workspace vectors used in matrix inversion
 pbest - optimal noncentrality parameter after post-processing.
 tem1, tem2 -  workspace variables storing the Id and transformed phenotypic residual of the individual added or removed in some intermediate step.
 *****************************************************/
int uphill(struct FAM *Info, int f, double *powersubR, double *famcomp_sum, double *famcomp_sumnew, double **famcomp_pre,double *tem_m[], double *tem_x[],double *pbest, int *tem1, double *tem2){
    
    int fam_r,ind_r,fam_a,ind_a,k,index=0;
    
    
    double de, powersubR_new;
    
    
    for (fam_r=1;fam_r<=f;fam_r++){
        
        for (ind_r=1;ind_r<=Info[fam_r].famn-Info[fam_r].nfix;ind_r++){
            ind_r=ind_r+Info[fam_r].nfix;
            for (fam_a=1;fam_a<=f;fam_a++){
                
                for (ind_a=1;ind_a<=Info[fam_a].c_famn;ind_a++){
                    
                    
                    neigh(Info, &ind_r, &fam_r, &ind_a, &fam_a, tem1, tem2);
                    
                    de=diffpower(Info, famcomp_sum,famcomp_sumnew,famcomp_pre, &ind_r, &fam_r, &ind_a, &fam_a, tem_m,tem_x);
                    
                    
                    
                    powersubR_new=*powersubR-de;
                    
                    if(powersubR_new>*pbest){
                        
                        
                        *pbest=powersubR_new;
                        update(Info,f,0);
                        
                        if (fam_r==fam_a){
                            Info[fam_a].c_subR_best[ind_a]=*tem1;
                            Info[fam_a].c_subD_best[ind_a]=*tem2;
                        }
                        else {
                            
                            
                            Info[fam_r].c_subR_best[Info[fam_r].c_famn]=*tem1;
                            Info[fam_r].c_subD_best[Info[fam_r].c_famn]=*tem2;
                            
                            if (ind_a<Info[fam_a].c_famn+1){
                                Info[fam_a].c_subR_best[ind_a]=Info[fam_a].c_subR_best[Info[fam_a].c_famn+1];
                                Info[fam_a].c_subD_best[ind_a]=Info[fam_a].c_subD_best[Info[fam_a].c_famn+1];
                            }
                            
                            
                        }
                        index++;
                    }
                    
                    unneigh(Info,famcomp_pre,&fam_r,&ind_r,&fam_a,&ind_a,tem1, tem2);
                    
                    
                }
            }
        }
    }
    
    return (index!=0);
}



/***************************************************
 This function rejects the proposed move and undoes any changes of workspace variable in the proposal.
 fam_r, ind_r - family and individual Id removed in the proposal
 fam_a, ind_r - family and individual Id added in the proposal
 tem1 and tem2 - workspace variables storing the Id and transformed phenotypic residual of the removed individual in the proposal
 famcomp_pre - workspace variable storing the family component for the family removed in the proposal
 **********************************************/
void unneigh(struct FAM *Info, double **famcomp_pre,int *fam_r,int *ind_r,int *fam_a,int *ind_a,int *tem1, double *tem2){
    
    int i,k;
    
    
    
    for (i=1;i<=3;i++){
        Info[*fam_r].famcomp[i]=famcomp_pre[1][i];
        Info[*fam_a].famcomp[i]=famcomp_pre[2][i];
    }
    
    
    if (*fam_r==*fam_a){
        
        
        Info[*fam_r].subR[*ind_r]=*tem1;
        Info[*fam_r].subD[*ind_r]=*tem2;
        
        Info[*fam_r].E[*ind_r]=0;
        
        
        for (k=1;k<=Info[*fam_r].famn;k++){
            Info[*fam_r].phiN[*ind_r][k]=Info[*fam_r].phiN[k][*ind_r]=Info[*fam_r].matrix[*tem1][Info[*fam_r].subR[k]];
        }
        for (k=1;k<=Info[*fam_r].r;k++){
            Info[*fam_r].E[*ind_r]+=Info[*fam_r].matrix[*tem1][k]*Info[*fam_r].D[k];
        }
        
        
    }
    else {
        
        
        Info[*fam_r].subR[Info[*fam_r].famn+1]=Info[*fam_r].subR[*ind_r];
        Info[*fam_r].subD[Info[*fam_r].famn+1]=Info[*fam_r].subD[*ind_r];
        
        
        
        
        if (*ind_r<Info[*fam_r].famn+1){
            Info[*fam_r].subR[*ind_r]=*tem1;
            Info[*fam_r].subD[*ind_r]=*tem2;
        }
        
        
        
        Info[*fam_r].famn=Info[*fam_r].famn+1;
        Info[*fam_a].famn=Info[*fam_a].famn-1;
        Info[*fam_r].c_famn=Info[*fam_r].c_famn-1;
        Info[*fam_a].c_famn=Info[*fam_a].c_famn+1;
        Info[*fam_r].E[Info[*fam_r].famn]=0;
        Info[*fam_r].E[*ind_r]=0;
        
        
        
        
        for (k=1;k<=Info[*fam_r].famn;k++){
            Info[*fam_r].phiN[Info[*fam_r].famn][k]=Info[*fam_r].phiN[k][Info[*fam_r].famn]=Info[*fam_r].matrix[Info[*fam_r].subR[Info[*fam_r].famn]][Info[*fam_r].subR[k]];
            Info[*fam_r].phiN[*ind_r][k]=Info[*fam_r].phiN[k][*ind_r]=Info[*fam_r].matrix[Info[*fam_r].subR[*ind_r]][Info[*fam_r].subR[k]];
        }
        
        
        
        
        
        if (*ind_r<Info[*fam_r].famn){
            for (k=1;k<=Info[*fam_r].r;k++){
                Info[*fam_r].E[*ind_r]+=Info[*fam_r].matrix[Info[*fam_r].subR[*ind_r]][k]*Info[*fam_r].D[k];
                Info[*fam_r].E[Info[*fam_r].famn]+=Info[*fam_r].matrix[Info[*fam_r].subR[Info[*fam_r].famn]][k]*Info[*fam_r].D[k];
            }
        }
        else{
            for (k=1;k<=Info[*fam_r].r;k++){
                Info[*fam_r].E[Info[*fam_r].famn]+=Info[*fam_r].matrix[Info[*fam_r].subR[Info[*fam_r].famn]][k]*Info[*fam_r].D[k];
            }
        }
        
        
        
        
    }
    
    
}

/***************************************************
 This function moves from candidate subset to its neighboring candidate subset, and updates the workspace variable "Info".
 fam_r, ind_r - family and individual Id to be removed
 fam_a, ind_a - family and individual Id to be added
 tem1 and tem2 - workspace variables storing the Id and transformed phenotypic residual for the removed individual
 **********************************************/

void neigh(struct FAM *Info, int *ind_r, int *fam_r, int *ind_a, int *fam_a, int *tem1, double *tem2){
    
    
    *tem1=Info[*fam_r].subR[*ind_r];
    *tem2=Info[*fam_r].subD[*ind_r];
    
    if (*fam_r==*fam_a){
        Info[*fam_r].subR[*ind_r]=Info[*fam_a].c_subR[*ind_a];
        Info[*fam_r].subD[*ind_r]=Info[*fam_a].c_subD[*ind_a];
        
    }
    else {
        
        if (*ind_r<Info[*fam_r].famn){
            Info[*fam_r].subR[*ind_r]=Info[*fam_r].subR[Info[*fam_r].famn];
            Info[*fam_r].subD[*ind_r]=Info[*fam_r].subD[Info[*fam_r].famn];
        }
        
        Info[*fam_a].subR[Info[*fam_a].famn+1]=Info[*fam_a].c_subR[*ind_a];
        Info[*fam_a].subD[Info[*fam_a].famn+1]=Info[*fam_a].c_subD[*ind_a];
        
        Info[*fam_r].famn=Info[*fam_r].famn-1;
        Info[*fam_a].famn=Info[*fam_a].famn+1;
        
        Info[*fam_r].c_famn=Info[*fam_r].c_famn+1;
        Info[*fam_a].c_famn=Info[*fam_a].c_famn-1;
        
    }
}






/***************************************************
 This function randomly selects r individuals from a collection of n individuals without replacement.
 temR - a vector of n elements to be selected from.
 tem_subR - the subset of individuals to be selected: tem_subR[1],...,tem_subR[totalfix] are "initially genotyped individuals"; tem_subR[totlafix+1],...,tem_subR[n] are "additionally genotyped individuals".
 tem_subR2 - a 0-1 vector labeling whether the individual is in tem_subR.
 **********************************************/

void sample(int *tem_subR, int *tem_subR2, int r, int n, int *temR,int totalfix){
    
    int z,z1,i,j=1;
    
    for (i=1;i<=n;i++){
        j=i-1+equal(r-i+1);
         //printf("test %d\n",j);
        
        z=temR[j];
        tem_subR[i+totalfix]=z;
        z1=temR[i];
        temR[i]=z;
        temR[j]=z1;
        tem_subR2[tem_subR[i+totalfix]]=1;
    }
    
    
}




/***************************************************
 Simple random sampling: select one individual from a collection of n individuals using uniform distribution.
 **********************************************/
int equal(int n){
    int k=0, choose;
    int i=1;
    double p,frac;
    
    while(k<1){
        frac=random2();
        p=1/(n-i+1.0);
        if(frac<p){
            choose=i;
            k++;
        }
        i++;
    }
    return choose;
}





/***************************************************
 This function computes the inverse of a symmetric positive definite matrix.
 y - inverse matrix
 a - kinship matrix to be inverted
 n - dimension of the matrix
 m, x - workspace vectors storing the elements of kinship matrix and identity matrix
 **********************************************/

void solve(double **y, double **a,int n, double *m,double *x){
    
    //printf("okay?\n");
    int i,j;
    
    
    for (i=1;i<=n;i++){
        for (j=1;j<=n;j++)
            m[(i-1)*n+j-1]=a[i][j];
    }
    
    for (i=1;i<=n;i++){
        for (j=1;j<=n;j++)
            x[(i-1)*n+j-1]=(i==j);
    }
    
    //clapack_dposv(CblasRowMajor, CblasUpper, n, n, m, n, x, n);
    
    integer tmp;
    integer nn=n;
    dposv_("U",&nn,&nn,m,&nn,x,&nn,&tmp);
    
    for (i=0;i<=(n-1);i++){
        for (j=0;j<=(n-1);j++)
            y[i+1][j+1]=x[i*n+j];
    }
    
    
   // printf("okay2?\n");
    
}





/***************************************************
 This function writes the current candidate subset into a file called fp
 f - number of families
 Info[fam].subR -  individual Id's in the candidate subset
 **********************************************/

void writevector(FILE *fp, struct FAM *Info, int f){
    
    int maxfamilysize=0;

    
   // int *tem=ivector(0,MAXFAMILYSIZE-1);
    
    int i,j;
    for (i=1;i<=f;i++){
        
        if (Info[i].famn!=0){
            
             int *tem=ivector(0,Info[i].famn-1);
            
            
            for (j=1;j<=Info[i].famn;j++){
                tem[j-1]=Info[i].subR[j];
            }
            
            qsort(tem,Info[i].famn,sizeof(int),comp);
            
            
            for (j=1;j<=Info[i].famn;j++){
                fprintf(fp,"%d\t%d\t%lf\t%d\n",Info[i].famID,Info[i].indID[tem[j-1]], Info[i].E[tem[j-1]],Info[i].fix[tem[j-1]]);
                
                
                // printf("%d\t%d\t%d\n",Info[i].famID,Info[i].indID[tem[j-1]],Info[i].fix[tem[j-1]]);
                fflush(fp);
            }
            
            free_ivector(tem,0,Info[i].famn-1);
            
        }
    }
    
}


