/* OneResLtnt.c Fits univariate response models with latent variable representation
 * Copyright (C) 2014  Georgios Papageorgiou, gpapageo@gmail.com
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details. 
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
 
// NOTES: 1. MAKE LATENT = 0 IF COMPALLOC > NCOMP, 2. MAKE U = 0 IF COMPALLOC > NCOMP
// CHECK CASE OF NO CONTINUOUS COVARIATES  

#define GSL_RANGE_CHECK_OFF
#define HAVE_INLINE
#include <R.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>  
#include <gsl/gsl_roots.h>
#include "matalg.h"
#include "pdfs.h"
#include "sampling.h"
#include "other.functions.h"
#include "mathm.h"
#include "SpecOneRL.h"
//#include "bnsp.h"
#include "cubature.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y)) 
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

extern int (*p_pcubature)(unsigned, integrand, void *, unsigned, const double *, const double *,
	                    size_t, double, double, error_norm, double *, double *);
extern int (*p_hcubature)(unsigned, integrand, void *, unsigned, const double *, const double *,
	                    size_t, double, double, error_norm, double *, double *);


void OneResLtnt(int *seed1, double *X, int *Y, double *H,
                int *sweeps1, int *burn1, int *thin1, int *ncomp1, 
                int *n1, int *p1, int *NBC1, 
                double *Vprior, double *eta1,
                double *mumu, double *sigmamu,
                double *alphaXi, double *betaXi,
                double *Alphaa1, double *Alphab1, double *TruncAlpha1,
                double *xbar, double *xsd, double *Ymean,
                int *family1, int *sampler1,
                int *npred1, double *Xpred, double *Hpred, int *maxy1, int *c,
                double *meanReg, double *medianReg, double *q1Reg, double *q3Reg, double *modeReg,
                double *denReg, double *denVar,
                char **WorkingDir, int *WF1)
{
    gsl_set_error_handler_off();

    // Random number generator initialization
    int seed = seed1[0]; //random number generator seed
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);
    unsigned long int s; //for calling random number generating functions

    // Sweeps, burn-in period, thin
    int sweeps = sweeps1[0]; //number of posterior samples
    int burn = burn1[0]; // integer burn-in period
    int thin = thin1[0]; // integer thin
    int nSamples = (sweeps - burn)/thin;

    //Family
    int family = family1[0]; //1=poisson, 2=binomial, 3=negative binomial, 4=beta binomial, 5=Generalized Poisson, 6=COM-Poisson
    int sampler = sampler1[0]; //1=slice, 2=truncated

    //Declare variables for loops: h for components, k,j for covariates, i for subjects, sw for sweeps
    int h, i, j, k, l, sw;

    // Dimensions
    int n = n1[0]; //number of observations
    int p = p1[0]; //number of covariates
    int NDC = NBC1[0]; //number of discrete covariates
    int NDV = NDC + 1; //number of discrete variables    
    int NCC = p - NDC; //number of continuous covariates
    int totran = p + 1; // number of random variables
    
    int NDCTemp, NCCTemp;
    NDCTemp = NDC;
    if (NDCTemp == 0) NDCTemp += 1;
    NCCTemp = NCC;
    if (NCCTemp == 0) NCCTemp += 1;  
    
    int ncomp = ncomp1[0]; //number of clusters/components
    int npred = npred1[0]; //number of predictions
    int nRespPars = 1; // number of parameters in response pmf
    if (family > 2) nRespPars = 2;
    if (family == 8) nRespPars = 3;
    int XiLoop = 1;
    if (family==5) XiLoop = 2;
    int move, NCCPi, NDCPi, NDVPi; 
    NDVPi = 0; 
    //Tolerance level
    double tol = 0.00001; //tolerance level for eigen values when inverting matrices
    double tol2 = 12.0; //tolerance for integral limits
     
    // Prior parameters
    // - 1 - Sigman^*: V, V^{-1}, df: Eh ~ Wishart(eta,V)
    gsl_matrix *V = gsl_matrix_alloc(totran,totran); //scale matrix of prior of E_h
    gsl_matrix *Vinv = gsl_matrix_alloc(totran,totran); //inverse scale matrix of prior of E_h
    for (j = 0; j < totran; j++)
        for (k = 0; k < totran; k++)
            gsl_matrix_set(V,j,k,Vprior[j*totran+k]);
    gsl_matrix_memcpy(Vinv,V);
    Inverse(totran,Vinv);
    double eta = eta1[0];
    // - 2 - Muh: mean  and Sigma_nu
    gsl_vector *MuMu = gsl_vector_alloc(p);
    gsl_matrix *SigmaMuInv = gsl_matrix_alloc(p,p);
    gsl_matrix *SigmaMuHalf = gsl_matrix_calloc(p,p);
    gsl_vector *SIMmu = gsl_vector_calloc(p);
    for (k = 0; k < p; k++)
        gsl_vector_set(MuMu,k,mumu[k]);
    for (j = 0; j < p; j++)
        for (k = 0; k < p; k++)
            gsl_matrix_set(SigmaMuInv,j,k,sigmamu[j*p+k]);
    Inverse(p,SigmaMuInv);
    gsl_matrix_memcpy(SigmaMuHalf,SigmaMuInv);
    gHalfInv(p,tol,SigmaMuHalf);
    gsl_blas_dgemv(CblasNoTrans,1.0,SigmaMuInv,MuMu,0.0,SIMmu); // Sigma_mu^{-1} Mu_mu
    // - 3 - alphaXi & betaXi: directly use arguments of the function
    // - 4 - concentration parameter alpha ~ Gamma(a,b) I[alpha>c]
    double Alphaa = Alphaa1[0];
    double Alphab = Alphab1[0];
    double TruncAlpha = TruncAlpha1[0];

    // Specify directory
    int WF = WF1[0]; // indicator: 1 = write files, 0 = no files.
    const int MAX_PATH = 300;
    char path_name[MAX_PATH + 1];

    // Open files
    FILE *out_file1, *out_file2, *out_file3, *out_file4, *out_file5, *out_file6, *out_file7, *out_file8, *out_file9,
         *out_file10, *out_file11, *out_file12, *out_file13, *out_file14, *out_file15, *out_file16, *out_file17, 
         *out_file18, *out_file19;

    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Sigmah.txt");
    out_file1 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.muh.txt");
    out_file2 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.xih.txt");
    out_file3 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.alpha.txt");
    out_file4 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.compAlloc.txt");
    out_file5 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.nmembers.txt");
    out_file6 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Updated.txt");
    out_file7 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.MeanReg.txt");
    out_file8 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.MedianReg.txt");
    out_file9 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q05Reg.txt");
    out_file10 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q10Reg.txt");
    out_file11 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q15Reg.txt");
    out_file12 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q20Reg.txt");
    out_file13 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q25Reg.txt");
    out_file14 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q75Reg.txt");
    out_file15 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q80Reg.txt");
    out_file16 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q85Reg.txt");
    out_file17 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q90Reg.txt");
    out_file18 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q95Reg.txt");
    out_file19 = fopen(path_name, "wt");

    // Define quantities of interest
    double Sigmah[ncomp][totran*totran]; //Sigma^*_h
    double muh[ncomp][p]; //mu_h
    double Xi[ncomp][nRespPars]; // response pmf parameters
    double alpha; //concentration parameter
    int compAlloc[n]; // allocation
    int nmembers[ncomp]; // number of cluster members
    double pdih[ncomp][n]; // matrix to store P(delta_j = k_j)

    // Define other quantities
    double pi[ncomp]; //stick breaking probabilities
    double Vdp[ncomp]; //V[i] ~ Beta(1,alpha)
    double newalpha; //store concentration parameter
    int nmb; //number of members
    double u[n]; // for implementing a slice sampler - Walker 2006
    double minU, sumPIh;
    int nUpdated;
    double latentYX[n][NDV]; // continuous latext variables underlying discrete ones
    double Dh[ncomp][totran*totran]; //to store diagonal matrices D_h
    double Eh[ncomp][totran*totran]; //to store E_h^{(t)}
    double XiC[nRespPars];
    double XiP[nRespPars];
    double XiLC, XiLP, cyC, cyP, cymoC, cymoP, Acp, NAcp, DAcp; //MH  step
    NAcp = 0.0; DAcp = 0.0;
    double nuh[ncomp][p]; // to store nu
    double nuSInu, nuSIxmm; //nu^T S^-1 nu and nu^T SI (x-mu)
    double storenuSInu[ncomp]; // to store nu^T Sigma^(-1) nu
    double storenuSI[ncomp][p]; // to store nu^T Sigma^(-1)
    double toths[p];  // total of clusters
    double temp, temp2;
    int temp3;
    double lower[NDV];
    double upper[NDV];    
    double lowerX[NDCTemp];
    double upperX[NDCTemp];
    double Lower, Upper;    
    int dimP = NDV*(NDV+3)/2+1;//params of conditional distribution of latent | continuous variables +1 for det
    double params[dimP]; 
    double StoreParams[ncomp][dimP];
    double val, err; //value and error from cubature function call
    int nocuba; //if limits are two far from zero, integral = 0 and nocubature is needed
    double logPostPdfDSigmaP, logPostPdfDSigmaT, logTransR; //MH CV mat step
    double logAcp, unifRV; //Any MH step
    double lmt = 9999.99;

    // Declare and allocate gsl vectors and matrices
    gsl_matrix *Sh = gsl_matrix_alloc(totran,totran);
    gsl_vector *nuSI = gsl_vector_alloc(p);
    gsl_matrix *Dg1 = gsl_matrix_alloc(p,p);
    gsl_matrix *Dg3 = gsl_matrix_alloc(p,p);
    gsl_matrix *nTDi = gsl_matrix_alloc(p, p);
    gsl_vector *vecZero = gsl_vector_calloc(totran);
    gsl_vector *vecZero2 = gsl_vector_calloc(NCC);
    gsl_vector *Z = gsl_vector_alloc(p); //sampling standard normal
    gsl_vector *U1 = gsl_vector_alloc(p); // matrices to carry out blas multiplications
    gsl_vector *U2 = gsl_vector_alloc(p); //
    gsl_vector *U3 = gsl_vector_alloc(p); //
    gsl_vector *yxstar = gsl_vector_alloc(NDV); //next to impute latent var in initialization
    gsl_vector *CondMean = gsl_vector_calloc(NDV);
    gsl_matrix *PartMean = gsl_matrix_alloc(NDV,NCC); //to become the mean of y^*|(cotinious var)
    gsl_matrix *CondCov = gsl_matrix_alloc(NDV,NDV); //cov mat of y^*|(cotinious var)
    gsl_matrix *Ehp = gsl_matrix_alloc(totran,totran);//proposed random Wishart matrices
    gsl_matrix *Dhp = gsl_matrix_alloc(totran,totran);//D_h^(p) from above decomposition
    gsl_matrix *SigmaShp = gsl_matrix_alloc(totran,totran);  //Sigma_h^*^(p) from above decomposition
    gsl_matrix *SigmaSH = gsl_matrix_alloc(totran,totran);  //Sigma_h^* after updating step
    gsl_matrix *Th = gsl_matrix_alloc(p,p);  //Cop Th
    //gsl_vector *CondMean2 = gsl_vector_calloc(NDCTemp); 
    //gsl_matrix *PartMean2 = gsl_matrix_alloc(NDCTemp,NCC); //to become the mean of x^*|(cotinious covs)
    //gsl_matrix *CondCov2 = gsl_matrix_alloc(NDCTemp,NDCTemp); //cov mat of x^*|(cotinious var)
     
    // GSL Matrix and vector views
    gsl_matrix_view Dg2, nu, Eht, Dht, SigmaSht, SigmaSht2, SigmaSht3, ThCopy;
    gsl_vector_view nu2, tothv, XmMview, Ci, MUSTAR;
    
    //Make adaptive
    int batchL = 50; //batch length
    double WB; //which batch
      //Sigma
    double acceptSigma = 0.0;    
    double SLL = totran+2;
    double SUL = 2*n;
    double DFEhp = SLL+50;  
      //Xi
    double acceptXi[nRespPars];
    for (j = 0; j < nRespPars; j++) 
        acceptXi[j] = 0.0;    
    double XLL = 0.01;
    double XUL = 200;
    double prec[nRespPars];
    for (j = 0; j < nRespPars; j++) 
        prec[j] = 1.0;        
    //Bases for vector and matrix views
    double baseXmM[NCCTemp]; // to view x-muh
    double baseXmM2[NCCTemp]; // to view x-muh
    double baseSigmaSh[totran*totran]; //to view one row of SigmaSh at a time
    double baseDh[totran*totran]; //to view one row of Dh at a time
    double baseEh[totran*totran]; //to view one row of Eh at a time
    double baseC[NCCTemp]; // to view all continuous variables
    double baseMUS[totran]; // to view means of all variables

    //Predictions
    int maxy = maxy1[0];
    double CPDF[npred][ncomp];
    double denomPred[npred];
    double StorePD[npred][maxy];
    double sumProb;
    int mode;
    double cdfy[npred]; // cdf of response to truncate predictions before maxy
    int start;
    double PredTrunc = 0.99999;
    double normConst;
    double StoreUpper[ncomp][maxy];
    double mu, f, B;
    double saveCdf[ncomp];
    int Iend = 0; // end sweep indicator
    double TMR[npred]; //mean regression from each iteration
    double TMdR[npred]; //median regression from each iteration
    double TQ05R[npred]; //05th percentile regression from each iteration
    double TQ10R[npred]; //10th percentile regression from each iteration
    double TQ15R[npred]; //15th percentile regression from each iteration
    double TQ20R[npred]; //20th percentile regression from each iteration
    double TQ25R[npred]; //25th percentile regression from each iteration
    double TQ75R[npred]; //75th percentile regression from each iteration
    double TQ80R[npred]; //80th percentile regression from each iteration
    double TQ85R[npred]; //85th percentile regression from each iteration
    double TQ90R[npred]; //90th percentile regression from each iteration
    double TQ95R[npred]; //95th percentile regression from each iteration

    // Sampler initialization

    // Step 1: Initial allocation to the ncomp components
    for (i = 0; i < n; i++)
        for (h = 0; h < ncomp; h++)
            pdih[h][i] = (double) 1/ncomp; //(ncomp-h) or (ncomp-h)*(ncomp-h);
    s = gsl_ran_flat(r,1.0,100000);
    allocation(s,n,ncomp,pdih,compAlloc,1);

    // Step 2: Find nmembers
    findNmembers(n,ncomp,compAlloc,nmembers);
    int active;
    active = 0;
    for (h = 0; h < ncomp; h++) if (nmembers[h] > 0) active += 1;

    // Step 3: Initialize alpha of the DP
    alpha = Alphaa*Alphab;

    // Step 4: Initialize V_h ~ Beta(nmembers[h]+1,n-nmembers[1:h]+alpha);
    if (sampler == 1){
        nmb = 0;
        for (h = 0; h < ncomp; h++){
            nmb += nmembers[h];
            Vdp[h] = gsl_ran_beta(r, (double) nmembers[h]+1.0, (double) n-nmb+alpha);
            if (h == 0)
                pi[h] = Vdp[h];
            else if (h > 0 && h < (ncomp-1))
                pi[h] = Vdp[h] * pi[h-1] * (1-Vdp[h-1]) / Vdp[h-1];
            else pi[h] = pi[h-1] * (1-Vdp[h-1]) / Vdp[h-1];
        }
    }

    // Step 5: Update the u_i's from unif(0,pi_k_i)
    if (sampler == 1){
        for (i = 0; i < n; i++)
            u[i] = gsl_ran_flat(r,0.0,pi[compAlloc[i]]);
        minU = 1.0;
        for (i = 0; i < n; i++)
            if (u[i] < minU) minU = u[i];
	}
    else{ 
        for (i = 0; i < n; i++)
            u[i] = 0.0;
        minU = -1.0;
    }

    // Step 6: mu_h: initialize covariate cluster means
    s = gsl_ran_flat(r,1.0,100000);
    initMuh(s,X,n,p,ncomp,xbar,xsd,compAlloc,nmembers,muh);

    // Step 7: Initialize parameters Xi_h of the response variables
    s = gsl_ran_flat(r,1.0,100000);
    if (family == 1 || family == 2 || family == 5 || family == 6 || family == 7 || family == 8)
        initGLMOneResLtnt1(s,Y,H,n,ncomp,nRespPars,nmembers,compAlloc,Xi,Ymean[0],family);
    if (family == 3 || family == 4)
        initGLMOneResLtnt2(s,Y,H,n,ncomp,nRespPars,nmembers,compAlloc,Xi,family);

    // Step 8: Inpute latentYX from truncated normal
    gsl_matrix_set_identity(CondCov);
    for (h = 0; h < ncomp; h++){
        for (j = 0; j < nRespPars; j++) XiC[j] = Xi[h][j];
        for (i = 0; i < n; i++){
            if (compAlloc[i] == h){  
                calcGLMLimitsYX(Y,H,X,i,NDC,n,XiC,lower,upper,family);
                gsl_vector_set_zero(yxstar);
                s = gsl_ran_flat(r,1.0,100000);
                sampleTMN(s,NDV,lower,upper,yxstar,CondMean,CondCov,tol);
                for (k = 0; k < NDV; k++)
                    latentYX[i][k] = gsl_vector_get(yxstar,k);
            } 
        } 
    }      

    //Step 9: Initialize E_h, D_h, Sigma_h^* for all h to identity
    initEDS(ncomp,totran,Sigmah,Eh,Dh);
    
    //Rprintf("%f %f %f %f \n",
    //cdf_beta_binomial_P2(10,1,1,1), 
    //cdf_beta_binomial_P2(10,0,1,1)+
    //exp(gsl_sf_lnchoose(10,1) + gsl_sf_lnbeta(1+1,10-1+1) - gsl_sf_lnbeta(1,1)),
    //cdf_beta_binomial_P2(10,4,1,10),
    //cdf_beta_binomial_P2(10,3,1,10)+exp(gsl_sf_lnchoose(10,4) + gsl_sf_lnbeta(4+1,10-4+10) - gsl_sf_lnbeta(1,10)));
    
     //Rprintf("%s %f %f \n", "oops:", cdf_generalized_poisson_P3(0,4.7,0.006,&temp),exp(-4.7/0.006));  
    
    //for (sw = 0; sw < -2; sw++) puts("this");      
    //#############################################SAMPLER 
    for (sw = 0; sw < sweeps; sw++){ 

        if (sw==0) Rprintf("%i %s \n",sw, "posterior samples...");
        if (((sw+1) % 500)==0) Rprintf("%i %s \n",sw+1, "posterior samples...");
        
        modf(sw/batchL,&WB);
        
        //Generate all cluster parameters
        nmb = 0; // counts nmembers up to cluster h
        sumPIh = 0.0; // monitor sum pi_h and stop when it is > 1-minU
        h = 0;
        while ((h < ncomp) && (sumPIh < (1-minU))){

            // - 1 - Update V_h ~ Beta(nmembers[h]+1,n-nmembers[1:h]+alpha);
            //puts("1");
            nmb += nmembers[h];
            Vdp[h] = gsl_ran_beta(r, (double) nmembers[h]+1.0, (double) n-nmb+alpha);

            // - 2 - Update pi_h = V[h] prod_{j < h} 1-V[j]
            //puts("2");
            if (h == 0)
                pi[h] = Vdp[h];
            else if (h > 0 && h < (ncomp-1))
                pi[h] = Vdp[h] * pi[h-1] * (1-Vdp[h-1]) / Vdp[h-1];
            else pi[h] = pi[h-1] * (1-Vdp[h-1]) / Vdp[h-1];
            
            // - 3 - Update Sigma_h^*, D_h and E_h
            //puts("3");
            if ((sw % batchL)==0 && WB > 0 && h==0){ 
	            if (acceptSigma > 0.25 && DFEhp > SLL) DFEhp -= MIN(0.01,1/sqrt(WB)) * DFEhp; 
	            if (acceptSigma <= 0.25 && DFEhp < SUL) DFEhp += MIN(0.01,1/sqrt(WB)) * DFEhp;
		        acceptSigma = 0.0;   
	        }	     
            SetShOneResLtntYX(p,n,NDV,h,ncomp,X,muh,latentYX,compAlloc,Sh);
            for (k = 0; k < totran*totran; k++){
                baseEh[k] = Eh[h][k];
                baseDh[k] = Dh[h][k]; 
                baseSigmaSh[k] = Sigmah[h][k];
            }
            Eht = gsl_matrix_view_array(baseEh,totran,totran);
            Dht = gsl_matrix_view_array(baseDh,totran,totran);
            SigmaSht = gsl_matrix_view_array(baseSigmaSh,totran,totran);            
            gsl_matrix_scale(&Eht.matrix, (double) 1/DFEhp); //proposal
            s = gsl_ran_flat(r,1.0,100000);
            rwish(s,totran,DFEhp,&Eht.matrix,Ehp);
            gsl_matrix_scale(&Eht.matrix, (double) DFEhp);
            decomposeEtoDS(NDV,NCC,Ehp,Dhp,SigmaShp);
            logPostPdfDSigmaP = logPostPdfDSigma(Dhp,SigmaShp,Ehp,Vinv,Sh,NDV,NCC,nmembers[h],eta);
            logPostPdfDSigmaT = logPostPdfDSigma(&Dht.matrix,&SigmaSht.matrix,&Eht.matrix,Vinv,Sh,NDV,NCC,nmembers[h],eta);
            logTransR = logtrnsR(&Eht.matrix,Ehp,NDV,NCC,DFEhp);
            logAcp = logPostPdfDSigmaP - logPostPdfDSigmaT + logTransR;
            Acp = exp(logAcp);
            if (Acp > 1.0) Acp = 1.0;
            unifRV = gsl_ran_flat(r,0.0,1.0);
            if (Acp > unifRV){
                if (nmembers[h] > 0) acceptSigma += 1/((double)(batchL*active));
                for (k = 0; k < totran; k++){
                    for (i = 0; i < totran; i++){
                        Sigmah[h][k*totran+i] = gsl_matrix_get(SigmaShp,k,i);
                        Eh[h][k*totran+i] = gsl_matrix_get(Ehp,k,i);
                        Dh[h][k*totran+i] = gsl_matrix_get(Dhp,k,i);
                    }
                }
            }//else they remain as from previous iteration
            if (Acp > unifRV)
                gsl_matrix_memcpy(SigmaSH,SigmaShp);
            else
                gsl_matrix_memcpy(SigmaSH,&SigmaSht.matrix);//current covariance matrix - for next step

            // - 4 - mu_h (assumes Sigma_{11} = 1 i.e. one count or one binomial response and continuous or binary covariates)           
            //puts("4");
            nu = gsl_matrix_submatrix(SigmaSH,1,0,p,1); //Sigma_{21} 
            for (k = 0; k < p; k++)
                nuh[h][k] = gsl_matrix_get(&nu.matrix,k,0); 
            ThCopy = gsl_matrix_submatrix(SigmaSH,1,1,p,p); //Sigma_{22}
            gsl_matrix_memcpy(Th,&ThCopy.matrix);
            gsl_blas_dgemm(CblasNoTrans,CblasTrans,-1.0,&nu.matrix,&nu.matrix,1.0,Th);
            ginv(p,tol,Th); // W_h^{-1}            
            if (nmembers[h] > 0){
                gsl_matrix_memcpy(nTDi,Th); 
                gsl_matrix_scale(nTDi,nmembers[h]);
                gsl_matrix_add(nTDi,SigmaMuInv); 
                gHalfInv(p,tol,nTDi); // (n_h W_h^-1 + D^-1)^{-1/2}
            }
            else
                gsl_matrix_memcpy(nTDi,SigmaMuHalf); //D^{1/2}                 
            for (k = 0; k < p; k++) // N(0,1)
	            gsl_vector_set(Z,k,gsl_ran_gaussian(r,1));
            gsl_blas_dgemv(CblasNoTrans,1.0,nTDi,Z,0.0,U1); // sample from MV normal; next add mean
            if (nmembers[h] > 0){
                SetSampleTotMuYX(p,NDV,n,h,ncomp,toths,compAlloc,X,nuh,latentYX);
                tothv = gsl_vector_view_array(toths,p);
                gsl_blas_dgemv(CblasNoTrans,1.0,Th,&tothv.vector,0.0,U2);
                gsl_vector_add(U2,SIMmu);
                gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,nTDi,nTDi,0.0,Dg1);                                 
                gsl_blas_dgemv(CblasNoTrans,1.0,Dg1,U2,0.0,U3);
            }
            else
                gsl_vector_memcpy(U3,MuMu);
            gsl_vector_add(U1,U3);
            for (k = 0; k < p; k++)
                muh[h][k] = gsl_vector_get(U1,k);

            // - 5 - Xi_h
            //puts("5");   
            if ((sw % batchL)==0 && WB > 0 && h==0){  	             
	            for (k = 0; k < nRespPars; k++){
	                if (acceptXi[k] > 0.25 && prec[k] > XLL) 	                 
	                    prec[k] -= MIN(0.01,1/sqrt(WB)) * prec[k]; 
	                if (acceptXi[k] <= 0.25 && prec[k] < XUL) 
	                    prec[k] += MIN(0.01,1/sqrt(WB)) * prec[k];
				}
		        for (k = 0; k < nRespPars; k++) 
		            acceptXi[k] = 0.0;   
	        }	                   
	        for (j = 0; j < XiLoop; j++){
                for (k = 0; k < nRespPars; k++)
                    XiC[k] = Xi[h][k];
                s = gsl_ran_flat(r,1.0,100000);                  
                propose(s,XiC,XiP,nRespPars,j,prec,family);           
                //puts("after"); 
                //if (family == 5) XiP[XiLoop-j-1] = XiC[XiLoop-j-1];
                temp3 = 0;
                for (k = 0; k < nRespPars; k++)
                    if (XiP[k] > 0.00001) temp3 +=1;                                      
                Dg2 = gsl_matrix_submatrix(SigmaSH,1,1,p,p);
                gsl_matrix_memcpy(Dg3,&Dg2.matrix);
                ginv(p,tol,Dg3); //Sigma_{xx}^{-1}
                nu2 = gsl_matrix_subcolumn(SigmaSH,0,1,p);
                gsl_blas_dgemv(CblasNoTrans,1.0,Dg3,&nu2.vector,0.0,nuSI);
                gsl_blas_ddot(nuSI,&nu2.vector,&nuSInu);
                if (nuSInu >= 1.0) nuSInu = 0.9999999999;
                storenuSInu[h] = nuSInu;
                for (k = 0; k < p; k++)
                    storenuSI[h][k] = gsl_vector_get(nuSI,k);
                //MH algorithm            
                XiLC = 0.0;  
                XiLP = 0.0;
                if (temp3 == nRespPars){
                    for (i = 0; i < n; i++){
                        if (compAlloc[i]==h && exp(XiLP) > 0){                        
                            nuSIxmm = 0.0;
                            for (k = 0; k < (NDV-1); k++)
                                nuSIxmm += storenuSI[h][k]*(latentYX[i][k+1]-muh[h][k]);                        
                            for (k = (NDV-1); k < p; k++)
                                nuSIxmm += storenuSI[h][k]*(X[k*n+i]-muh[h][k]);                                                   
                            calcGLMLimits(Y[i],H[i],XiC,&cymoC,&cyC,family);
                            calcGLMLimits(Y[i],H[i],XiP,&cymoP,&cyP,family);                            
                            XiLC += log(gsl_cdf_ugaussian_P((cyC-nuSIxmm)/sqrt(1-nuSInu)) - gsl_cdf_ugaussian_P((cymoC-nuSIxmm)/sqrt(1-nuSInu)));
                            XiLP += log(gsl_cdf_ugaussian_P((cyP-nuSIxmm)/sqrt(1-nuSInu)) - gsl_cdf_ugaussian_P((cymoP-nuSIxmm)/sqrt(1-nuSInu)));                       
                        }
                    }
                    if (family == 1){
                        NAcp = exp(XiLP)*gsl_ran_gamma_pdf(XiP[0],alphaXi[0],1/betaXi[0])*
                               gsl_ran_gamma_pdf(XiC[0],prec[0]*XiP[0]*XiP[0],1/(prec[0]*XiP[0]));
                        DAcp = exp(XiLC)*gsl_ran_gamma_pdf(XiC[0],alphaXi[0],1/betaXi[0])*
                               gsl_ran_gamma_pdf(XiP[0],prec[0]*XiC[0]*XiC[0],1/(prec[0]*XiC[0]));
                    }
                    else if (family == 2){
                        temp2 = XiP[0] - 1 + XiP[0]*(1-XiP[0])*(1-XiP[0])*prec[0];
                        temp = temp2 * XiP[0]/(1-XiP[0]);
                        NAcp = exp(XiLP)*gsl_ran_beta_pdf(XiP[0],alphaXi[0],betaXi[0])*
                               gsl_ran_beta_pdf(XiC[0],temp,temp2);
                        temp2 = XiC[0] - 1 + XiC[0]*(1-XiC[0])*(1-XiC[0])*prec[0];
                        temp = temp2 * XiC[0]/(1-XiC[0]);
                        DAcp = exp(XiLC)*gsl_ran_beta_pdf(XiC[0],alphaXi[0],betaXi[0])*
                               gsl_ran_beta_pdf(XiP[0],temp,temp2);
                    }
                    else if (family == 3 || family == 4 || family == 6 || family == 7){
                        NAcp = exp(XiLP)*gsl_ran_gamma_pdf(XiP[0],alphaXi[0],1/betaXi[0])*
                                         gsl_ran_gamma_pdf(XiP[1],alphaXi[1],1/betaXi[1])*
                                         gsl_ran_gamma_pdf(XiC[0],prec[0]*XiP[0]*XiP[0],1/(prec[0]*XiP[0]))*
                                         gsl_ran_gamma_pdf(XiC[1],prec[1]*XiP[1]*XiP[1],1/(prec[1]*XiP[1]));
                        DAcp = exp(XiLC)*gsl_ran_gamma_pdf(XiC[0],alphaXi[0],1/betaXi[0])*
                                         gsl_ran_gamma_pdf(XiC[1],alphaXi[1],1/betaXi[1])*
                                         gsl_ran_gamma_pdf(XiP[0],prec[0]*XiC[0]*XiC[0],1/(prec[0]*XiC[0]))*
                                         gsl_ran_gamma_pdf(XiP[1],prec[1]*XiC[1]*XiC[1],1/(prec[1]*XiC[1]));
                    }
                    else if (family == 5){
                        temp = 1;
                        if (j==0) temp = gsl_ran_gamma_pdf(XiP[0],alphaXi[0],1/betaXi[0])*
                                         gsl_ran_gamma_pdf(XiC[0],prec[0]*XiP[0]*XiP[0],1/(prec[0]*XiP[0]));
                        if (j==1) temp = gsl_ran_gaussian_pdf(XiP[1]-alphaXi[1],betaXi[1]);//random walk: gsl_ran_gaussian_pdf(XiP[1]-XiC[1],prec[1]);
                        NAcp = exp(XiLP)*temp;
                        temp = 1;
                        if (j==0) temp = gsl_ran_gamma_pdf(XiC[0],alphaXi[0],1/betaXi[0])*
                                         gsl_ran_gamma_pdf(XiP[0],prec[0]*XiC[0]*XiC[0],1/(prec[0]*XiC[0]));
                        if (j==1) temp = gsl_ran_gaussian_pdf(XiC[1]-alphaXi[1],betaXi[1]);//gsl_ran_gaussian_pdf(XiP[1]-XiC[1],prec[1]);
                        DAcp = exp(XiLC)*temp;
                    }
                    else if (family == 8){
                        NAcp = exp(XiLP)*gsl_ran_gamma_pdf(XiP[0],alphaXi[0],1/betaXi[0])*
                                         gsl_ran_gamma_pdf(XiP[1],alphaXi[1],1/betaXi[1])*
                                         gsl_ran_gaussian_pdf(XiP[2]-alphaXi[2],betaXi[2])*
                                         gsl_ran_gamma_pdf(XiC[0],prec[0]*XiP[0]*XiP[0],1/(prec[0]*XiP[0]))*
                                         gsl_ran_gamma_pdf(XiC[1],prec[1]*XiP[1]*XiP[1],1/(prec[1]*XiP[1]));
                                         //gsl_ran_gaussian_pdf(XiP[2]-XiC[2],prec[3]);
                        DAcp = exp(XiLC)*gsl_ran_gamma_pdf(XiC[0],alphaXi[0],1/betaXi[0])*
                                         gsl_ran_gamma_pdf(XiC[1],alphaXi[1],1/betaXi[1])*
                                         gsl_ran_gaussian_pdf(XiC[2]-alphaXi[2],betaXi[2])*
                                         gsl_ran_gamma_pdf(XiP[0],prec[0]*XiC[0]*XiC[0],1/(prec[0]*XiC[0]))*
                                         gsl_ran_gamma_pdf(XiP[1],prec[1]*XiC[1]*XiC[1],1/(prec[1]*XiC[1]));
                                         //gsl_ran_gaussian_pdf(XiP[2]-XiC[2],prec[3]);
                    }
                    Acp = NAcp/DAcp; 
                    //Rprintf("%i %i %i %i %i %f %f %f %f %f %i %f %f %i\n",sw,h,nmembers[h],j,XiLoop,XiC[0],XiC[1],XiP[0],XiP[1],Acp,nRespPars,prec[0],prec[1],family);
                    if (Acp > 1.0) Acp = 1.0;
                    unifRV = gsl_ran_flat(r,0.0,1.0);  
                    if (Acp > unifRV){
                        if (family < 5){ 
                            if (nmembers[h] > 0) 
                                for (k = 0; k < nRespPars; k++) 
                                    acceptXi[k] += 1/((double)(batchL*active));
                        }
                        if (family == 5){
                            if (nmembers[h] > 0) 
                                acceptXi[j] += 1/((double)(batchL*active));
                        } 
                        for (k = 0; k < nRespPars; k++) 
                            Xi[h][k] = XiP[k];
                    }  
                    else   
                        for (k = 0; k < nRespPars; k++) 
                            Xi[h][k] = XiC[k];
                }
                else
                    for (k = 0; k < nRespPars; k++) 
                        Xi[h][k] = XiC[k];
            }                         
            for (k = 0; k < nRespPars; k++) 
                XiC[k] = Xi[h][k];

            // - 6 - Impute latent y^* 
            //puts("6");
            MNCondParams1of2(NDV,NCC,SigmaSH,tol,PartMean,CondCov,params); 
            for (i = 0; i < n; i++){
                if (compAlloc[i]==h){                     
                    calcGLMLimitsYX(Y,H,X,i,NDC,n,XiC,lower,upper,family);                     
                    
                    baseMUS[0] = 0.0;
                    for (k = 0; k < p; k++)
                        baseMUS[k+1] = muh[h][k];
                    MUSTAR = gsl_vector_view_array(baseMUS,totran);                 
                    for (k = 0; k < NCC; k++) 
                        baseC[k] = X[NDC*n+i];                    
                    Ci = gsl_vector_view_array(baseC,NCC);                                        
                    MNCondParams2of2(NDV,NCC,&MUSTAR.vector,&Ci.vector,PartMean,CondMean,params);                                        
                    for (k = 0; k < NDV; k++)
                        gsl_vector_set(yxstar,k,latentYX[i][k]);
                    s = gsl_ran_flat(r,1.0,100000);
                    sampleTMN(s,NDV,lower,upper,yxstar,CondMean,CondCov,tol);
                    for (k = 0; k < NDV; k++)
                        latentYX[i][k] = gsl_vector_get(yxstar,k);                    
                } 
            }
            // Moving on to next cluster
            sumPIh += pi[h];  
            //if (sw ==337) Rprintf("%i %i %i %f %f \n",sw,nUpdated,h,sumPIh,minU);          
            h++;
        }// h loop
        nUpdated = h;        
        
        // - check -
        for (i = 0; i < n; i++){
            if (compAlloc[i] >= ncomp){
                for (k = 0; k < NDV; k++)
                    latentYX[i][k] = 0.0;
                }
        }
 
        // - 7 - Update the u_i's from unif(0,pi_k_i)
        //puts("7");
        if (sampler == 1){ 
            for (i = 0; i < n; i++){
                if (compAlloc[i] < ncomp) u[i] = gsl_ran_flat(r,0.0,pi[compAlloc[i]]);
                else u[i] = 0.0;                
            }
            minU = 1.0;
            for (i = 0; i < n; i++)
                if (u[i] < minU) minU = u[i];
        }  
 
        // - 8 - P(delta_j = k_j) // needs fixing
        //Rprintf("%s %i \n","pdelta",sw); 
        for (h = 0; h < nUpdated; h++){			
            for (k = 0; k < nRespPars; k++) XiC[k] = Xi[h][k];                            
            for (k = 0; k < totran*totran; k++) baseSigmaSh[k] = Sigmah[h][k];                        
            SigmaSht = gsl_matrix_view_array(baseSigmaSh,totran,totran);                        
            SigmaSht2 = gsl_matrix_submatrix(&SigmaSht.matrix,NDV,NDV,NCC,NCC);                                                           
            MNCondParams1of2b(NDV,NCC,&SigmaSht.matrix,tol,PartMean,&SigmaSht2.matrix,params);//inverses Sht2                                        
            for (i = 0; i < n; i++){                
                if (pi[h] > u[i]){                                                                                                  
                    baseMUS[0] = 0.0;  
                    for (k = 0; k < p; k++)
                        baseMUS[k+1] = muh[h][k];
                    MUSTAR = gsl_vector_view_array(baseMUS,totran);                                     
                    for (k = 0; k < NCC; k++)
                        baseC[k] = X[(NDC+k)*n+i];                    
                    Ci = gsl_vector_view_array(baseC,NCC);                                                                               
                    MNCondParams2of2(NDV,NCC,&MUSTAR.vector,&Ci.vector,PartMean,CondMean,params);                       
                    calcGLMLimitsYX(Y,H,X,i,NDC,n,XiC,lower,upper,family);                                                            
                    nocuba = 0;
                    for (k = 0; k < NDV; k++)
                        if (((lower[k] - params[k]) > tol2 / params[k+NDV]) || ((params[k] - upper[k]) > tol2 / params[k+NDV]))
                            nocuba += 1;                    
                    if (nocuba == 0){                                                               
                        for (k = 0; k < NDV; k++){
                            if (lower[k]==-lmt) lower[k] = params[k] - tol2 / params[NDV+k];
                            if (upper[k]==lmt) upper[k] = params[k] + tol2 / params[NDV+k];
                        }                    
                        if (NDV > 3) p_hcubature(1,MultiNormalPDF,params,NDV,lower,upper,0,0.0001,0,ERROR_INDIVIDUAL,&val,&err);
                        if (NDV <= 3) p_pcubature(1,MultiNormalPDF,params,NDV,lower,upper,0,0.0001,0,ERROR_INDIVIDUAL,&val,&err);                                
                    }else
                        val = 0.0; 
                    temp2 = 0;
                    if (NCC > 0){
                        for (k = 0; k < NCC; k++) 
                            baseXmM[k] = X[NDC*n+k*n+i]-muh[h][NDC+k];                                                            
                        XmMview = gsl_vector_view_array(baseXmM,NCC);                    
                        temp2 = logMVNormalpdf3(NCC,&XmMview.vector,vecZero2,&SigmaSht2.matrix);
				    }				    
                    if (sampler == 1) pdih[h][i] = exp(log(val)+temp2);
                    else if (sampler == 2) pdih[h][i] = exp(log(val)+temp2+log(pi[h]));
                    if (gsl_isnan(pdih[h][i])) {pdih[h][i] = 0.0;Rprintf("%s %i %i %i %f %f %i %f %f\n", "oops:",sw,h,i,XiC[0],XiC[1],Y[i],H[i],
						cdf_generalized_poisson_P3(18,4.7,0.006,&temp));}
                }
                else pdih[h][i] = 0.0;               		         
                if (sw==1 && h==3 && i==263 && 0==1) Rprintf(
                "%s %i %i %i %i %f | %f %f | %f %f %f | %f %f %f | %f %f %f %f | %f %f | %f %f %f %f %f \n",
                "prob",sw,h,i,Y[i],H[i], pi[h],u[i], val,temp2,log(pdih[h][i]),
                XiC[0],XiC[1],XiC[2],
                lower[0],upper[0],lower[1],upper[1],latentYX[i][0],latentYX[i][1],
                params[0],params[1],params[2],params[3],params[4]);                                                                
            }
        }
        for (h = nUpdated; h < ncomp; h++)
            for (i = 0; i < n; i++)
                pdih[h][i] = 0.0;

        // - 9 - Update allocation to
        //Rprintf("%s %i \n","alloc",sw);
        s = gsl_ran_flat(r,1.0,100000);
        allocation(s,n,ncomp,pdih,compAlloc,sw);

        // - 10 - Find nmembers
        //puts("10");
        findNmembers(n,ncomp,compAlloc,nmembers);
        active = 0;
        for (h = 0; h < ncomp; h++) if (nmembers[h] > 0) active += 1;
        
        // - 12a - Label switching
        //puts("11");
        s = gsl_ran_flat(r,1.0,100000);
        labelSwtchANEW(s,n,totran,p,ncomp,nRespPars,Sigmah,Dh,Eh,nuh,muh,Xi,storenuSInu,storenuSI,nmembers,compAlloc,pi);

        // - 12b - Label switching
        s = gsl_ran_flat(r,1.0,100000);
        labelSwtchBNEW(s,n,totran,p,ncomp,nRespPars,Sigmah,Dh,Eh,nuh,muh,Xi,storenuSInu,storenuSI,nmembers,compAlloc,Vdp);

        // - 13 - Update concentration parameter
        s = gsl_ran_flat(r,1.0,100000);
        newalpha = updateAlpha(s,n,ncomp,Alphaa,Alphab,TruncAlpha,nmembers,alpha);
        alpha = newalpha;

        // - 14 - Predictions
        if (((sw - burn) >= 0) && (((sw - burn ) % thin) == 0) && (npred > 0)){
            //puts("denominator");
            for (i = 0; i < npred; i++) denomPred[i] = 0.0;                                   
            for (i = 0; i < npred; i++){
                for (h = 0; h < ncomp; h++){                                                        
                    move = 0;
                    NDCPi = 0;
                    NCCPi = 0;
                    for (k = 1; k < totran; k++){
                        if (c[k*npred+i]==1){
                            if (k < NDV) NDCPi +=1;
                            else NCCPi += 1;
                            for (l = 1; l < totran; l++)
                                if (c[l*npred+i]==1) 
                                    baseSigmaSh[move++] = Sigmah[h][k*totran+l];                                                
						}                        				
				    }				    
                    SigmaSht = gsl_matrix_view_array(baseSigmaSh,NDCPi+NCCPi,NDCPi+NCCPi);                                                                                
                    SigmaSht3 = gsl_matrix_submatrix(&SigmaSht.matrix,NDCPi,NDCPi,NCCPi,NCCPi);                                                          
                    val = 1;
                    if (NDCPi > 0){                        
                        move = 0;
                        for (k = 1; k < totran; k++)
                            if (c[k*npred+i]==1) 
                                baseMUS[move++] = muh[h][k-1];
                        MUSTAR = gsl_vector_view_array(baseMUS,move);  
                        move = 0;
                        for (k = 0; k < NCC; k++)
                            if (c[(k+NDV)*npred+i]==1)
                                baseC[move++] = Xpred[(k+NDC)*npred+i];
                        Ci = gsl_vector_view_array(baseC,NCCPi);
                        MNCondParams2(NDCPi,NCCPi,&MUSTAR.vector,&Ci.vector,&SigmaSht.matrix,tol,params);                                                                          
                        move = 0;
                        for (k = 0; k < NDC; k++){
                            if (Xpred[k*npred+i]==0 && c[(k+1)*npred+i]==1) {lowerX[move] = -lmt; upperX[move] = 0; move++;}
                            if (Xpred[k*npred+i]==1 && c[(k+1)*npred+i]==1) {lowerX[move] = 0; upperX[move] = lmt; move++;}
	                    }                                        
                        nocuba = 0;
                        for (k = 0; k < NDCPi; k++)
                            if (((lowerX[k] - params[k]) > tol2 / params[k+NDCPi]) || ((params[k] - upperX[k]) > tol2 / params[k+NDCPi]))
                                nocuba += 1;                    
                        if (nocuba == 0){                                                               
                            for (k = 0; k < NDCPi; k++){
                                if (lowerX[k]==-lmt) lowerX[k] = params[k] - tol2 / params[NDCPi+k];
                                if (upperX[k]==lmt) upperX[k] = params[k] + tol2 / params[NDCPi+k];
                            }                    
                            if (NDCPi > 3) p_hcubature(1,MultiNormalPDF,params,NDCPi,lowerX,upperX,0,0.0001,0,ERROR_INDIVIDUAL,&val,&err);
                            if (NDCPi <= 3) p_pcubature(1,MultiNormalPDF,params,NDCPi,lowerX,upperX,0,0.0001,0,ERROR_INDIVIDUAL,&val,&err);                                
                        }else
                            val = 0.0;
					}					
					CPDF[i][h] = 1;
					if (NCCPi > 0){
					    ginv(NCCPi,tol,&SigmaSht3.matrix);                                                                                                                                                  
                        move = 0;
                        for (k = 0; k < NCC; k++)
                            if (c[(k+NDV)*npred+i]==1)
                                baseXmM2[move++] = Xpred[(NDC+k)*npred+i]-muh[h][NDC+k];                                        
                        XmMview = gsl_vector_view_array(baseXmM2,NCCPi);                    
                        gsl_vector_view ZeroSub = gsl_vector_subvector(vecZero2,0,NCCPi);                    
                        temp2 = logMVNormalpdf3(NCCPi,&XmMview.vector,&ZeroSub.vector,&SigmaSht3.matrix);                    
                        CPDF[i][h] = exp(temp2); 
					}                                                                                                
                    denomPred[i] += pi[h]*CPDF[i][h]*val;
                }
            }                                                            
            //puts("numerator");  
            Iend = 0;
            for (k = 0; k < ncomp; k++) for (l = 0; l < maxy; l++) StoreUpper[k][l] = -123456789.0;
			for (k = 0; k < maxy; k++) for (i = 0; i < npred; i++) StorePD[i][k] = 0.0;
            for (k = 0; k < ncomp; k++) saveCdf[k] = 0.0;
            for (i = 0; i < npred; i++){
                move = 0;
                for (k = 0; k < NCC; k++) 
                    if (c[(k+NDV)*npred+i]==1)
                        baseC[move++] = Xpred[(k+NDC)*npred+i]; 
                NCCPi = move;
                Ci = gsl_vector_view_array(baseC,NCCPi);        
                start = 0; cdfy[i] = 0;                              
                while(cdfy[i] < (PredTrunc+0) && start < maxy){					                   
                    for (h = 0; h < ncomp; h++){                                                                                                        
                        if (StoreUpper[h][start]==-123456789.0){
							for (k = 0; k < nRespPars; k++) XiC[k] = Xi[h][k];                            
                            if (start == 0)   
                                Lower = -lmt; 
                            else 
                                Lower = StoreUpper[h][start-1];                                                           
                            if (family == 1) saveCdf[h] += gsl_ran_poisson_pdf(start,Hpred[0]*XiC[0]);
                            else if (family == 2) saveCdf[h] += gsl_ran_binomial_pdf(start,XiC[0],Hpred[0]);
                            else if (family == 3) saveCdf[h] += gsl_ran_negative_binomial_pdf(start,XiC[1]/(Hpred[0]+XiC[1]),XiC[0]);
                            else if (family == 4) saveCdf[h] += exp(gsl_sf_lnchoose(Hpred[0],start) + gsl_sf_lnbeta(start+XiC[0],Hpred[0]-start+XiC[1]) - gsl_sf_lnbeta(XiC[0],XiC[1]));
                            else if (family == 5){		
		                        mu = Hpred[0]*XiC[0]; f = XiC[1];//sqrt(Hpred[0]);
		                        if (start == 0) saveCdf[h] = cdf_generalized_poisson_P3(0,mu,f,&normConst);		                        
		                        else if (start > 0){
		                            B = -mu/(f-1);
		                            if ((f>=1)||(f < 1 && start < B))
		                            saveCdf[h] += exp(log(mu)+(start-1)*log(mu+(f-1)*start)-start*log(f)-(mu+(f-1)*start)/f-gsl_sf_lnfact(start))/normConst;		
							    } 
	                        }	                        
	                        if (saveCdf[h] > 1.0) saveCdf[h] = 1.0;                    
                            Upper = gsl_cdf_ugaussian_Pinv(saveCdf[h]);
                            if (Upper < -lmt) Upper = -lmt;
                            if (Upper > lmt) Upper = lmt;                                                                                                                                                                                                                             
                            StoreUpper[h][start] = Upper;                             
                            if (family==5 && (gsl_isinf(normConst) || normConst==0.0)) {Iend=1; Rprintf("%s %i %i %i %i \n", "ooPs :",sw,h,i,seed);}                                                                                            
                        }
                        else{
                            if (start == 0) Lower = -lmt;
                            else Lower = StoreUpper[h][start-1];
                            Upper = StoreUpper[h][start];
                        }
                        lower[0] = Lower;
                        upper[0] = Upper;                        
                        move = 0;
                        for (k = 0; k < NDC; k++){
                            if (Xpred[k*npred+i]==0 && c[(k+1)*npred+i]==1) {lower[move+1] = -lmt; upper[move+1] = 0; move++;}
                            if (Xpred[k*npred+i]==1 && c[(k+1)*npred+i]==1) {lower[move+1] = 0; upper[move+1] = lmt; move++;}
	                    }	                     	                    
	                    if (start == 0){
	                        move = 0;
                            NDVPi = 0;  
                            for (k = 0; k < totran; k++){								
                                if (c[k*npred+i]==1){
                                    if (k < NDV) NDVPi +=1;
                                    for (l = 0; l < totran; l++)
                                        if (c[l*npred+i]==1) 
                                            baseSigmaSh[move++] = Sigmah[h][k*totran+l];                                                
						        }                        				
				            }
                            SigmaSht = gsl_matrix_view_array(baseSigmaSh,NDVPi+NCCPi,NDVPi+NCCPi);                                                                                          	                                                                               
                            baseMUS[0] = 0.0;
                            move = 1;
                            for (k = 1; k < totran; k++)
                                if (c[k*npred+i]==1) 
                                    baseMUS[move++] = muh[h][k-1];
                            MUSTAR = gsl_vector_view_array(baseMUS,move);                                                                                      
                            MNCondParams2(NDVPi,NCCPi,&MUSTAR.vector,&Ci.vector,&SigmaSht.matrix,tol,params);                                                                                                                                                                                                                                                                                                          
                            for (k = 0; k < dimP; k++) 
                                StoreParams[h][k] = params[k];   
                        }
                        else
                            for (k = 0; k < dimP; k++) 
                                params[k] = StoreParams[h][k];                                                                  
                        nocuba = 0; 
                        for (k = 0; k < NDVPi; k++)
                            if (((lower[k] - params[k]) > tol2 / params[k+NDVPi]) || ((params[k] - upper[k]) > tol2 / params[k+NDVPi]))
                                nocuba += 1;                    
                        if (nocuba == 0){                                                                
                            for (k = 0; k < NDVPi; k++){
                                if (lower[k]==-lmt) lower[k] = params[k] - tol2 / params[NDVPi+k];
                                if (upper[k]==lmt) upper[k] = params[k] + tol2 / params[NDVPi+k];
                            }                    
                            if (NDVPi > 3) p_hcubature(1,MultiNormalPDF,params,NDVPi,lower,upper,0,0.0001,0,ERROR_INDIVIDUAL,&val,&err);
                            if (NDVPi <= 3) p_pcubature(1,MultiNormalPDF,params,NDVPi,lower,upper,0,0.0001,0,ERROR_INDIVIDUAL,&val,&err);                                 
                            if (sw==3800 && i==4 && start==5 && h==20 && 1==0) Rprintf("%s %f %f %f %f %f %f %f %f %f \n","val:",
                            val,StoreUpper[h][0],StoreUpper[h][1],StoreUpper[h][2],StoreUpper[h][3],StoreUpper[h][4],StoreUpper[h][5],lower[0],upper[0]);                        
                        }else 
                            val = 0.0;                                                                                       
                        StorePD[i][start] += pi[h] * CPDF[i][h] * val / denomPred[i];
                        if (gsl_isnan(StorePD[i][start])) {StorePD[i][start] = 0.0;Rprintf("%s %i %i %i %i \n", "ooPDs:",sw,h,i,start);}
                        if (val>1 && 1==0) {Rprintf("%s %i %i %i %i %f %f %f %f %f \n", "val:",sw,h,i,start,lower[0],upper[0],params[0],params[1],val);}
                        if (sw==96 && 1==0) fprintf(out_file4, "%i %i %i %i %f %f %f %f %f %f %f %f %f %f %f \n",
                        sw,i,start,h,cdfy[i],pi[h],CPDF[i][h],params[0],params[1],params[2],lower[0],upper[0],val,StorePD[i][start],denomPred[i]);                        
					}
					cdfy[i] += StorePD[i][start];
					start++;
                }
                //Rprintf("%i %i %i %f \n",sw,i,start,cdfy);                                       
            }            
            // posterior summaries 
            if (Iend == 0){ 
                //Mean
                for (i = 0; i < npred; i++) TMR[i] = 0.0;
                for (i = 0; i < npred; i++)
                    for (k = 0; k < maxy; k++)
                        TMR[i] += (double) (StorePD[i][k] / cdfy[i] * k); //mean regression 
                for (i = 0; i < npred; i++)
                    meanReg[i] += (double) TMR[i]/(nSamples); //mean of mean regression
                //Q1,Q2,Q3                
                for (i = 0; i < npred; i++){
                    k=0;
                    sumProb = StorePD[i][k];
                    while (sumProb <= 0.05){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ05R[i] = k;                    
                    while (sumProb <= 0.10){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ10R[i] = k;                    
                    while (sumProb <= 0.15){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ15R[i] = k;                    
                    while (sumProb <= 0.20){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ20R[i] = k;                    
                    while (sumProb <= 0.25){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ25R[i] = k;
                    q1Reg[i] += (double) k/(nSamples); // mean q1 reg
                    while (sumProb <= 0.50){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TMdR[i] = k; //median regression
                    medianReg[i] += (double) k/(nSamples); //mean q2 reg                    
                    while (sumProb <= 0.75){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ75R[i] = k;                    
                    q3Reg[i] += (double) k/(nSamples); //mean q3 reg
                    while (sumProb <= 0.80){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ80R[i] = k;                  
                    while (sumProb <= 0.85){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ85R[i] = k;                    
                    while (sumProb <= 0.90){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ90R[i] = k;                  
                    while (sumProb <= 0.95){
                        k++;
                        sumProb += StorePD[i][k] / cdfy[i];
                    }
                    TQ95R[i] = k;                     
                }                               
                // Mode
                for (i = 0; i < npred; i++){
                    mode = gsl_stats_max_index(StorePD[i],1,maxy);
                    modeReg[i] += (double) (mode)/nSamples; 
                }
                //Mean and Var of the probabilities
                for (i = 0; i < npred; i++){
                    for (k = 0; k < maxy; k++){
                        denReg[k+maxy*i] += (double) (StorePD[i][k]/cdfy[i])/nSamples;
                        denVar[k+maxy*i] += (double) (StorePD[i][k]/cdfy[i])*(StorePD[i][k]/cdfy[i])/nSamples;                        
                    }
                }
            }
        }//end pred
        // - 15 - Write to files
        if (((sw - burn) >= 0) && (((sw - burn ) % thin) == 0) && (WF == 1) && (Iend==0)){

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < (totran*totran); j++)
                    fprintf(out_file1, "%f ", Sigmah[h][j]);
                fprintf (out_file1, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < p; j++)
                    fprintf(out_file2, "%f ", muh[h][j]);
                fprintf (out_file2, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < nRespPars; j++)
                    fprintf(out_file3, "%f ", Xi[h][j]);
                fprintf (out_file3, "\n");
            }

            fprintf(out_file4, "%f \n", alpha); 

            for (i = 0; i < n; i++)
                fprintf(out_file5, "%i ", compAlloc[i]);
            fprintf (out_file5, "\n");

            for (h = 0; h < ncomp; h++)
                fprintf(out_file6, "%i ", nmembers[h]);
            fprintf (out_file6, "\n");

            fprintf(out_file7, "%i \n", nUpdated);

            for (i = 0; i < npred; i++)
                fprintf(out_file8, "%f ", TMR[i]);
            fprintf (out_file8, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file9, "%f ", TMdR[i]);
            fprintf (out_file9, "\n");            

            for (i = 0; i < npred; i++)
                fprintf(out_file10, "%f ", TQ05R[i]);
            fprintf(out_file10, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file11, "%f ", TQ10R[i]);
            fprintf(out_file11, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file12, "%f ", TQ15R[i]);
            fprintf(out_file12, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file13, "%f ", TQ20R[i]);
            fprintf(out_file13, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file14, "%f ", TQ25R[i]);
            fprintf(out_file14, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file15, "%f ", TQ75R[i]);
            fprintf(out_file15, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file16, "%f ", TQ80R[i]);
            fprintf(out_file16, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file17, "%f ", TQ85R[i]);
            fprintf(out_file17, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file18, "%f ", TQ90R[i]);
            fprintf(out_file18, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file19, "%f ", TQ95R[i]);
            fprintf(out_file19, "\n");

        }
        if (Iend == 1) sw -= 1;
        if ((sw==(sweeps-1)) && (!((sw+1) % 500)==0)) Rprintf("%i %s \n",sw+1, "posterior samples...");
    }//end of sw
    //Free up random number
    gsl_rng_free (r);

    //Close files
    fclose(out_file1); fclose(out_file2); fclose(out_file3); fclose(out_file4); fclose(out_file5);
    fclose(out_file6); fclose(out_file7); fclose(out_file8); fclose(out_file9); fclose(out_file10); 
    fclose(out_file11); fclose(out_file12); fclose(out_file13); fclose(out_file14); fclose(out_file15); 
    fclose(out_file16); fclose(out_file17); fclose(out_file18); fclose(out_file19);

    //Free up gsl matrices
    gsl_matrix_free(V); gsl_matrix_free(Vinv);
    gsl_vector_free(MuMu); gsl_matrix_free(SigmaMuInv); gsl_matrix_free(SigmaMuHalf); gsl_vector_free(SIMmu);
    gsl_matrix_free(Sh);
    gsl_vector_free(nuSI);
    gsl_matrix_free(Dg1); gsl_matrix_free(Dg3);
    gsl_matrix_free(nTDi);
    gsl_vector_free(vecZero); gsl_vector_free(vecZero2);
    gsl_vector_free(Z); gsl_vector_free(U1); gsl_vector_free(U2); gsl_vector_free(U3);
    gsl_vector_free(CondMean); gsl_vector_free(yxstar);
    gsl_matrix_free(PartMean); gsl_matrix_free(CondCov); 
    gsl_matrix_free(Ehp); gsl_matrix_free(Dhp); gsl_matrix_free(SigmaShp);
    gsl_matrix_free(SigmaSH); gsl_matrix_free(Th);
    //gsl_vector_free(CondMean2); gsl_matrix_free(PartMean2); gsl_matrix_free(CondCov2);
}
