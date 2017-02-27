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
//include "cubature.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void OneResLtnt(int *seed1, double *X, int *Y, double *H,
                int *sweeps1, int *burn1, int *thin1, int *ncomp1, int *n1, int *p1,
                double *Vprior, double *Vdf1,
                double *munu, double *sigmanu,
                double *mumu, double *sigmamu,
                double *alphaXi, double *betaXi,
                double *Alphaa1, double *Alphab1, double *TruncAlpha1,
                double *xbar, double *xsd, double *Ymean, double *prec1,
                int *family1, int *sampler1,
                int *npred1, double *Xpred, double *Hpred, int *allEqlI, int *maxy1,
                double *meanReg, double *modeReg,
                double *Q05Reg, double *Q10Reg, double *Q15Reg, double *Q20Reg, double *Q25Reg,
                double *Q50Reg, double *Q75Reg, double *Q80Reg, double *Q85Reg, double *Q90Reg,
                double *Q95Reg, double *denReg, double *denVar,
                char **WorkingDir, int *WF1, int *compA)
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
    int family = family1[0]; //1=poisson, 2=binomial, 3=negative binomial, 4=beta binomial, 5=generalized poisson, 6=com-poisson, 7=hyper-poisson, 8=ctpd
    int sampler = sampler1[0]; //1=slice, 2=truncated

    //Declare variables for loops: h for components, k,j for covariates, i for subjects, sw for sweeps
    int h, i, j, k, sw;

    // Dimensions
    int n = n1[0]; // number of observations
    int p = p1[0]; // number of covariates
    int ncomp = ncomp1[0]; //number of clusters/components
    int npred = npred1[0]; //number of predictions
    int nRespPars = 1; // number of parameters in response pmf
    if (family > 2) nRespPars = 2;
    if (family == 8) nRespPars = 3;

    //Tolerance level
    double tol = 0.0001; //tolerance level for eigen values when inverting matrices

    // Prior parameters
    // - 1 - Th: V, V^{-1}, df
    gsl_matrix *V = gsl_matrix_alloc(p,p);
	gsl_matrix *Vinv = gsl_matrix_alloc(p,p);
	for (j = 0; j < p; j++)
        for (k = 0; k < p; k++)
            gsl_matrix_set(V,j,k,Vprior[j*p+k]);
    gsl_matrix_memcpy(Vinv,V);
    Inverse(p,Vinv);
    double Vdf = Vdf1[0];
    // - 2 - Nuh: mu_nu and Sigma_nu
    gsl_vector *MuNu = gsl_vector_alloc(p);
    gsl_matrix *SigmaNuInv = gsl_matrix_alloc(p,p);
    gsl_matrix *SigmaNuHalf = gsl_matrix_alloc(p,p);
    gsl_vector *SIMnu = gsl_vector_calloc(p);
    for (j = 0; j < p; j++)
        gsl_vector_set(MuNu,j,munu[j]);
    for (j = 0; j < p; j++)
        for (k = 0; k < p; k++)
            gsl_matrix_set(SigmaNuInv,j,k,sigmanu[j*p+k]);
    Inverse(p,SigmaNuInv);
    gsl_matrix_memcpy(SigmaNuHalf,SigmaNuInv);
    gHalfInv(p,tol,SigmaNuHalf);
    gsl_blas_dgemv(CblasNoTrans,1.0,SigmaNuInv,MuNu,0.0,SIMnu); // Sigma_nu^{-1} Mu_nu
    // - 3 - Muh: mean  and Sigma_nu
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
    // - 4 - A. gamma ~ Gamma(alpha,beta) for poisson. B. gamma ~ Beta(alpha,beta) for binomial
    // Directly use arguments of the function
    // - 5 - concentration parameter alpha ~ Gamma(a,b) I[alpha>c]
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
         *out_file18, *out_file19, *out_file20, *out_file21, *out_file22;

    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Th.txt");
    out_file1 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Sigmah.txt");
    out_file2 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.SigmahI.txt");
    out_file3 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.nuh.txt");
    out_file4 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.muh.txt");
    out_file5 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.xih.txt");
    out_file6 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.alpha.txt");
    out_file7 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.compAlloc.txt");
    out_file8 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.nmembers.txt");
    out_file9 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Updated.txt");
    out_file10 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.MeanReg.txt");
    out_file11 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q05Reg.txt");
    out_file12 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q10Reg.txt");
    out_file13 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q15Reg.txt");
    out_file14 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q20Reg.txt");
    out_file15 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q25Reg.txt");
    out_file16 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q50Reg.txt");
    out_file17 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q75Reg.txt");
    out_file18 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q80Reg.txt");
    out_file19 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q85Reg.txt");
    out_file20 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q90Reg.txt");
    out_file21 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Q95Reg.txt");
    out_file22 = fopen(path_name, "wt");

    // Define quantities of interest
    double Th1[ncomp][p*p]; // to store Th
    double Sigmah[ncomp][p*p]; // to store \Sigma_{h}
    double SigmahI[ncomp][p*p]; // to store \Sigma_{h}^{-1}
    double nuh[ncomp][p]; // to store nu
    double muh[ncomp][p]; // mu_h
    double Xi[ncomp][nRespPars]; // response pmf parameters
    double alpha; //concentration parameter
    int compAlloc[n]; // allocation
    int nmembers[ncomp]; // number of cluster members
    double pdjkj[ncomp][n]; // matrix to store P(delta_j = k_j)

    // Define other quantities
    double pi[ncomp]; //stick breaking probabilities
    double Vdp[ncomp]; //V[i] ~ Beta(1,alpha)
    double newalpha; //store concentration parameter
    int nmb;
    double u[n]; // for implementing a slice sampler - Walker 2006
    double minU, sumPIh;
    int nUpdated, upperL;
    double latentx[n]; // continuous latext variables underlying discrete ones
    double XiC[nRespPars];
    double XiP[nRespPars];
    double XiLC, XiLP, cyC, cyP, cymoC, cymoP, Acp, NAcp, DAcp; //MH  step
    double nuSInu, nuSIxmm; //nu^T S^-1 nu and nu^T SI (x-mu)
    double storenuSInu[ncomp]; // to store nu^T Sigma^(-1) nu
    double storenuSI[ncomp][p]; // to store nu^T Sigma^(-1)
    double totprodhs[p];  // total of clusters: (x-mu)*x1
    double sumxsq; // sum of cluster latentx[i]^2
    double toths[p];  // total of clusters: x-mu-y1^* * nu
    double temp, temp2;
    int temp3;
    double lower, upper, lower2, upper2;
    double prec[nRespPars];
    for (j = 0; j < nRespPars; j++) prec[j] = prec1[j];

    // Declare and allocate gsl vectors and matrices
    gsl_matrix *Th = gsl_matrix_alloc(p, p);
    gsl_matrix *Sh = gsl_matrix_alloc(p,p);
    gsl_vector *nu = gsl_vector_alloc(p);
    gsl_vector *nuSI = gsl_vector_alloc(p);
    gsl_matrix *Dg1 = gsl_matrix_alloc(p,p);
    gsl_matrix *Dg2 = gsl_matrix_alloc(p,p);
    gsl_matrix *Dg3 = gsl_matrix_alloc(p,p);
    gsl_matrix *nTDi = gsl_matrix_alloc(p, p);
    gsl_vector *vecZero = gsl_vector_calloc(p);
    gsl_vector *Z = gsl_vector_alloc(p); //sampling standard normal
    gsl_vector *U1 = gsl_vector_alloc(p); // matrices to carry out blas multiplications
    gsl_vector *U2 = gsl_vector_alloc(p); //
    gsl_vector *U3 = gsl_vector_alloc(p); //

    // GSL Matrix and vector views
    gsl_matrix_view nuMatrix, SigmahIview;
    gsl_vector_view totprodhv, tothv, XmMview;

    //Bases for vector and matrix views
    double baseXmM[p]; // to view x-muh
    double baseSigmahI[p*p];// to view Sigmah^{-1}

    //Predictions
    double PredProb;
    int newClusI;
    double denomPred[npred];
    int maxy = maxy1[0];
    int maxy2 = maxy;
    double Disthi[ncomp][maxy];
    double StorePD[npred][maxy];
    double sumProb;
    int mode;
    double cdfy; // cdf of response to truncate predictions before maxy
    int start, move;
    double PredTrunc = 0.9999;
    double CDFL[ncomp];
    double CDFU[ncomp];
    double normConst[ncomp];
    for (j = 0; j < ncomp; j++){
        CDFL[j] = 0.0;
        CDFU[j] = 0.0;
        normConst[j] = 0.0;
    }
    double StoreLower[ncomp][maxy];
    double StoreUpper[ncomp][maxy];
    int ll[ncomp];
    int ul[ncomp];
    int Iend = 0;
    double CPDF[npred][ncomp];
    double TMR[npred]; //mean regression from each iteration
    double TQ05R[npred]; //05th percentile regression from each iteration
    double TQ10R[npred]; //10th percentile regression from each iteration
    double TQ15R[npred]; //15th percentile regression from each iteration
    double TQ20R[npred]; //20th percentile regression from each iteration
    double TQ25R[npred]; //25th percentile regression from each iteration
    double TQ50R[npred]; //50th percentile regression from each iteration
    double TQ75R[npred]; //75th percentile regression from each iteration
    double TQ80R[npred]; //80th percentile regression from each iteration
    double TQ85R[npred]; //85th percentile regression from each iteration
    double TQ90R[npred]; //90th percentile regression from each iteration
    double TQ95R[npred]; //95th percentile regression from each iteration

    // Sampler initialization

    // Step 1: Initial allocation to the ncomp components
    for (i = 0; i < n; i++)
        for (h = 0; h < ncomp; h++)
            pdjkj[h][i] = (double) 1/ncomp; //(ncomp-h) or (ncomp-h)*(ncomp-h);

    s = gsl_ran_flat(r,1.0,100000);
    allocation(s,n,ncomp,pdjkj,compAlloc,1);
    for (i = 0; i < n; i++)
        compAlloc[i] = compA[i];

    // Step 2: Find nmembers
    findNmembers(n,ncomp,compAlloc,nmembers);

    // Step 3: Initialize alpha of the DP
    alpha = 2.0;

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

    // Step 4: Update the u_i's from unif(0,pi_k_i)
    if (sampler == 1)
        for (i = 0; i < n; i++)
            u[i] = gsl_ran_flat(r,0.0,pi[compAlloc[i]]);
    else
        for (i = 0; i < n; i++)
            u[i] = 0.0;
    minU = 1.0;
    for (i = 0; i < n; i++)
        if (u[i] < minU) minU = u[i];

    // Step 5: mu_h: initialize covariate cluster means
    s = gsl_ran_flat(r,1.0,100000);
    initMuh(s,X,n,p,ncomp,xbar,xsd,compAlloc,nmembers,muh);

    // Step 6: Correlation vectors unu_h are set to zero
    for (h = 0; h < ncomp; h++)
        for (k = 0; k < p; k++)
            nuh[h][k] = 0.0;

    // Step 7: Initialize parameters Xi_h of the response variables
    s = gsl_ran_flat(r,1.0,100000);
    if (family == 1 || family == 2 || family == 5 || family == 6 || family == 7 || family == 8)
        initGLMOneResLtnt1(s,Y,H,n,ncomp,nRespPars,nmembers,compAlloc,Xi,Ymean[0],family);
    if (family == 3 || family == 4)
        initGLMOneResLtnt2(s,Y,H,n,ncomp,nRespPars,nmembers,compAlloc,Xi,family);

    // Step 8: Inpute latentx from N(0,1)*I
    for (i = 0; i < n; i++){
        for (j = 0; j < nRespPars; j++) XiC[j] = Xi[compAlloc[i]][j];
        calcGLMLimits(Y,H,i,XiC,&lower,&upper,family);
        lower = -9999; upper=9999;
        s = gsl_ran_flat(r,1.0,100000);
        //Rprintf("%s %i %i %f %f %f %f %f \n","test:",i,Y[i],H[i],XiC[0],XiC[1],lower,upper);
        sampleTUN(s,i,0.0,1.0,lower,upper,latentx);
        //Rprintf("%s %f \n","test:",latentx[i]);
    }

    //test:
/*
    Hpred[0]=1; XiC[0] = 2.9; XiC[1] = 0.64;
    calcGLMLimitsPredL(Hpred,5,0,XiC,&lower,5);
    Rprintf("%s %f \n","original:", lower);
    calcGLMLimitsPredL(Hpred,4,0,XiC,&lower,5);
    Rprintf("%s %f \n","original2:", lower);
    calcGLMLimitsPredGP(Hpred,5,0,XiC,&lower,&upper,&CDFL,&CDFU,&normConst);
    Rprintf("%s %f %f %f %f %f \n","new1:",lower,upper,CDFL,CDFU,normConst);
    calcGLMLimitsPredLGP(Hpred,4,0,XiC,&lower,&CDFL,normConst);

    //gsl_cdf_ugaussian_Pinv(cdf_tri_parametric_P2(Y[i]-1,H[i]*Xi[0],Xi[1],Xi[2]));
    Rprintf("%s %f %f %f \n","new2:",cdf_tri_parametric_P2(0,1,3,0.1),
                              cdf_tri_parametric_P(1,2,3,0.2),
                              cdf_tri_parametric_P(1,3,3,0.3));
*/

//gsl_cdf_ugaussian_Pinv(cdf_tri_parametric_P2(Y[i]-1,H[i]*Xi[0],Xi[1],Xi[2]));

//Rprintf("%s %f \n","new2:",cdf_tri_parametric_P2(137,26.071804*84.860290,1.622364,-140.385474));

//#############################################SAMPLER
    for (sw = 0; sw < sweeps; sw++){

        if (sw==0) Rprintf("%i %s \n",sw, "posterior samples...");
        if (((sw+1) % 500)==0) Rprintf("%i %s \n",sw+1, "posterior samples...");

        //Generate all cluster parameters
        nmb = 0; // counts nmembers up to cluster h
        sumPIh = 0.0; // monitor sum pi_h and stop when it is > 1-minU
        h = 0;
        while ((h < ncomp) && (sumPIh <= (1-minU))){
//Rprintf("%s %i %i \n","V...",sw,h);
            // - 1 - Update V_h ~ Beta(nmembers[h]+1,n-nmembers[1:h]+alpha);
            nmb += nmembers[h];
            Vdp[h] = gsl_ran_beta(r, (double) nmembers[h]+1.0, (double) n-nmb+alpha);

            // - 2 - Update pi_h = V[h] prod_{j < h} 1-V[j]
            if (h == 0)
                pi[h] = Vdp[h];
            else if (h > 0 && h < (ncomp-1))
                pi[h] = Vdp[h] * pi[h-1] * (1-Vdp[h-1]) / Vdp[h-1];
            else pi[h] = pi[h-1] * (1-Vdp[h-1]) / Vdp[h-1];
//Rprintf("%s %i %i \n","Th...",sw,h);
            // - 3 - Th
            if (nmembers[h] > 0){
                SetShOneResLtnt(p,n,h,ncomp,X,muh,latentx,nuh,compAlloc,Sh);
                gsl_matrix_add(Sh,Vinv);
                Inverse(p,Sh);
                s = gsl_ran_flat(r,1.0,100000);
                rwish(s,p,Vdf+nmembers[h],Sh,Th);
            }
            else{
                s = gsl_ran_flat(r,1.0,100000);
                rwish(s,p,Vdf,V,Th);
            }

            for (j = 0; j < p; j++)
                for (k = 0; k < p; k++)
                    Th1[h][j*p+k] = gsl_matrix_get(Th,j,k);

            //if (sw==131) Rprintf("%i %i %i %f %f\n",sw,h,nmembers[h],gsl_matrix_get(Sh,0,0),Th1[h][0]);
            //if (sw==131) Rprintf("%i %i %i %i %f %f %lf %lf \n",sw,h,nmembers[h],nmb,Vdp[h],pi[h],sumPIh,1-minU);
//Rprintf("%s %i %i \n","nh...",sw,h);
            // - 4 - nh
            if (nmembers[h] > 0){
                sumxsq = SSQh(n,h,compAlloc,latentx); //set sample totals
                gsl_matrix_memcpy(nTDi,Th); // copy sample from rwish
                gsl_matrix_scale (nTDi,sumxsq); // scale by sum y_i*^2
                gsl_matrix_add(nTDi,SigmaNuInv); // add prior mat
                gHalfInv(p,tol,nTDi); // nTDi^{-1/2}
            }
            else
                gsl_matrix_memcpy(nTDi,SigmaNuHalf); // if n_h = 0, use prior cov mat
            for (k = 0; k < p; k++) // N(0,1)
	            gsl_vector_set(Z,k,gsl_ran_gaussian(r,1));
            gsl_blas_dgemv(CblasNoTrans,1.0,nTDi,Z,0.0,U1); // sample from MV normal
            if (nmembers[h] > 0){
                gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,nTDi,nTDi,0.0,Dg1); // obtains nTDi^{-1}
                SetSampleTotNu(p,n,h,ncomp,totprodhs,compAlloc,X,muh,latentx); //set sample totals
                totprodhv = gsl_vector_view_array(totprodhs,p); //view sample totals
                gsl_blas_dgemv(CblasNoTrans,1.0,Th,&totprodhv.vector,0.0,U2); // part of mean
                gsl_vector_add(U2,SIMnu); // add prior contribution
                gsl_blas_dgemv(CblasNoTrans,1.0,Dg1,U2,0.0,U3); // mean
            }
            else
                gsl_vector_memcpy(U3,MuNu);
            gsl_vector_add(U1,U3); //posterior sample
            //gsl_vector_set_all(U1,0.0);//PRMM
            gsl_vector_memcpy(nu,U1); //store 1
            for (k = 0; k < p; k++)
                nuh[h][k] = gsl_vector_get(nu,k); // store 2
//Rprintf("%s %i %i \n","mh...",sw,h);
            // - 5 - mu_h
            if (nmembers[h] > 0){
                gsl_matrix_memcpy(nTDi,Th); // copy sample from rwish
                gsl_matrix_scale (nTDi,nmembers[h]); // scale by n_h
                gsl_matrix_add(nTDi,SigmaMuInv); // add inv of prior cov mat
                gHalfInv(p,tol,nTDi); // nTDi^{-1/2}
            }
            else
                gsl_matrix_memcpy(nTDi,SigmaMuHalf); // if n_h = 0, use prior cov mat
            for (k = 0; k < p; k++) // N(0,1)
	            gsl_vector_set(Z,k,gsl_ran_gaussian(r,1));
            gsl_blas_dgemv(CblasNoTrans,1.0,nTDi,Z,0.0,U1);  // sample from MV normal
            if (nmembers[h] > 0){
                SetSampleTotMu(p,n,h,ncomp,toths,compAlloc,X,nuh,latentx);
                tothv = gsl_vector_view_array(toths,p);
                gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,nTDi,nTDi,0.0,Dg1); // obtains nTDi^{-1}
                gsl_blas_dgemv(CblasNoTrans,1.0,Th,&tothv.vector,0.0,U2); // T toth
                gsl_vector_add(U2,SIMmu);
                gsl_blas_dgemv(CblasNoTrans,1.0,Dg1,U2,0.0,U3); // (nT+D^{-1})^{-1} (T toth + D^{-1} xbar)
            }
            else
                gsl_vector_memcpy(U3,MuMu);
            gsl_vector_add (U1,U3);
            for (k = 0; k < p; k++)
                muh[h][k] = gsl_vector_get(U1,k);
//Rprintf("%s %i %i \n","xih...",sw,h);
            // - 6 - Xi_h
            // Current and proposed parameters
            for (k = 0; k < nRespPars; k++)
                XiC[k] = Xi[h][k];

            s = gsl_ran_flat(r,1.0,100000);
            propose(s,XiC,XiP,nRespPars,prec,family);

            XiLC = 0.0;
            XiLP = 0.0;

//Rprintf("%i %i %f %f %f %f %f %f \n",sw,h,XiC[0],XiC[1],XiC[2],XiP[0],XiP[1],XiP[2]);

            gsl_matrix_memcpy(Dg2,Th); // copy Th: sample from rwish
            Inverse(p,Dg2);
            nuMatrix = gsl_matrix_view_vector(nu,p,1);
            gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,&nuMatrix.matrix,&nuMatrix.matrix,1.0,Dg2); // Sigma_h
            gsl_matrix_memcpy(Dg3,Dg2);
            Inverse(p,Dg3); // Dg3 = Sigma^{-1}

            // Store Dg2 in Sigmah and Dg3 in SigmahI
            for (j = 0; j < p; j++){
                for (k = 0; k < p; k++){
                    Sigmah[h][j*p+k] = gsl_matrix_get(Dg2,j,k);
                    SigmahI[h][j*p+k] = gsl_matrix_get(Dg3,j,k);
                }
            }

            // Calculate nu^T Sigma^{-1} and nu^T Sigma^{-1} nu
            gsl_blas_dgemv(CblasNoTrans,1.0,Dg3,nu,0.0,nuSI);
            gsl_blas_ddot(nuSI,nu,&nuSInu);
            if (nuSInu >= 1.0) nuSInu = 0.9999999999;

            // Store for use in p(dj=kj)
            storenuSInu[h] = nuSInu;
            for (k = 0; k < p; k++)
                storenuSI[h][k] = gsl_vector_get(nuSI,k);

//Rprintf("%s %i %i\n","mh algo...",sw,h);
            //MH algorithm
            temp3 = 0;
            if (family < 8){
                for (k = 0; k < nRespPars; k++)
                    if (XiP[k] > 0.00001) temp3 +=1;
            }
            else if (family == 8) temp3 = nRespPars;

//if (sw==54364) Rprintf("%i %i %f %f %f %f %f %f\n",sw,h,XiC[0],XiC[1],XiC[2],XiP[0],XiP[1],XiP[2]);
            if (temp3 == nRespPars){
                for (i = 0; i < n; i++){
                    if (compAlloc[i]==h && exp(XiLP) > 0){
                        nuSIxmm = 0.0;
                        for (k = 0; k < p; k++)
                            nuSIxmm += storenuSI[h][k]*(X[k*n+i]-muh[h][k]);
                        calcGLMLimits(Y,H,i,XiC,&cymoC,&cyC,family);
                        calcGLMLimits(Y,H,i,XiP,&cymoP,&cyP,family);
                        XiLC += log(gsl_cdf_ugaussian_P((cyC-nuSIxmm)/sqrt(1-nuSInu)) - gsl_cdf_ugaussian_P((cymoC-nuSIxmm)/sqrt(1-nuSInu)));
                        XiLP += log(gsl_cdf_ugaussian_P((cyP-nuSIxmm)/sqrt(1-nuSInu)) - gsl_cdf_ugaussian_P((cymoP-nuSIxmm)/sqrt(1-nuSInu)));
                        //Rprintf("%i %i %i %i %f %f %f %f %f %f %f %f %f %f %f %f \n",
                        //sw,h,i,Y[i],H[i],XiC[0],XiC[1],cymoC,cyC,XiP[0],XiP[1],XiP[2],cymoP,cyP,XiLC,XiLP);
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
                    NAcp = exp(XiLP)*gsl_ran_gamma_pdf(XiP[0],alphaXi[0],1/betaXi[0])*
                                     gsl_ran_gaussian_pdf(XiP[1]-alphaXi[1],betaXi[1])*
                                     gsl_ran_gamma_pdf(XiC[0],prec[0]*XiP[0]*XiP[0],1/(prec[0]*XiP[0]));//random walk for other parameter
                                     //gsl_ran_gaussian_pdf(XiP[1]-XiC[1],prec[1]);
                    DAcp = exp(XiLC)*gsl_ran_gamma_pdf(XiC[0],alphaXi[0],1/betaXi[0])*
                                     gsl_ran_gaussian_pdf(XiC[1]-alphaXi[1],betaXi[1])*
                                     gsl_ran_gamma_pdf(XiP[0],prec[0]*XiC[0]*XiC[0],1/(prec[0]*XiC[0]));
                                     //gsl_ran_gaussian_pdf(XiP[1]-XiC[1],prec[1]);
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
                temp = NAcp/DAcp;
                if (temp > 1.0) Acp = 1.0;else Acp = temp;
                temp = gsl_ran_flat(r,0.0,1.0);
                if (Acp > temp)
                    for (k = 0; k < nRespPars; k++) Xi[h][k] = XiP[k];
                else
                    for (k = 0; k < nRespPars; k++) Xi[h][k] = XiC[k];
            }
            else
                for (k = 0; k < nRespPars; k++) Xi[h][k] = XiC[k];
            for (k = 0; k < nRespPars; k++) XiC[k] = Xi[h][k];

//if (sw==54364) Rprintf("%s %i %i \n","impute latent...",sw,h);
            // - 7 - Impute latent y^*
            for (i = 0; i < n; i++){
                if (compAlloc[i]==h){
                    nuSIxmm = 0.0;
                    for (k = 0; k < p; k++)
                        nuSIxmm += storenuSI[h][k]*(X[k*n+i]-muh[h][k]);
                    calcGLMLimits(Y,H,i,XiC,&lower,&upper,family);
                    s = gsl_ran_flat(r,1.0,100000);
                    sampleTUN(s,i,nuSIxmm,sqrt(1-nuSInu),lower,upper,latentx);
                }
            }
            sumPIh += pi[h];
            h++;
        }// h loop
        nUpdated = h;

        // - 8 - Update the u_i's from unif(0,pi_k_i)
        if (sampler == 1){
            for (i = 0; i < n; i++)
                u[i] = gsl_ran_flat(r,0.0,pi[compAlloc[i]]);
            minU = 1.0;
            for (i = 0; i < n; i++)
                if (u[i] < minU) minU = u[i];
        }

//if (sw==54364) Rprintf("%s %i %i\n","pdelta...",sw,h);
        // - 9 - P(delta_j = k_j)
        for (h = 0; h < nUpdated; h++){
            // View Sigmah^{-1}
            for (k = 0; k < p*p; k++)
                baseSigmahI[k] = SigmahI[h][k];
            SigmahIview = gsl_matrix_view_array(baseSigmahI,p,p);
            for (k = 0; k < nRespPars; k++) XiC[k] = Xi[h][k];
//Rprintf("%i %i %f %f  \n",sw,h,XiC[0],XiC[1]);
            for (i = 0; i < n; i++){
                if (pi[h] > u[i]){
                    nuSIxmm = 0.0;
                    for (k = 0; k < p; k++)
                        nuSIxmm += storenuSI[h][k]*(X[k*n+i]-muh[h][k]);
//Rprintf("%i %i %i %i %f %f %f  \n",sw,h,i,Y[i],H[i],XiC[0],XiC[1]);
//if (XiC[1]>0.23 || i==0)
                    calcGLMLimits(Y,H,i,XiC,&lower,&upper,family);
//if (sw==8 && h>35 && h<41)
//Rprintf("%i %i %i %i %f %f %f %f %f  \n",sw,h,i,Y[i],H[i],XiC[0],XiC[1],lower,upper);
                    temp = gsl_cdf_ugaussian_P((upper-nuSIxmm)/sqrt(1-storenuSInu[h])) - gsl_cdf_ugaussian_P((lower-nuSIxmm)/sqrt(1-storenuSInu[h]));
                    for (k = 0; k < p; k++)
                        baseXmM[k] = X[k*n+i]-muh[h][k];
                    XmMview = gsl_vector_view_array(baseXmM,p);
                    temp2 = logMVNormalpdf3(p,&XmMview.vector,vecZero,&SigmahIview.matrix);
                    if (sampler == 1) pdjkj[h][i] = exp(log(temp)+temp2);
                    else if (sampler == 2) pdjkj[h][i] = exp(log(temp)+temp2+log(pi[h]));
//if (sw==12 && h==9 && i==499) Rprintf("%i %i %i %f %f %f %f %f %f %f %f \n",
//                                        sw,h,i,XiC[0],XiC[1],temp,pdjkj[h][i],lower,upper,nuSIxmm,storenuSInu[h]);
                }
                else
                    pdjkj[h][i] = 0.0;
                //if (sw==2847 && i==149 && h==3)
                //Rprintf("%i %i %i %i %f %f %f %f %f %f %f %f %f %f %f\n",
                //sw,h,i,Y[i],H[i],XiC[0],XiC[1],XiC[2],nuSIxmm,storenuSInu[h],upper,lower,temp,temp2,pdjkj[h][i]);


                //Rprintf("%i %f ", h, pdjkj[h][i]);
            }
        }
        for (h = nUpdated; h < ncomp; h++)
            for (i = 0; i < n; i++)
                pdjkj[h][i] = 0.0;
//if (sw==8 && h>35 && h<41)
//if (sw==54364) Rprintf("%s %i \n","alloc...",sw);
        // - 10 - Update allocation to
        s = gsl_ran_flat(r,1.0,100000);
        allocation(s,n,ncomp,pdjkj,compAlloc,sw);
//if (sw==8 && h>35 && h<41)
//Rprintf("%s %i \n","nmem...",sw);
        // - 11 - Find nmembers
        findNmembers(n,ncomp,compAlloc,nmembers);

        // - 12a - Label switching
        s = gsl_ran_flat(r,1.0,100000);
        labelSwtchA(s,n,p,ncomp,nRespPars,Th1,Sigmah,SigmahI,nuh,muh,Xi,storenuSInu,storenuSI,nmembers,compAlloc,pi);

        // - 12b - Label switching
        s = gsl_ran_flat(r,1.0,100000);
        //if (sw==131) Rprintf("%s %li \n","BL",s);
        labelSwtchB(s,n,p,ncomp,nRespPars,Th1,Sigmah,SigmahI,nuh,muh,Xi,storenuSInu,storenuSI,nmembers,compAlloc,Vdp);

        // - 13 - Update concentration parameter
        s = gsl_ran_flat(r,1.0,100000);
        newalpha = updateAlpha(s,n,ncomp,Alphaa,Alphab,TruncAlpha,nmembers,alpha);
        alpha = newalpha;
//if (sw==54364) Rprintf("%s %i \n","pred...",sw);
        // - 14 - Predictions

        if (((sw - burn) >= 0) && (((sw - burn ) % thin) == 0) && (npred > 0)){

            newClusI = 0;
            for (i = 0; i < npred; i++) denomPred[i] = 0.0;
            for (k = 0; k < maxy; k++) for (i = 0; i < npred; i++) StorePD[i][k] = 0.0;
            for (h = 0; h < ncomp; h++) if (nmembers[h] > 0) upperL = h;
            upperL += 2; //one to work ok and one for the empty cluster; check for in between zeros

            for (h = 0; h < upperL; h++){
                if (nmembers[h] > 0 || newClusI == 0){
                    if (nmembers[h] == 0) newClusI = 1;
                    for (k = 0; k < p*p; k++)
                        baseSigmahI[k] = SigmahI[h][k];
                    SigmahIview = gsl_matrix_view_array(baseSigmahI,p,p);
                    for (i = 0; i < npred; i++){
                        for (k = 0; k < p; k++)
                            baseXmM[k] = Xpred[k*npred+i]-muh[h][k];
                        XmMview = gsl_vector_view_array(baseXmM,p);
                        temp2 = logMVNormalpdf3(p,&XmMview.vector,vecZero,&SigmahIview.matrix);
                        CPDF[i][h] = exp(temp2);
                        if (nmembers[h] > 0) PredProb = CPDF[i][h]*nmembers[h];
                        else if (nmembers[h] == 0) PredProb = CPDF[i][h]*alpha;
                        denomPred[i] += PredProb;
                    }
                }
            }

            // posterior predictive distribution
            Iend = 0;
            for (i = 0; i < npred; i++){
                for (h = 0; h < ncomp; h++) for (k = 0; k < maxy; k++) Disthi[h][k] = 0.0;
                newClusI = 0;
                for (h = 0; h < upperL; h++){
                    if (nmembers[h] > 0 || newClusI == 0){
                        if (nmembers[h] == 0) newClusI = 1;
                        for (k = 0; k < nRespPars; k++)
                            XiC[k] = Xi[h][k];
                        nuSIxmm = 0.0;
                        for (k = 0; k < p; k++)
                            nuSIxmm += storenuSI[h][k]*(Xpred[k*npred+i]-muh[h][k]);

                        //start = Hpred[i]*XiC[0];
                        //if (family == 3) start = Hpred[i]*XiC[0]/XiC[1];
                        //if (family == 4) start = Hpred[i]*XiC[0]/(XiC[0] + XiC[1]);
                        //while (start >= maxy) start = maxy/2;

                        start = 0;

                        if (i == 0 || allEqlI[0] == 0){
                            if (family < 5) calcGLMLimitsPred(Hpred,start,i,XiC,&lower,&upper,family);
                            else if (family==5) calcGLMLimitsPredGP(Hpred,start,i,XiC,&lower,&upper,&CDFL[h],&CDFU[h],&normConst[h]);
                            else if (family==6) calcGLMLimitsPredCP(Hpred,start,i,XiC,&lower,&upper,&CDFL[h],&CDFU[h],&normConst[h]);
                            else if (family==7) calcGLMLimitsPredHP(Hpred,start,i,XiC,&lower,&upper,&CDFL[h],&CDFU[h],&normConst[h]);
                            else if (family==8) calcGLMLimitsPredCTP(Hpred,start,i,XiC,&lower,&upper,&CDFL[h],&CDFU[h],&normConst[h]);
                            StoreLower[h][start] = lower;
                            StoreUpper[h][start] = upper;
                            if (family>=5) if (gsl_isinf(normConst[h]) || normConst[h]==0.0 || gsl_isnan(CDFL[h]) || gsl_isnan(-CDFL[h])
                                || gsl_isnan(CDFU[h]) || gsl_isnan(-CDFU[h])) Iend=1;
                        }
                        else{
                            lower = StoreLower[h][start];
                            upper = StoreUpper[h][start];
                        }
                        temp = gsl_cdf_ugaussian_P((upper-nuSIxmm)/sqrt(1-storenuSInu[h]))-
                               gsl_cdf_ugaussian_P((lower-nuSIxmm)/sqrt(1-storenuSInu[h]));
                        upper2 = upper;
                        cdfy = temp;
                        if (nmembers[h] > 0) Disthi[h][start] = temp * CPDF[i][h] * nmembers[h];
                        else if (nmembers[h] == 0) Disthi[h][start] = temp * CPDF[i][h] * alpha;
                        move = 1;

//if (sw==40064) fprintf(out_file7,"%s %i %i %i %i %i %i %f %f %f %f %f %f %f %f %f %f %f %f \n","start",
//sw,i,h,start,move,k,XiC[0],XiC[1],XiC[2],lower,upper,CDFL[h],CDFU[h],normConst[h],cdfy,nuSIxmm,storenuSInu[h],temp);

                        if (family==2 || family==4) maxy2 = Hpred[i]+1;
                        while(cdfy<PredTrunc && ((start-move)>=0 || (start+move)<maxy2)){
                            if (start >= move){
                                k = start - move;
                                upper = lower;
                                if ((i == 0) || (k < ll[h]) || (allEqlI[0] == 0)){
                                    if (family < 5) calcGLMLimitsPredL(Hpred,k,i,XiC,&lower,family);
                                    else if (family == 5) calcGLMLimitsPredLGP(Hpred,k,i,XiC,&lower,&CDFL[h],normConst[h]);
                                    else if (family == 6) calcGLMLimitsPredLCP(Hpred,k,i,XiC,&lower,&CDFL[h],normConst[h]);
                                    else if (family == 7) calcGLMLimitsPredLHP(Hpred,k,i,XiC,&lower,&CDFL[h],normConst[h]);
                                    else if (family == 8) calcGLMLimitsPredLCTP(Hpred,k,i,XiC,&lower,&CDFL[h],normConst[h]);
                                    StoreLower[h][k] = lower;
                                    ll[h] = k;
                                }
                                else lower = StoreLower[h][k];
                                temp = gsl_cdf_ugaussian_P((upper-nuSIxmm)/sqrt(1-storenuSInu[h]))-
                                       gsl_cdf_ugaussian_P((lower-nuSIxmm)/sqrt(1-storenuSInu[h]));
                                cdfy += temp;
                                if (nmembers[h] > 0) Disthi[h][k] = temp * CPDF[i][h] * nmembers[h];
                                else if (nmembers[h] == 0) Disthi[h][k] = temp * CPDF[i][h] * alpha;

//if (sw==40064) fprintf(out_file7,"%s %i %i %i %i %i %i %f %f %f %f %f %f %f %f %f %f %f %f \n","start-move",
//sw,i,h,start,move,k,XiC[0],XiC[1],XiC[2],lower,upper,CDFL[h],CDFU[h],normConst[h],cdfy,nuSIxmm,
//storenuSInu[h],temp);

                            }
                            if ((start+move) < maxy2){
                                k = start + move;
                                lower2 = upper2;
                                if ((i == 0) || (k > ul[h]) || (allEqlI[0] == 0)){
                                    if (family < 5) calcGLMLimitsPredU(Hpred,k,i,XiC,&upper2,family);
                                    else if (family == 5) calcGLMLimitsPredUGP(Hpred,k,i,XiC,&upper2,&CDFU[h],normConst[h]);
                                    else if (family == 6) calcGLMLimitsPredUCP(Hpred,k,i,XiC,&upper2,&CDFU[h],normConst[h]);
                                    else if (family == 7) calcGLMLimitsPredUHP(Hpred,k,i,XiC,&upper2,&CDFU[h],normConst[h]);
                                    else if (family == 8) calcGLMLimitsPredUCTP(Hpred,k,i,XiC,&upper2,&CDFU[h],normConst[h]);
                                    StoreUpper[h][k] = upper2;
                                    ul[h] = k;
                                }
                                else upper2 = StoreUpper[h][k];
                                temp = gsl_cdf_ugaussian_P((upper2-nuSIxmm)/sqrt(1-storenuSInu[h]))-
                                       gsl_cdf_ugaussian_P((lower2-nuSIxmm)/sqrt(1-storenuSInu[h]));
                                cdfy += temp;
                                if (nmembers[h] > 0) Disthi[h][k] = temp * CPDF[i][h] * nmembers[h];
                                else if (nmembers[h] == 0) Disthi[h][k] = temp * CPDF[i][h] * alpha;
//if (sw==40064) fprintf(out_file7,"%s %i %i %i %i %i %i %f %f %f %f %f %f %f %f %f %f %f %f \n","start+move",
//sw,i,h,start,move,k,XiC[0],XiC[1],XiC[2],lower2,upper2,CDFL[h],CDFU[h],normConst[h],cdfy,
//nuSIxmm,storenuSInu[h],temp);

                            }
                            move++;
                        }
                        for (k = 0; k <= (start-move); k++) Disthi[h][k] = 0.0;
                        for (k = (start+move); k < maxy; k++) Disthi[h][k] = 0.0;
/*
                        cdfy = 0.0;
                        for (k = 0; k < maxy; k++){

                            if(cdfy <=0.9999999999999999999999999999999999999999999999999999999 && k < maxy){

                                calcGLMLimitsPred(Hpred,k,i,XiC,&lower,&upper,family);
                                temp = gsl_cdf_ugaussian_P((upper-nuSIxmm)/sqrt(1-storenuSInu[h]))-
                                       gsl_cdf_ugaussian_P((lower-nuSIxmm)/sqrt(1-storenuSInu[h]));
                                cdfy += temp;
                            }
                            else temp=0;
                            Disthi[h][k] = temp * PredProb[i][h]/denomPred[i];
                        }
*/
                    }
                }
                for (k = 0; k < maxy; k++)
                   for (h = 0; h < upperL; h++)
                       StorePD[i][k] += Disthi[h][k]/denomPred[i];
            }
//if (sw >= 592 && sw < 600) Rprintf("%s %i %i %i %f %f %f %f %f %f %f %f %f %i %i \n","pred trunc2...",sw,i,h,XiC[0],Hpred[i-1],Xpred[i-1],upper,lower,
//nuSIxmm,sqrt(1-storenuSInu[h]),
//gsl_cdf_ugaussian_P((upper-nuSIxmm)/sqrt(1-storenuSInu[h]))-gsl_cdf_ugaussian_P((lower-nuSIxmm)/sqrt(1-storenuSInu[h])),
//cdfy,k,maxy);
            // posterior summaries
            if (Iend == 0){
                // Mean
                for (i = 0; i < npred; i++) TMR[i] = 0.0;

                for (i = 0; i < npred; i++)
                    for (k = 0; k < maxy; k++)
                        TMR[i] += (double) (StorePD[i][k] * k);

                for (i = 0; i < npred; i++)
                    meanReg[i] += (double) TMR[i]/(nSamples);

                //Quantiles
                for (i = 0; i < npred; i++){
                    k=0;
                    sumProb = StorePD[i][k];
                    while (sumProb <= 0.05){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ05R[i] = k;
                    Q05Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.10){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ10R[i] = k;
                    Q10Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.15){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ15R[i] = k;
                    Q15Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.20){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ20R[i] = k;
                    Q20Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.25){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ25R[i] = k;
                    Q25Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.50){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ50R[i] = k;
                    Q50Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.75){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ75R[i] = k;
                    Q75Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.80){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ80R[i] = k;
                    Q80Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.85){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ85R[i] = k;
                    Q85Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.90){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ90R[i] = k;
                    Q90Reg[i] += (double) k/nSamples;
                    while (sumProb <= 0.95){
                        k++;
                        sumProb += StorePD[i][k];
                    }
                    TQ95R[i] = k;
                    Q95Reg[i] += (double) k/nSamples;
                }
                // Mode
                for (i = 0; i < npred; i++){
                    mode = gsl_stats_max_index(StorePD[i],1,maxy);
                    modeReg[i] += (double) (mode)/nSamples;
                }
                //Mean and Var of the probabilities
                for (i = 0; i < npred; i++){
                    for (k = 0; k < maxy; k++){
                        denReg[k+maxy*i] += (double) StorePD[i][k]/nSamples;
                        denVar[k+maxy*i] += (double) StorePD[i][k]*StorePD[i][k]/nSamples;
                        //if (k<10) Rprintf("%s %i %i %f %f \n","pred,num,den,store",i,k,denReg[k+maxy*i], StorePD[i][k]);
                    }
                }
            }
        }
        // - 15 - Write to files
        if (((sw - burn) >= 0) && (((sw - burn ) % thin) == 0) && (WF == 1) && (Iend==0)){

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < (p*p); j++)
                    fprintf(out_file1, "%f ", Th1[h][j]);
                fprintf(out_file1, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < (p*p); j++)
                    fprintf(out_file2, "%f ", Sigmah[h][j]);
                fprintf(out_file2, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < (p*p); j++)
                    fprintf(out_file3, "%f ", SigmahI[h][j]);
                fprintf(out_file3, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < p; j++)
                    fprintf(out_file4, "%f ", nuh[h][j]);
                fprintf(out_file4, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < p; j++)
                    fprintf(out_file5, "%f ", muh[h][j]);
                fprintf(out_file5, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < nRespPars; j++)
                    fprintf(out_file6, "%f ", Xi[h][j]);
                fprintf(out_file6, "\n");
            }

            fprintf(out_file7, "%f \n", alpha);

            for (i = 0; i < n; i++)
                fprintf(out_file8, "%i ", compAlloc[i]);
            fprintf(out_file8, "\n");

            for (h = 0; h < ncomp; h++)
                fprintf(out_file9, "%i ", nmembers[h]);
            fprintf(out_file9, "\n");

            fprintf(out_file10, "%i \n", nUpdated);

            for (i = 0; i < npred; i++)
                fprintf(out_file11, "%f ", TMR[i]);
            fprintf(out_file11, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file12, "%f ", TQ05R[i]);
            fprintf(out_file12, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file13, "%f ", TQ10R[i]);
            fprintf(out_file13, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file14, "%f ", TQ15R[i]);
            fprintf(out_file14, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file15, "%f ", TQ20R[i]);
            fprintf(out_file15, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file16, "%f ", TQ25R[i]);
            fprintf(out_file16, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file17, "%f ", TQ50R[i]);
            fprintf(out_file17, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file18, "%f ", TQ75R[i]);
            fprintf(out_file18, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file19, "%f ", TQ80R[i]);
            fprintf(out_file19, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file20, "%f ", TQ85R[i]);
            fprintf(out_file20, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file21, "%f ", TQ90R[i]);
            fprintf(out_file21, "\n");

            for (i = 0; i < npred; i++)
                fprintf(out_file22, "%f ", TQ95R[i]);
            fprintf(out_file22, "\n");

        }
        if (Iend == 1) sw -= 1;
        if (sw==(sweeps-1) && (!((sw+1) % 500)==0)) Rprintf("%i %s \n",sw+1, "posterior samples...");
    }//end of sw

    //Free up random number
    gsl_rng_free (r);

    //Close files
    fclose(out_file1); fclose(out_file2); fclose(out_file3); fclose(out_file4); fclose(out_file5);
    fclose(out_file6); fclose(out_file7); fclose(out_file8); fclose(out_file9); fclose(out_file10);
    fclose(out_file11); fclose(out_file12); fclose(out_file13); fclose(out_file14); fclose(out_file15);
    fclose(out_file16); fclose(out_file17); fclose(out_file18); fclose(out_file19); fclose(out_file20);
    fclose(out_file21); fclose(out_file22);

    //Free up gsl matrices
    gsl_matrix_free(V);gsl_matrix_free(Vinv);
    gsl_vector_free(MuNu); gsl_matrix_free(SigmaNuInv); gsl_matrix_free(SigmaNuHalf); gsl_vector_free(SIMnu);
    gsl_vector_free(MuMu); gsl_matrix_free(SigmaMuInv); gsl_matrix_free(SigmaMuHalf); gsl_vector_free(SIMmu);
    gsl_matrix_free(Th); gsl_matrix_free(Sh);
    gsl_vector_free(nu);
    gsl_vector_free(nuSI);
    gsl_matrix_free(Dg1); gsl_matrix_free(Dg2); gsl_matrix_free(Dg3);
    gsl_matrix_free(nTDi);
    gsl_vector_free(vecZero);
    gsl_vector_free(Z); gsl_vector_free(U1); gsl_vector_free(U2); gsl_vector_free(U3);
}
