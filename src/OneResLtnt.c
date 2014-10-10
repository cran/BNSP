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
#include "matalg.h"
#include "pdfs.h"
#include "sampling.h"
#include "other.functions.h"
#include "mathm.h"
#include "SpecOneRL.h"
//include "cubature.h"

void OneResLtnt(int *seed1, double *X, int *Y, double *H,
                int *sweeps1, int *burn1, int *ncomp1, int *n1, int *p1,
                double *Vprior, double *Vdf1,
                double *munu, double *sigmanu,
                double *mumu, double *sigmamu,
                double *alphagamma1, double *betagamma1,
                double *Alphaa1, double *Alphab1, double *TruncAlpha1,
                double *xbar, double *xsd, double *Ymean,
                int *family1, int *sampler1,
                int *npred1, double *Xpred, double *Hpred, int *maxy1, double *meanReg, double *medianReg,
                double *q1Reg, double *q3Reg, double *modeReg, char **WorkingDir, int *WF1)
{
    gsl_set_error_handler_off();

    // Random number generator initialization
    int seed = seed1[0]; //random number generator seed
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);
    unsigned long int s; //for calling random number generating functions

    // Sweeps, burn-in period
    int sweeps = sweeps1[0]; //number of posterior samples
    int burn = burn1[0]; // integer burn-in period

    //Family
    int family = family1[0]; //1=poisson, 2=binomial
    int sampler = sampler1[0]; //1=slice, 2=truncated

    //Declare variables for loops: h for components, k,j for covariates, i for subjects, sw for sweeps
    int h, i, j, k, sw;

    // Dimensions
    int n = n1[0]; // number of observations
    int p = p1[0]; // number of covariates
    int ncomp = ncomp1[0]; //number of clusters/components
    int npred = npred1[0]; //number of predictions

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
    double alphagamma = alphagamma1[0];
    double betagamma = betagamma1[0];
    // - 5 - concentration parameter alpha ~ Gamma(a,b) I[alpha>c]
    double Alphaa = Alphaa1[0];
    double Alphab = Alphab1[0];
    double TruncAlpha = TruncAlpha1[0];

    // Specify directory
    int WF = WF1[0]; // indicator: 1 = write files, 0 = no files.
    const int MAX_PATH = 300;
    char path_name[MAX_PATH + 1];

    // Open files
    FILE *out_file1, *out_file2, *out_file3, *out_file4, *out_file5, *out_file6, *out_file7, *out_file8, *out_file9, *out_file10,
         *out_file11;

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
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.gammah.txt");
    out_file6 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.alpha.txt");
    out_file7 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.compAlloc.txt");
    out_file8 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.nmembers.txt");
    out_file9 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.Updated.txt");
    out_file10 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.PD.txt");
    out_file11 = fopen(path_name, "wt");

    // Define quantities of interest
    double Th1[ncomp][p*p]; // to store Th
    double Sigmah[ncomp][p*p]; // to store \Sigma_{h}
    double SigmahI[ncomp][p*p]; // to store \Sigma_{h}^{-1}
    double nuh[ncomp][p]; // to store nu
    double muh[ncomp][p]; // mu_h
    double gamma[ncomp]; // GLM mean
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
    double gammaC, gammaP, gammaLC, gammaLP, cyC, cyP, cymoC, cymoP, Acp, NAcp, DAcp; //MH  step
    double nuSInu, nuSIxmm; //nu^T S^-1 nu and nu^T SI (x-mu)
    double storenuSInu[ncomp]; // to store nu^T Sigma^(-1) nu
    double storenuSI[ncomp][p]; // to store nu^T Sigma^(-1)
    double totprodhs[p];  // total of clusters: (x-mu)*x1
    double sumxsq; // sum of cluster latentx[i]^2
    double toths[p];  // total of clusters: x-mu-y1^* * nu
    double temp, temp2;
    int temp3, temp4;
    double lower, upper;

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
    double PredProb[npred][ncomp];
    int newClusI[npred];
    double denomPred[npred];
    int maxy = maxy1[0];
    double Disthi[ncomp][maxy];
    double StorePD[npred][maxy];
    double sumProb;
    int mode;

    // Sampler initialization

    // Step 1: Initial allocation to the ncomp components
    for (i = 0; i < n; i++)
        for (h = 0; h < ncomp; h++)
            pdjkj[h][i] = (double) 1/ncomp; //(ncomp-h) or (ncomp-h)*(ncomp-h);

    s = gsl_ran_flat(r,1.0,100000);
    allocation(s,n,ncomp,pdjkj,compAlloc);

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

    // Step 7: Initialize parameters gamma_h of the response variables
    s = gsl_ran_flat(r,1.0,100000);
    initGLMOneResLtnt(s,Y,H,n,p,ncomp,nmembers,compAlloc,gamma,Ymean[0],family);

    // Step 8: Inpute latentx from N(0,1)*I
    for (i = 0; i < n; i++){
        gammaC = gamma[compAlloc[i]];
        calcGLMLimits(Y,H,i,gammaC,&lower,&upper,family);
        s = gsl_ran_flat(r,1.0,100000);
        sampleTUN(s,i,0.0,1.0,lower,upper,latentx);
    }

    //#############################################SAMPLER
    for (sw = 0; sw < sweeps; sw++){

        if (sw==0) Rprintf("%i %s \n",sw, "posterior samples...");
        if (((sw+1) % 1000)==0) Rprintf("%i %s \n",sw+1, "posterior samples...");

        //Generate all cluster parameters
        nmb = 0; // counts nmembers up to cluster h
        sumPIh = 0.0; // monitor sum pi_h and stop when it is > 1-minU
        h = 0;
        while ((h < ncomp) && (sumPIh <= (1-minU))){

            // - 1 - Update V_h ~ Beta(nmembers[h]+1,n-nmembers[1:h]+alpha);
            nmb += nmembers[h];
            Vdp[h] = gsl_ran_beta(r, (double) nmembers[h]+1.0, (double) n-nmb+alpha);

            // - 2 - Update pi_h = V[h] prod_{j < h} 1-V[j]
            if (h == 0)
                pi[h] = Vdp[h];
            else if (h > 0 && h < (ncomp-1))
                pi[h] = Vdp[h] * pi[h-1] * (1-Vdp[h-1]) / Vdp[h-1];
            else pi[h] = pi[h-1] * (1-Vdp[h-1]) / Vdp[h-1];

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
            gsl_vector_memcpy(nu,U1); //store 1
            for (k = 0; k < p; k++)
                nuh[h][k] = gsl_vector_get(nu,k); // store 2

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

            // - 6 - gamma_h
            // Current and proposed gammas and their likelihoods
            gammaC = gamma[h];
            if (family == 1) gammaP = gsl_ran_gamma(r,gammaC*gammaC,1/gammaC); //has mean gammaC and variance=1
            else if (family == 2) gammaP = gsl_ran_beta(r,gammaC+1.0,2.0-gammaC); // has mode gammaC
            gammaLC = 0.0;
            gammaLP = 0.0;

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

            //MH algorithm
            if (gammaP > 0.00001){
                for (i = 0; i < n; i++){
                    if (compAlloc[i]==h){
                        nuSIxmm = 0.0;
                        for (k = 0; k < p; k++)
                            nuSIxmm += gsl_vector_get(nuSI,k)*(X[k*n+i]-muh[h][k]);
                        calcGLMLimits(Y,H,i,gammaC,&cymoC,&cyC,family);
                        calcGLMLimits(Y,H,i,gammaP,&cymoP,&cyP,family);
                        gammaLC += log(gsl_cdf_ugaussian_P((cyC-nuSIxmm)/sqrt(1-nuSInu)) - gsl_cdf_ugaussian_P((cymoC-nuSIxmm)/sqrt(1-nuSInu)));
                        gammaLP += log(gsl_cdf_ugaussian_P((cyP-nuSIxmm)/sqrt(1-nuSInu)) - gsl_cdf_ugaussian_P((cymoP-nuSIxmm)/sqrt(1-nuSInu)));
                    }
                }
                if (family == 1){
                    NAcp = exp(gammaLP)*gsl_ran_gamma_pdf(gammaP,alphagamma,1/betagamma)*
                           gsl_ran_gamma_pdf(gammaC,gammaP*gammaP,1/gammaP);
                    DAcp = exp(gammaLC)*gsl_ran_gamma_pdf(gammaC,alphagamma,1/betagamma)*
                           gsl_ran_gamma_pdf(gammaP,gammaC*gammaC,1/gammaC);
                }
                else if (family == 2){
                    NAcp = exp(gammaLP)*gsl_ran_beta_pdf(gammaP,alphagamma,betagamma)*
                           gsl_ran_beta_pdf(gammaC,gammaP+1,2-gammaP);
                    DAcp = exp(gammaLC)*gsl_ran_beta_pdf(gammaC,alphagamma,betagamma)*
                           gsl_ran_beta_pdf(gammaP,gammaC+1,2-gammaC);
                }
                temp = NAcp/DAcp;
                if (temp > 1.0) Acp = 1.0;else Acp = temp;
                temp = gsl_ran_flat(r,0.0,1.0);
                if (Acp > temp) gamma[h] = gammaP;else gamma[h] = gammaC;
            }else gamma[h] = gammaC;
            gammaC = gamma[h];

            // - 7 - Impute latent y^*
            for (i = 0; i < n; i++){
                if (compAlloc[i]==h){
                    nuSIxmm = 0.0;
                    for (k = 0; k < p; k++)
                        nuSIxmm += gsl_vector_get(nuSI,k)*(X[k*n+i]-muh[h][k]);
                    calcGLMLimits(Y,H,i,gammaC,&lower,&upper,family);
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

        // - 9 - P(delta_j = k_j)
        for (h = 0; h < nUpdated; h++){
            // View Sigmah^{-1}
            for (k = 0; k < p*p; k++)
                baseSigmahI[k] = SigmahI[h][k];
            SigmahIview = gsl_matrix_view_array(baseSigmahI,p,p);
            gammaC = gamma[h];
            for (i = 0; i < n; i++){
                if (pi[h] > u[i]){
                    nuSIxmm = 0.0;
                    for (k = 0; k < p; k++)
                        nuSIxmm += storenuSI[h][k]*(X[k*n+i]-muh[h][k]);
                    calcGLMLimits(Y,H,i,gammaC,&lower,&upper,family);
                    temp = gsl_cdf_ugaussian_P((upper-nuSIxmm)/sqrt(1-storenuSInu[h])) - gsl_cdf_ugaussian_P((lower-nuSIxmm)/sqrt(1-storenuSInu[h]));
                    for (k = 0; k < p; k++)
                        baseXmM[k] = X[k*n+i]-muh[h][k];
                    XmMview = gsl_vector_view_array(baseXmM,p);
                    temp2 = logMVNormalpdf3(p,&XmMview.vector,vecZero,&SigmahIview.matrix);
                    if (sampler == 1) pdjkj[h][i] = exp(log(temp)+temp2);
                    else if (sampler == 2) pdjkj[h][i] = exp(log(temp)+temp2+log(pi[h]));
                }
                else
                    pdjkj[h][i] = 0.0;
            }
        }
        for (h = nUpdated; h < ncomp; h++)
            for (i = 0; i < n; i++)
                pdjkj[h][i] = 0.0;

        // - 10 - Update allocation to
        s = gsl_ran_flat(r,1.0,100000);
        allocation(s,n,ncomp,pdjkj,compAlloc);

        // - 11 - Find nmembers
        findNmembers(n,ncomp,compAlloc,nmembers);

        // - 12a - Label switching
        s = gsl_ran_flat(r,1.0,100000);
        labelSwtchA(s,n,p,ncomp,Th1,Sigmah,SigmahI,nuh,muh,gamma,nmembers,compAlloc,pi);

        // - 12b - Label switching
        s = gsl_ran_flat(r,1.0,100000);
        labelSwtchB(s,n,p,ncomp,Th1,Sigmah,SigmahI,nuh,muh,gamma,nmembers,compAlloc,Vdp);

        // - 13 - Update concentration parameter
        s = gsl_ran_flat(r,1.0,100000);
        newalpha = updateAlpha(s,n,ncomp,Alphaa,Alphab,TruncAlpha,nmembers,alpha);
        alpha = newalpha;

        // - 14 - Predictions
        if (sw >= burn && npred > 0){
            // posterior weights
            for (i = 0; i < npred; i++) {newClusI[i] = 0; denomPred[i] = 0.0;}
            for (k = 0; k < maxy; k++) for (i = 0; i < npred; i++) StorePD[i][k] = 0.0;
            for (h = 0; h < ncomp; h++) for (i = 0; i < npred; i++) PredProb[i][h] = 0.0;
            if (sampler == 1){
                if (nUpdated < ncomp) upperL = nUpdated+1;
                else if (nUpdated == ncomp) upperL = nUpdated;
            }
            else if (sampler == 2) upperL = nUpdated;
            for (h = 0; h < upperL; h++){
                for (k = 0; k < p*p; k++)
                    baseSigmahI[k] = SigmahI[h][k];
                SigmahIview = gsl_matrix_view_array(baseSigmahI,p,p);
                for (i = 0; i < npred; i++){
                    for (k = 0; k < p; k++)
                        baseXmM[k] = Xpred[k*n+i]-muh[h][k];
                    XmMview = gsl_vector_view_array(baseXmM,p);
                    if (nmembers[h] > 0)
                        temp2 = logMVNormalpdf3(p,&XmMview.vector,vecZero,&SigmahIview.matrix);
                    else
                        temp2 = 0.0;
                    PredProb[i][h] = exp(temp2)*nmembers[h]/(alpha+n);
                    if (nmembers[h] == 0 && newClusI[i] == 0){
                        temp2 = logMVNormalpdf3(p,&XmMview.vector,vecZero,&SigmahIview.matrix);
                        PredProb[i][h] = exp(temp2)*alpha/(alpha+n);
                        newClusI[i] = 1;
                    }
                    if ((h == nUpdated) && (newClusI[i] == 0) && (nUpdated < ncomp)){
                        for (k = 0; k < p; k++) // N(0,1)
	                        gsl_vector_set(Z,k,gsl_ran_gaussian(r,1));
                        gsl_blas_dgemv(CblasNoTrans,1.0,SigmaMuHalf,Z,0.0,U1);  // sample from MV normal
                        gsl_vector_add(U1,MuMu);
                        for (k = 0; k < p; k++)
                            baseXmM[k] = Xpred[k*n+i]-gsl_vector_get(U1,k);
                        for (k = 0; k < p; k++) // N(0,1)
	                        gsl_vector_set(Z,k,gsl_ran_gaussian(r,1));
                        gsl_blas_dgemv(CblasNoTrans,1.0,SigmaNuHalf,Z,0.0,U1);  // sample from MV normal
                        gsl_vector_add(U1,MuNu);
                        rwish(s,p,Vdf,V,Th);
                        Inverse(p,Th);
                        nuMatrix = gsl_matrix_view_vector(U1,p,1);
                        gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,&nuMatrix.matrix,&nuMatrix.matrix,1.0,Th);
                        temp2 = logMVNormalpdf2(p,&XmMview.vector,vecZero,Th);
                        PredProb[i][h] = exp(temp2)*alpha/(alpha+n);
                    }
                    denomPred[i] += PredProb[i][h];
                }
            }
            // posterior predictive distribution
            for (i = 0; i < npred; i++){
                for (h = 0; h < ncomp; h++){
                    if (PredProb[i][h]>0){
                        gammaC = gamma[h];
                        nuSIxmm = 0.0;
                        for (k = 0; k < p; k++)
                            nuSIxmm += storenuSI[h][k]*(Xpred[k*n+i]-muh[h][k]);
                        for (k = 0; k < maxy; k++){
                            calcGLMLimitsPred(Hpred,k,i,gammaC,&lower,&upper,family);
                            temp = gsl_cdf_ugaussian_P((upper-nuSIxmm)/sqrt(1-storenuSInu[h]))-
                                   gsl_cdf_ugaussian_P((lower-nuSIxmm)/sqrt(1-storenuSInu[h]));
                            Disthi[h][k] = temp * PredProb[i][h]/denomPred[i];
                        }
                    }
                }
                for (k = 0; k < maxy; k++)
                    for (h = 0; h < ncomp; h++)
                        if (PredProb[i][h]>0) StorePD[i][k] += Disthi[h][k];
            }
            // posterior summaries
            // Mean
            for (i = 0; i < npred; i++)
                for (k = 0; k < maxy; k++)
                    meanReg[i] += (double) (StorePD[i][k] * k)/((sweeps-burn)*Hpred[i]);
            // Median
            for (i = 0; i < npred; i++){
                sumProb = 0.0;
                k=0;
                while (sumProb <= 0.5){
                    sumProb += StorePD[i][k];
                    k++;
                }
                medianReg[i] += (double) (k-1)/((sweeps-burn)*Hpred[i]);
            }
            // Q1
            for (i = 0; i < npred; i++){
                sumProb = 0.0;
                k=0;
                while (sumProb <= 0.25){
                    sumProb += StorePD[i][k];
                    k++;
                }
                q1Reg[i] += (double) (k-1)/((sweeps-burn)*Hpred[i]);
            }
            // Q3
            for (i = 0; i < npred; i++){
                sumProb = 0.0;
                k=0;
                while (sumProb <= 0.75){
                    sumProb += StorePD[i][k];
                    k++;
                }
                q3Reg[i] += (double) (k-1)/((sweeps-burn)*Hpred[i]);
            }
            // Mode
            for (i = 0; i < npred; i++){
                mode = gsl_stats_max_index(StorePD[i],1,maxy);
                modeReg[i] += (double) (mode)/((sweeps-burn)*Hpred[i]);
            }
        }

        // - 15 - Write to files
        if ((sw >= burn) && (WF == 1)){

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < (p*p); j++)
                    fprintf(out_file1, "%f ", Th1[h][j]);
                fprintf (out_file1, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < (p*p); j++)
                    fprintf(out_file2, "%f ", Sigmah[h][j]);
                fprintf (out_file2, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < (p*p); j++)
                    fprintf(out_file3, "%f ", SigmahI[h][j]);
                fprintf (out_file3, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < p; j++)
                    fprintf(out_file4, "%f ", nuh[h][j]);
                fprintf (out_file4, "\n");
            }

            for (h = 0; h < ncomp; h++){
                for (j = 0; j < p; j++)
                    fprintf(out_file5, "%f ", muh[h][j]);
                fprintf (out_file5, "\n");
            }

            for (h = 0; h < ncomp; h++)
                fprintf(out_file6, "%f \n", gamma[h]);

            fprintf(out_file7, "%f \n", alpha);

            for (i = 0; i < n; i++)
                fprintf(out_file8, "%i ", compAlloc[i]);
            fprintf (out_file8, "\n");

            for (h = 0; h < ncomp; h++)
                fprintf(out_file9, "%i ", nmembers[h]);
            fprintf (out_file9, "\n");

            fprintf(out_file10, "%i \n", nUpdated);

            for (i = 0; i < npred; i++){
                for (k = 0; k < maxy; k++)
                    fprintf(out_file11, "%f ", StorePD[i][k]);
                fprintf (out_file11, "\n");
            }
        }
        if (sw==(sweeps-1) && (!((sw+1) % 1000)==0)) Rprintf("%i %s \n",sw+1, "posterior samples...");
    }//end of sw
    //Free up random number
    gsl_rng_free (r);

    //Close files
    fclose(out_file1); fclose(out_file2); fclose(out_file3); fclose(out_file4); fclose(out_file5);
    fclose(out_file6); fclose(out_file7); fclose(out_file8); fclose(out_file9); fclose(out_file10);
    fclose(out_file11);

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
    //free(path); free(copypath);
}
