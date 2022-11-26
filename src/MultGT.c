/* multgT.c Bayesian Nonparametric Multivariate Regression 
 * Copyright (C) 2022  Georgios Papageorgiou, gpapageo@gmail.com
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
#include <gsl/gsl_multiset.h>
#include <gsl/gsl_sf_psi.h>
#include "cubature.h"
/*
#include "matalg.h"
#include "pdfs.h"
#include "sampling.h"
#include "other.functions.h"
#include "mathm.h" 
#include "spec.BCM.h"
#include "spec.T.h"
*/
extern void   computeStStar(double *Y, int *time, int N, int t, int p, gsl_matrix *StStar); 
extern void   ginv(int p, double tol, gsl_matrix *A);
extern double FisherTr(double r, int I);
extern double ScalcMult(int p, int m, int LG, double tol, double ceta, int Ngamma, double *Ytilde, double sigma2ij[m][p], double *X, int gamma[p][LG], gsl_matrix *Ri, gsl_matrix *St, double *qf2, double U[m][p], int mcm);
extern void   proposeBlockInd(unsigned long int s, int *vecInd, int L, int B, int BS, int *shufInd, double c, double d, int *vecIndP);
extern void   postMeanVarEta2(int p, int m, int LG, double tol, double ceta, int Ngamma, double *Ytilde, double sigma2ij[m][p], double *X, int gamma[p][LG], gsl_matrix *Ri, gsl_vector *MeanEta, gsl_matrix *varEta, double U[m][p], int mcm);
extern void   cSqRes2(int p, int m, int LG, int gamma[p][LG], int Ngamma, double *X, gsl_vector *MeanEta, double *Y, double *sqRes);
extern void   DeltaAlphaHatExp(int m, int p, int l, double tol, double LPV[m][p], double *sqRes, int *delta, int Ndelta, int start, int end, double *AllBases, double sigma2, double sigma2ij[m][p], double calpha, gsl_matrix *D, gsl_vector *alphaHat, double U[m][p], int mcm);
extern void   sampleMN(unsigned long int s, int p, gsl_vector *y, gsl_vector *mu, gsl_matrix *Sigma, double tol);
extern double logMVNormalpdf(int dim, gsl_vector *x, gsl_vector *mu, gsl_matrix *S, double tol);
extern void   rwish(unsigned long int s, int p, double n, gsl_matrix *scale, gsl_matrix *rw);
extern void   decomposeEtoDS(int nres, int nconf, gsl_matrix *E, gsl_matrix *D, gsl_matrix *S);
extern double logPropPdfDR(gsl_matrix *D, gsl_matrix *E, gsl_matrix *M, gsl_matrix *K, int p, double num1, double num2, double num3, double num4);
extern double det(int p, gsl_matrix *E);
extern void allocation(unsigned long int s, int n, int ncomp, double Prob[ncomp][n], int *compAlloc, int sw);
extern void findNmembers(int n, int ncomp, int *compAlloc, int *nmembers);
extern double updateAlpha(unsigned long int s, int n, int ncomp, double a, double b, double TruncAlpha, int *nmembers, double alpha);
extern double NormalQuadr(int p, int m, int LG, int Ngamma, double *Ytilde, double sigma2ij[m][p], double *X, int gamma[p][LG], gsl_matrix *Ri, double *beta, double U[m][p], int mcm);
extern void   DeltaPsiHat(int m, int p, int l, int LC, double tol, double omegaij[m][p], double U[m][p], double psi[p][LC], double cpsi, int start, int end, double *AllBases, int *delta, int Ndelta, gsl_vector *alphaHat, gsl_matrix *D);
extern int MultiNormalGammaPDF(unsigned dim, const double *u, void *parameters, unsigned fdim, double *fval);

extern int (*p_pcubature)(unsigned, integrand, void *, unsigned, const double *, const double *, size_t, double, double, error_norm, double *, double *);
extern int (*p_hcubature)(unsigned, integrand, void *, unsigned, const double *, const double *, size_t, double, double, error_norm, double *, double *);

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))   
#define MAX_PATH 300

void multgT(int *seed1, char **WorkingDir, 
            int *sbt1,
            double *Y, double *X, double *Z, double *W,
            int *LGDC, 
            double *blockSizeProb1, double *blockSizeProb2, double *blockSizeProb3, 
            int *maxBS, int *NGDC, int *vecLG, int *vecLD, int *vecLC,
            int *cusumVecLG, int *cusumVecLD, int *cusumVecLC, 
            double *tuneCa, double *tuneSigma, double *tuneCb, double *tuneAlpha, 
            double *tuneSigma2R, double *tuneR, double *tunePsi, double *tuneFi2, double *tuneCpsi,            
            double *pimu, double *cetaParams, double *pisigma, int *HNca, double *calphaParams,            
            int *HNcp, double *cpsiParams, int *HNfi, double *fiParams,            
            double *Rparams, int *HNszk, double *szkParams, double *piomega,             
            double *tau1, int *FT, double *dev, double *Res, int *H1, double *DPparams,            
            int *isDz, 
            int *cont, int *LASTgamma, int *LASTdelta, double *LASTalpha,  
            double *LASTsigma2zk, double *LASTR, double *LASTDE,
            double *LASTmuRsigma2R, 
            int *LASTsw, int *LASTWB,            
            double *LASTceta, double *LASTcalpha,             
            int *LASTcompAlloc, double *LASTalphaDP, 
            int *LASTksi, double *LASTpsifi2cpsi, int *isSin)         
{
    gsl_set_error_handler_off();
 
    // Random number generator initialization
    int seed = seed1[0]; //random number generator seed
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);
    unsigned long int s; //for calling random number generating functions

    // Specify directory
    char path_name[MAX_PATH + 1];

    // Open files
    FILE *out_file1, *out_file2, *out_file3, *out_file4, *out_file5, *out_file6,
         *out_file7, *out_file8, *out_file9, *out_file10, *out_file11, *out_file12, 
         *out_file13, *out_file14, *out_file15, *out_file16, *out_file17, *out_file18, 
         *out_file19;

    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.gamma.txt");
    out_file1 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.cbeta.txt");
    out_file2 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.delta.txt");
    out_file3 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.alpha.txt");
    out_file4 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.R.txt");
    out_file5 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.DE.txt");
    out_file6 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.muR.txt");
    out_file7 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.sigma2R.txt");
    out_file8 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.calpha.txt");
    out_file9 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.sigma2.txt");
    out_file10 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.beta.txt");
    out_file11 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.compAlloc.txt");
    out_file12 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.nmembers.txt");
    out_file13 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.deviance.txt");
    out_file14 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.DPconc.txt");
    out_file15 = fopen(path_name, "a");    
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.psi.txt");
    out_file16 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.ksi.txt");
    out_file17 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.phi2.txt");
    out_file18 = fopen(path_name, "a");    
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.cpsi.txt");
    out_file19 = fopen(path_name, "a");
     
    // Sweeps, burn-in period, thin
    int sweeps = sbt1[0]; //number of posterior samples
    int burn = sbt1[1]; // integer burn-in period
    int thin = sbt1[2]; // integer thin

    // Dimensions
    int m = sbt1[3];
    int p = sbt1[4]; 
    int mcm = sbt1[5];
    int d = p*(p-1)/2;
    if (d==0) d=1;
    int LG = LGDC[0];
    int LD = LGDC[1];
    int LC = LGDC[2];
    int NG = NGDC[0];
    int ND = NGDC[1];
    int NC = NGDC[2];
    int H = H1[0];//number of clusters
    int move; 
    //Dim Fix
    int MVLD = sbt1[6];
    
    //Tolerance level  
    double tol = 0.000000015;//0.00000001; //tolerance level for eigen values when inverting matrices
    
    //Declare variables for loops:
    int sw, h1, i, j, k, l;
     
    // Prior parameters
    //c.beta
    double alphaeta = cetaParams[0];
    double betaeta = cetaParams[1];
    //Rprintf("%s %f %f \n","c.beta: ",alphaeta,betaeta);
    
    //pi_mu, pi_sigma, pi_omega
    double cmu[p][NG];
    double dmu[p][NG];
    double csigma[p][ND];
    double dsigma[p][ND]; 
    double comega[p][NC];
    double domega[p][NC];     
	move = 0;
	for (k = 0; k < p; k++){
	    for (l = 0; l < NG; l++){
            cmu[k][l] = pimu[move++];
            dmu[k][l] = pimu[move++];
	    }    
    }
    move = 0;
    for (k = 0; k < p; k++){
        for (l = 0; l < ND; l++){
            csigma[k][l] = pisigma[move++];
            dsigma[k][l] = pisigma[move++];
		}
	}
	move = 0;
    for (k = 0; k < p; k++){
        for (l = 0; l < NC; l++){
            comega[k][l] = piomega[move++];
            domega[k][l] = piomega[move++];
		}
	}
	//Rprintf("%s %f %f %f %f %f %f %f %f \n","c-d mu: ",cmu[0][0],dmu[0][0],cmu[0][1],dmu[0][1],cmu[1][0],dmu[1][0],cmu[1][1],dmu[1][1]);
	//Rprintf("%s %f %f %f %f %f %f %f %f \n","c-d sg",csigma[0][0],dsigma[0][0],csigma[0][1],dsigma[0][1],csigma[1][0],dsigma[1][0],csigma[1][1],dsigma[1][1]);
    
    //c.alpha & sigma^2_{zero k}
    double phi2G;  //for half normal prior
    double alphaG, betaG; //for IG prior
    //varmuR & varsigma2R 
    double varmuR, varsigma2R;  
    varmuR = Rparams[1];
    varsigma2R = Rparams[2];       
    //Rprintf("%s %f %f \n","R's: ",varmuR,varsigma2R);
    // DP alpha ~ Gamma(a,b) I[alpha>c]
    double Alphaa = DPparams[0];  
    double Alphab = DPparams[1];
    double TruncAlpha = DPparams[2];
    
    //tau2 
    double tau2 = tau1[0]*tau1[0];

    // Declare quantities of interest
    double Rvec[d];
    double grt[d];
    double eta[p*(LG+1)];
    int gamma[p][LG];
    int Ngamma;   
    int delta[p][LD];
    int ksi[p][LC];
    double alpha[p][LD];
    double calpha[p];
    double ceta;
    double muR[H];
    double sigma2R;
    double sigma2zk[p];
    double sigma2ij[m][p];
    double psi[p][LC];
    double cpsi[p];
    double fi2[p];
    double omegaij[m][p];
    double alphaDP; //concentration parameter
    int compAlloc[d];
    double pdh[H][d];
    int nmembers[H];
        
    int *vecGamma[NG];
    for (j = 0; j < NG; j++)
        vecGamma[j] = malloc(sizeof(int) * vecLG[j]);
    int *vecDelta[ND];
    for (j = 0; j < ND; j++)
        vecDelta[j] = malloc(sizeof(int) * vecLD[j]);
    int *vecKsi[NC];
    for (j = 0; j < NC; j++)
        vecKsi[j] = malloc(sizeof(int) * vecLC[j]);
         
    // Declare other quantities
    double w[H]; //prior probabilities
    double Vdp[H-1]; //V[i] ~ Beta(1,alpha)
    double newAlphaDP; //store concentration parameter
    int nmb; //cumulative number of members
    double theta[d];
    double vizTheta[d];
    double LPM[m][p];
    double LPV[m][p];
    double LPVP[m][p];
    double LPS[m][p];
    double LPSP[m][p];
    double psiP[p][LC];
    double psiPK[LC];
    double alphaPD[LD];
    double BaseSubAlpha[MVLD]; 
    double sqResC[m*p];
    double sqResP[m*p];
    double residuals[m][p];
    double U[m][p];
    double logPropDRP, logPropDRC, priorLogR, logAcp, Acp, unifRV, QFC, QFD, detR,
           logMVNormC, logMVNormP, SPC, SPP, temp, sigma2RP, logLikP, logLikC, dev0;
    int NPJ, NCJ; 
    double cetahat, Sprime, SDprime, elPrime, elDPrime, Q2, cetaP, calphaP;
    double rtilde[m*p];
    double Ytilde[m*p];
    double YtildeP[m*p];
    double sumtheta, sstheta;
    double sigma2ijP[m][p];
    double sumlnj, Qa, Qb, Qc, mode, curv, UijP, p1, p2;
    double omegaijP[m][p];
    double val, err; //value and error from cubature function call
    double mplog = m*p*log(2*M_PI);
    double params[p*(p+5)/2]; 
    int compDev = sbt1[7];
    double lower[p];
    double upper[p];
    for (j = 0; j < p; j++){
        lower[j] = 0.0;
        upper[j] = 50.0;
	}
    
    //For selecting block size gamma_B or delta_B
    int block, blockSize;
    double nBlocks;
    int maxBSG = maxBS[0];
    double blockSizeProbG[maxBSG];  
    for (k = 0; k < maxBSG; k++)
        blockSizeProbG[k] = blockSizeProb1[k];
    unsigned int vecBSG[maxBSG];
    int maxBSD = maxBS[1];  
    double blockSizeProbD[maxBSD];
    for (k = 0; k < maxBSD; k++)
        blockSizeProbD[k] = blockSizeProb2[k];
    unsigned int vecBSD[maxBSD];
    int maxBSC = maxBS[2];
    double blockSizeProbC[maxBSC];
    for (k = 0; k < maxBSC; k++)
        blockSizeProbC[k] = blockSizeProb3[k];
    unsigned int vecBSC[maxBSC];
      
    //Gamma indeces: for shuffling
    int *indexG[NG];
    for (j = 0; j < NG; j++)
        indexG[j] = malloc(sizeof(int) * vecLG[j]);
    for (j = 0; j < NG; j++)
        for (k = 0; k < vecLG[j]; k++)
            indexG[j][k] = k;
            
    //Delta indeces: for shuffling
    int *indexD[ND];
    for (j = 0; j < ND; j++)
        indexD[j] = malloc(sizeof(int) * vecLD[j]);
    for (j = 0; j < ND; j++)
        for (k = 0; k < vecLD[j]; k++)
            indexD[j][k] = k;
            
    //Xi indeces: for shuffling
    int *indexX[NC];
    for (j = 0; j < NC; j++)
        indexX[j] = malloc(sizeof(int) * vecLC[j]);
    for (j = 0; j < NC; j++)
        for (k = 0; k < vecLC[j]; k++)
            indexX[j][k] = k;

    //Proposed gamma & delta
    int NgammaP;
    int *vecGammaP[NG];
    for (j = 0; j < NG; j++)
        vecGammaP[j] = malloc(sizeof(int) * vecLG[j]);        
    int *vecDeltaP[ND];
    for (j = 0; j < ND; j++)
        vecDeltaP[j] = malloc(sizeof(int) * vecLD[j]);    
    int *vecKsiP[NC];
    for (j = 0; j < NC; j++)
        vecKsiP[j] = malloc(sizeof(int) * vecLC[j]);
    int vecGammaSin[2];//for sinusoidals
    int vecGammaSinP[2];//for sinusoidals fix to pointers
        
    // Declare and allocate gsl vectors and matrices
    gsl_matrix *PSIt = gsl_matrix_alloc(p,p);
    gsl_matrix *CopyPSIt = gsl_matrix_alloc(p,p);
    gsl_matrix *EtC = gsl_matrix_alloc(p,p);
    gsl_matrix *DtC = gsl_matrix_alloc(p,p);
    gsl_matrix *RtC = gsl_matrix_alloc(p,p);
    gsl_matrix *RtCinv = gsl_matrix_alloc(p,p);
    gsl_matrix *EtP = gsl_matrix_alloc(p,p);
    gsl_matrix *DtP = gsl_matrix_calloc(p,p);
    gsl_matrix *RtP = gsl_matrix_alloc(p,p);
    gsl_matrix *RtPinv = gsl_matrix_alloc(p,p);
    gsl_matrix *A = gsl_matrix_alloc(d,d);
    gsl_matrix *D = gsl_matrix_alloc(MVLD,MVLD);
    gsl_matrix *varEta = gsl_matrix_alloc(p*(LG+1),p*(LG+1));
    gsl_vector *meanEta = gsl_vector_alloc(p*(LG+1));
    gsl_vector *alphaHat = gsl_vector_alloc(MVLD);
    gsl_vector *alphaP = gsl_vector_alloc(MVLD);
    gsl_matrix *St = gsl_matrix_alloc(p,p);
    gsl_matrix *St2 = gsl_matrix_alloc(p,p);
    gsl_matrix *StP = gsl_matrix_alloc(p,p);    
   
    // GSL Matrix and vector views
    gsl_matrix_view subD, subVarEta;
    gsl_vector_view vecTheta, Ar, subAlphaHat,subAlphaP, subMeanEta, vecEta;
     
    //Make adaptive
    int batchL = 50; //batch length
    double WB; //which batch
    
    double acceptR = 0.0; 
    double RtLL = p+2;
    double RtUL = 999999999999999;
    
    double acceptCeta = 0.0;
    double GLL = 0.01;
    double GUL = 100;
    
    double acceptSigma2R = 0.0;
    double F1LL = 0.001;//0.00000000001;
    double F1UL = 200;
    
    double acceptAlpha[p][ND];
    for (l = 0; l < p; l++)
	    for (j = 0; j < ND; j++)
	        acceptAlpha[l][j] = 0.0;	
    double HLL = 0.5;//0.00000000001;
    double HUL = 200;
    
    double acceptSZK[p];
    for (k = 0; k < p; k++)
	    acceptSZK[k] = 0.0;
    double SZKLL = 0.001;
    double SZKUL = 200;
    
    double acceptCa[p];
    for (k = 0; k < p; k++)
	    acceptCa[k] = 0.0;
    double FCaLL = 0.01;
    double FCaUL = 200;
    
    double acceptPsi[p][NC];
    for (l = 0; l < p; l++)
	    for (j = 0; j < NC; j++)
	        acceptPsi[l][j] = 0.0;
    double LLL = 0.5;//0.00000000001;
    double LUL = 200;
    
    double acceptFi2[p];
    for (k = 0; k < p; k++)
        acceptFi2[k] = 0.0;
    double Fi2LL = 0.001;
    double Fi2UL = 200;
    
    double acceptCpsi[p];
    for (k = 0; k < p; k++)
	    acceptCpsi[k] = 0.0;
    double CpLL = 0.01;
    double CpUL = 200;
    
    // Sampler initialization           
    int time[m*p];   
    for (i = 0; i < (m*p); i++)    
        time[i] = 0;
    
    // - 1 - Gamma      
    if (cont[0]==1){
        for (j = 0; j < p; j++)
            for (k = 0; k < LG; k++)
                gamma[j][k] = LASTgamma[k+j*LG];
    }else{    
        for (j = 0; j < p; j++)
            for (k = 0; k < LG; k++)
                gamma[j][k] = 1;
    }
          
    Ngamma = 0;
    for (j = 0; j < p; j++)
        for (k = 0; k < LG; k++)
            Ngamma += gamma[j][k];
    
    // - 2 - Delta & Alpha
    if (cont[0]==1){
        for (j = 0; j < p; j++){
            for (k = 0; k < LD; k++){
                delta[j][k] = LASTdelta[k+j*LD];
                alpha[j][k] = LASTalpha[k+j*LD]; 
		    }
	    }
    }else{
        for (j = 0; j < p; j++){
            for (k = 0; k < LD; k++){
                delta[j][k] = 0;
                alpha[j][k] = 0; 
		    }
	    }
	}
	
	// - 3 - Sigma2zk
    //if (cont[0]==1){
        for (j = 0; j < p; j++)
            sigma2zk[j] = LASTsigma2zk[j];
	//}else{    
    //    for (j = 0; j < p; j++)
    //        sigma2zk[j] = 1.0;
    //} 
    
    // - 4 - U
    for (i = 0; i < m; i++)
		for (j = 0; j < p; j++)
		    U[i][j] = gsl_ran_flat(r,0.9,1.1);  
    
    // - 5 - LPV, sigma2ij, Ytilde, St          
    for (i = 0; i < m; i++){
		for (j = 0; j < p; j++){
            LPV[i][j] = 0;
            for (k = 0; k < LD; k++)
                LPV[i][j] += alpha[j][k]*Z[m+k*m+i];         
            sigma2ij[i][j] = sigma2zk[j]*exp(LPV[i][j]);
            sigma2ijP[i][j] = sigma2ij[i][j];
            Ytilde[i*p+j] = Y[i*p+j]*sqrt(U[i][j]/sigma2ij[i][j]); 
            YtildeP[i*p+j] = Ytilde[i*p+j];
        }
    }
    computeStStar(Ytilde,time,m*p,0,p,St); 
    
    // - 6 - R related         
    //if (cont[0]==1){
        for (i = 0; i < p; i++){
		    for (j = 0; j < p; j++){
                gsl_matrix_set(RtC,i,j,LASTR[j+i*p]);
                gsl_matrix_set(DtC,i,j,LASTDE[j+i*p]);
                gsl_matrix_set(EtC,i,j,LASTDE[p*p+j+i*p]);
			}
		}
	//}else{
	//	  gsl_matrix_set_identity(RtC);
    //    gsl_matrix_set_identity(DtC);
    //    gsl_matrix_set_identity(EtC);	
	//}
        
    gsl_matrix_memcpy(RtCinv,RtC);
    ginv(p,tol,RtCinv);

    //if (cont[0]==0){
    //    for (j = 0; j < H; j++) 
    //        muR[j] = 0.0; 
    //    sigma2R = 1.0; 
    //}else{
    for (j = 0; j < H; j++) 
        muR[j] = LASTmuRsigma2R[j]; 
    sigma2R = LASTmuRsigma2R[H]; 
    //}
    
    move = 0;
    for (k = 0; k < (p-1); k++)
        for (l = (k+1); l < p; l++)
            Rvec[move++] = gsl_matrix_get(RtC,k,l);
    
    for (j = 0; j < d; j++)
        theta[j] = FisherTr(Rvec[j],FT[0]) + gsl_ran_gaussian(r,tau1[0]);
             
    // - 7 - c_eta    
    if (cont[0]==0){
        ceta = 1;
        cetahat = ceta; //starting value is the current value
        //Rprintf("%s %i %i %i %i \n","SPC & ceta ",p,m,LG,Ngamma);
        SPC = ScalcMult(p,m,LG,tol,ceta,Ngamma,Ytilde,sigma2ij,X,gamma,RtCinv,St,&Q2,U,mcm);
        elPrime = 99.9;
        while(elPrime > 0.000000001 || -elPrime > 0.000000001){
            Sprime = -Q2/pow(cetahat+1,2);
            SDprime = 2*Q2/pow(cetahat+1,3);
            elPrime = -0.5*(Ngamma+p)/(cetahat+1)-Sprime/2-(alphaeta+1)/cetahat+betaeta/pow(cetahat,2);
            elDPrime = 0.5*(Ngamma+p)/(pow(cetahat+1,2))-SDprime/2+(alphaeta+1)/pow(cetahat,2)-2*betaeta/pow(cetahat,3);
            cetahat -= elPrime / elDPrime;   
        }
        ceta = cetahat;
    }else{
		ceta = LASTceta[0]; 
	}
    
    //Rprintf("%s %f %f \n","SPC & ceta ",SPC,ceta);
        
    // - 8 - calpha
    move = 0;
    for (k = 0; k < p; k++){
        if (HNca[k]==0){
            alphaG = calphaParams[move++]; betaG = calphaParams[move++];
            if (alphaG > 1.0) calpha[k] = betaG/(alphaG-1);
            else calpha[k] = 110.0;
        }
        if (HNca[k]==1) {
			phi2G = calphaParams[move++];
			calpha[k] = sqrt(2*phi2G/M_PI);
		}
		if (cont[0]==1) calpha[k] = LASTcalpha[k];
    }
    
    // - 9 - Cluster allocation, alphaDP, Vdp & w
    for (i = 0; i < d; i++)
        for (h1 = 0; h1 < H; h1++)
            pdh[h1][i] = 1/((double)H);    
    s = gsl_ran_flat(r,1.0,100000);
    allocation(s,d,H,pdh,compAlloc,0);    
    if (cont[0]==1)
        for (i = 0; i < d; i++)
            compAlloc[i] = LASTcompAlloc[i];    
    findNmembers(d,H,compAlloc,nmembers);    
    alphaDP = 0.5;
    if (cont[0]==1)
        alphaDP = LASTalphaDP[0];
    nmb = 0;
    for (h1 = 0; h1 < H-1; h1++){
        nmb += nmembers[h1];
        Vdp[h1] = gsl_ran_beta(r, (double) nmembers[h1]+1.0, (double) d-nmb+alphaDP);
        if (h1 == 0)
            w[h1] = Vdp[h1];
        else 
            w[h1] = Vdp[h1] * w[h1-1] * (1-Vdp[h1-1]) / Vdp[h1-1];
    }
    w[H-1] = w[H-2] * (1-Vdp[H-2]) / Vdp[H-2];
    
    // - 10 - Residuals
    for (j = 0; j < p; j++)
        for (i = 0; i < m; i++)
            residuals[i][j] = Res[i+m*j];    
    
    // - 11 - Ksi & Psi
    if (cont[0]==1){
        for (j = 0; j < p; j++){
            for (k = 0; k < LC; k++){
                ksi[j][k] = LASTksi[k+j*LC];
                psi[j][k] = LASTpsifi2cpsi[k+j*LC];
		    }
	    }
    }else{    
        for (j = 0; j < p; j++){
            for (k = 0; k < LC; k++){
                ksi[j][k] = 0;
                psi[j][k] = 0;
		    }
	    }
	}

	// - 12 - fi2
	if (cont[0]==1){
        for (j = 0; j < p; j++)
            fi2[j] = LASTpsifi2cpsi[LC*p+j];
    }else{        
        for (j = 0; j < p; j++)
            fi2[j] = 2.0;
	}

    // - 13 - LPS, omegaij
    for (i = 0; i < m; i++){
		for (j = 0; j < p; j++){
            LPS[i][j] = 0;
            for (k = 0; k < LC; k++)
                LPS[i][j] += psi[j][k]*W[m+k*m+i];
            omegaij[i][j] = fi2[j]*exp(LPS[i][j]);
        }
    }
    
    // - 14 - cpsi
    move = 0;
    for (k = 0; k < p; k++){
        if (HNcp[k]==0){
            alphaG = cpsiParams[move++]; betaG = cpsiParams[move++];
            if (alphaG > 1.0) cpsi[k] = betaG/(alphaG-1);
            else cpsi[k] = 110.0;
        }
        if (HNcp[k]==1) {
			phi2G = cpsiParams[move++];
			cpsi[k] = sqrt(2*phi2G/M_PI);
		}
		if (cont[0]==1) cpsi[k] = LASTpsifi2cpsi[LC*p+p+j];
    }
                       
    //#############################################SAMPLER
    for (sw = 0; sw < sweeps; sw++){
        if (sw==0) Rprintf("%i %s \n",sw+1, "posterior sample...");
        if (((sw+1) % 1000)==0) Rprintf("%i %s \n",sw+1, "posterior samples...");
        
	    modf((sw+1)/batchL,&WB);
	    WB += LASTWB[0];
	    
	    // - 1 - Uij  
	    //Rprintf("%i %s \n",sw,"Uij");      
        for (i = 0; i < m; i++){
            for (j = 0; j < p; j++){
				sumlnj = 0;
	            for (l = 0; l < p; l++) 
	                if (l != j) sumlnj += residuals[i][l]*sqrt(U[i][l])*gsl_matrix_get(RtCinv,l,j)/sqrt(sigma2ij[i][l]);				
				Qa = 0.5*(omegaij[i][j]-1);
				Qb = -residuals[i][j]*sumlnj/sqrt(sigma2ij[i][j]);
				Qc = -0.5*(pow(residuals[i][j],2)*gsl_matrix_get(RtCinv,j,j)/sigma2ij[i][j]+omegaij[i][j]);
				if (omegaij[i][j] <= 1){
					p1 = 1; 
					p2 = -1/Qc; 
				}else{
					mode = 1/pow((-0.5*Qb + sqrt(0.25*pow(Qb,2)-4*Qa*Qc))/(2*Qa),2);
	                curv = -0.25*Qb*pow(mode,-1.5) - Qa*pow(mode,-2);
				    p1 = 1-curv*mode*mode; 
				    p2 = -1/(curv*mode);
				}				
				UijP = gsl_ran_gamma(r,p1,p2);                
                Acp = exp(Qc*(UijP-U[i][j])+Qb*(sqrt(UijP)-sqrt(U[i][j]))+Qa*(log(UijP)-log(U[i][j])))*
                      gsl_ran_gamma_pdf(U[i][j],p1,p2)/gsl_ran_gamma_pdf(UijP,p1,p2);
                unifRV = gsl_ran_flat(r,0.0,1.0);
	            if (Acp > unifRV){
	                U[i][j] = UijP;			    
	                Ytilde[i*p+j] = Y[i*p+j]*sqrt(U[i][j]/sigma2ij[i][j]);
	                YtildeP[i*p+j] = Ytilde[i*p+j]; 
				}
			    //fprintf(out_file17,"%s %i %i %i | %f %f %f | %f %f %f | %f %f %f | %f \n","Uij: ",
			    //                    sw,i,j,omegaij[i][j],sigma2ij[i][j],residuals[i][j],Qa,Qb,Qc,sumlnj,mode,curv,U[i][j]);	
			                
			}
		}	    
      
        // - 2 - gamma  
        //Rprintf("%i %s \n",sw,"gamma");				        
        computeStStar(Ytilde,time,m*p,0,p,St);                
		SPC = ScalcMult(p,m,LG,tol,ceta,Ngamma,Ytilde,sigma2ij,X,gamma,RtCinv,St,&Q2,U,mcm);		
		
		//Rprintf("%s %i %i %f %f %f \n","before gamma: ",sw,Ngamma,sigma2zk[0]*SPC,sigma2zk[0]*Q2,ceta);				
		
		for (l = 0; l < p; l++){		    
		    for (j = 0; j < NG; j++){	
				//Rprintf("%i %i %i %s \n",sw,l,j,"gamma");			        		        		        
		        for (k = 0; k < vecLG[j]; k++) 
		            vecGamma[j][k] = gamma[l][cusumVecLG[j]+k];			        		        		        
		        if (!isSin[j]){
		            gsl_ran_multinomial(r,maxBSG,1,blockSizeProbG,vecBSG);
                    blockSize = 0;
	                while(vecBSG[blockSize]==0) blockSize++;
                    blockSize += 1;
                    nBlocks = ceil((double)vecLG[j]/blockSize);
                    gsl_ran_shuffle(r,indexG[j],vecLG[j],sizeof(int));
                }
     		    else{
			        blockSize=1; 
			        nBlocks=1;
			    }	            	            
	            for (block = 0; block < nBlocks; block++){                    
                    s = gsl_ran_flat(r,1.0,100000);	
                    if (!isSin[j])                	                
	                    proposeBlockInd(s,vecGamma[j],vecLG[j],block,blockSize,indexG[j],cmu[l][j],dmu[l][j],vecGammaP[j]);
	                else{
	                    vecGammaSin[0]=vecGamma[j][0];
	                    proposeBlockInd(s,vecGammaSin,1,block,blockSize,indexG[j],cmu[l][j],dmu[l][j],vecGammaSinP); 	                
	                    vecGammaP[j][0]=vecGammaSinP[0];
	                    vecGammaP[j][1]=vecGammaSinP[0];
				    } 	                    
	                for (k = 0; k < vecLG[j]; k++)
	                    gamma[l][cusumVecLG[j]+k] = vecGammaP[j][k];	                    	            	            	                
	                NPJ = 0;  
                    NCJ = 0;
                    for (k = 0; k < vecLG[j]; k++){
                        NPJ += vecGammaP[j][k];
                        NCJ += vecGamma[j][k];
					} 					
					NgammaP = Ngamma - NCJ + NPJ;    						       	            	     	                	                	                
	                SPP = ScalcMult(p,m,LG,tol,ceta,NgammaP,Ytilde,sigma2ij,X,gamma,RtCinv,St,&Q2,U,mcm);                                        
                    Acp = exp((SPC-SPP)/2)*pow(ceta+1,0.5*(NCJ-NPJ));                    
                    //Rprintf("%s %i %i %i %i %f %f %f %f %i %i\n","gamma: ",sw,l,j,block,Acp,SPC,SPP,ceta,NCJ,NPJ);				
                    unifRV = gsl_ran_flat(r,0.0,1.0);                                
	                if (Acp > unifRV){
                        for (k = 0; k < vecLG[j]; k++)	                  
	                        vecGamma[j][k] = vecGammaP[j][k];
                        Ngamma = NgammaP; 
                        SPC = SPP;
                    } 
                    else{  
                        for (k = 0; k < vecLG[j]; k++)
	                        gamma[l][cusumVecLG[j]+k] = vecGamma[j][k];                        
				    }				    
                }  	    	   
	        }
	    }
	    //Rprintf("%s %i %i %f %f \n","after gamma: ",sw,Ngamma,sigma2zk[0]*SPC,sigma2zk[0]*Q2);
	    
	    // - 3 - Update delta and alpha
        //Rprintf("%i %s \n",sw,"delta & alpha"); 
        subMeanEta = gsl_vector_subvector(meanEta,0,Ngamma+p);
        subVarEta = gsl_matrix_submatrix(varEta,0,0,Ngamma+p,Ngamma+p);
        
        //Rprintf("%s %i %i %i %f %f %i \n","before alpha: ",p,m,LG,tol,ceta,Ngamma);	                        
        
        postMeanVarEta2(p,m,LG,tol,ceta,Ngamma,Ytilde,sigma2ij,X,gamma,RtCinv,&subMeanEta.vector,&subVarEta.matrix,U,mcm); 
        
        //puts("var");
        //print_matrix(&subVarEta.matrix);
          
        //for (k=0;k<(Ngamma+p);k++)
        //    Rprintf("%i %s %f \n",sw,"mean_eta:",gsl_vector_get(meanEta,k));
        
        cSqRes2(p,m,LG,gamma,Ngamma,X,&subMeanEta.vector,Y,sqResC);                
        
        //for (k=0;k<(m*p);k++)
            //Rprintf("%i %s %f \n",sw,"sqResC",sqResC[k]);                     
         
        //for (l = 0; l < p; l++) 
        //    for (j = 0; j < LG; j++)
		//	    Rprintf("%i ",gamma[l][j]);
		//Rprintf("%i \n",Ngamma); 
         
        for (l = 0; l < p; l++){ 			            
            for (j = 0; j < ND; j++){
				//Rprintf("%i %i %i %s \n",sw,l,j,"delta & alpha");		        
		        if ((sw % batchL)==0 && WB > 0){ 
	                if (acceptAlpha[l][j] > 0.25 && tuneAlpha[ND*l+j] < HUL) tuneAlpha[ND*l+j] += MIN(0.01,1/sqrt(WB)) * tuneAlpha[ND*l+j]; 
	                if (acceptAlpha[l][j] <= 0.25 && tuneAlpha[ND*l+j] > HLL) tuneAlpha[ND*l+j] -= MIN(0.01,1/sqrt(WB)) * tuneAlpha[ND*l+j];
	                acceptAlpha[l][j] = 0.0;   
	                if (tuneAlpha[ND*l+j] < HLL) tuneAlpha[ND*l+j] = HLL;
	                if (tuneAlpha[ND*l+j] > HUL) tuneAlpha[ND*l+j] = HUL;
                }                            
		        
		        for (k = 0; k < vecLD[j]; k++) 
                    vecDelta[j][k] = delta[l][cusumVecLD[j] + k];                   
		        
		        gsl_ran_multinomial(r,maxBSD,1,blockSizeProbD,vecBSD);
                blockSize = 0;
                while(vecBSD[blockSize]==0) blockSize++;
                blockSize += 1; 
                nBlocks = ceil((double)vecLD[j]/blockSize);
                
                gsl_ran_shuffle(r,indexD[j],vecLD[j],sizeof(int));                                 
                
                for (block = 0; block < nBlocks; block++){
                    s = gsl_ran_flat(r,1.0,100000);                    
                    proposeBlockInd(s,vecDelta[j],vecLD[j],block,blockSize,indexD[j],csigma[l][j],dsigma[l][j],vecDeltaP[j]);
                    
                    //vecDeltaP[j][0]=1;vecDeltaP[j][1]=0;vecDeltaP[j][2]=1;
                    //vecDeltaP[j][3]=1;vecDeltaP[j][4]=1;vecDeltaP[j][5]=0;
                    //vecDeltaP[j][6]=0;vecDeltaP[j][7]=0;vecDeltaP[j][8]=0;
                    //vecDeltaP[j][9]=0; vecDeltaP[j][10]=0;               
                    
                    NPJ = 0;  
                    NCJ = 0;                       
                    for (k = 0; k < vecLD[j]; k++){
                        NPJ += vecDeltaP[j][k];
                        NCJ += vecDelta[j][k];
                    }                      
                    if (NPJ > 0){                 
                        subAlphaHat = gsl_vector_subvector(alphaHat,0,NPJ);
                        subD = gsl_matrix_submatrix(D,0,0,NPJ,NPJ);                        
                        if (isDz[j] == 0){
                            DeltaAlphaHatExp(m,p,l,tol,LPV,sqResC,vecDeltaP[j],NPJ,cusumVecLD[j],cusumVecLD[j+1], 
                                             Z,sigma2zk[l],sigma2ij,calpha[l],&subD.matrix,&subAlphaHat.vector,U,mcm);  
                        }else{ 
                            gsl_vector_set(&subAlphaHat.vector,0,alpha[l][cusumVecLD[j]]);
                            gsl_matrix_set(&subD.matrix,0,0,1);
					    }                      
                        subAlphaP = gsl_vector_subvector(alphaP,0,NPJ);                        
                        gsl_matrix_scale(&subD.matrix,tuneAlpha[ND*l+j]);                        
                        s = gsl_ran_flat(r,1.0,100000);
                        sampleMN(s,NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);                        
                        
                        //for (k=0;k<NPJ;k++)
                        //    gsl_vector_set(&subAlphaP.vector,k,gsl_vector_get(&subAlphaHat.vector,k));
                        
                        logMVNormP = logMVNormalpdf(NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);                                               
                        
                        //print_matrix(&subD.matrix);
                        //for (k=0;k<NPJ;k++)
                        //    Rprintf("%s %i %i %i %i %i %f %f\n","hat and prop: ",
                        //        sw,l,j,block,NPJ,gsl_vector_get(alphaHat,k),gsl_vector_get(alphaP,k));
                    }else logMVNormP = 0.0;                                                                                     
                    move = 0;
                    for (k = 0; k < vecLD[j]; k++){
                        if (vecDeltaP[j][k]==1) {alphaPD[cusumVecLD[j]+k] = gsl_vector_get(&subAlphaP.vector,move++);
						//	Rprintf("%s %i %i %f \n","from if 1: ",k,vecDeltaP[j][k],alphaPD[cusumVecLD[j]+k]); 
					    }else {alphaPD[cusumVecLD[j]+k] = 0;
							   //Rprintf("%s %i %i %f \n","from if 0: ",k,vecDeltaP[j][k],alphaPD[cusumVecLD[j]+k]);
							   }
                        //if (sw==0 && l==0 && j==0 && block==0) 
                        //Rprintf("%s %i %i %i %i | %f %i %i %i %i %i | %f %f %f\n","alpha P: ",
                        //sw,l,j,block,nBlocks,vecLD[j],cusumVecLD[j],k,NPJ,vecDeltaP[j][k],
                        //gsl_vector_get(alphaHat,k),alphaPD[cusumVecLD[j]+k],gsl_vector_get(alphaP,k));
                        //if (k < NPJ) Rprintf("%f \n",gsl_vector_get(&subAlphaP.vector,k));
                    }                                        
                    for (i = 0; i < m; i++){
                        LPVP[i][l] = LPV[i][l];
                        for (k = cusumVecLD[j]; k < cusumVecLD[j+1]; k++){
                            LPVP[i][l] += (alphaPD[k]-alpha[l][k])*Z[m+k*m+i]; 
                            //if (1==0 && sw==0 && l==0 && j==0 && block==0) 
                            //Rprintf("%s %i %i %i %i %f %f %f %f \n","new:", sw,l,j,block,alphaPD[k],alpha[l][k],Z[m+k*m+i],LPVP[i][l]);
						}
                        sigma2ijP[i][l] = sigma2zk[l]*exp(LPVP[i][l]);
                        YtildeP[i*p+l] = Y[i*p+l]*sqrt(U[i][l]/sigma2ijP[i][l]);                                                                            
                    }               
                    computeStStar(YtildeP,time,m*p,0,p,StP);  
                    
                    //if (1==0){ 						    
					//    for (i = 0; i < m; i++){ 
					//	    for (k = 0; k < p; k++)
					//	        Rprintf("%s %i %i %f %f %f %f \n","s-t:",
					//	             i,k,sigma2zk[l],LPVP[i][k],sigma2ijP[i][k],YtildeP[i*p+k]); 				
					//	}
					//}
					//print_matrix(StP);                                           
                    					
                    SPP = ScalcMult(p,m,LG,tol,ceta,Ngamma,YtildeP,sigma2ijP,X,gamma,RtCinv,StP,&Q2,U,mcm);                    
                    //if (1==0 && sw==3 && l==1 && j==0 && block==0){ 
                    //    Rprintf("%s %i %i %i %i %f %i %i %i %f %f %i \n","spp int: ",sw,l,j,block,SPP,p,m,LG,tol,ceta,Ngamma);
                    //    for (k = 0; k < (m*p); k++)
                    //       Rprintf("%s %f %f ","tilde: ",Ytilde[k],YtildeP[k]);
                    //    for (i = 0; i < m; i++)
                    //    Rprintf("%s %f %f %f %f ","sig lpv: ",sigma2ij[0][l],sigma2ijP[0][l],LPV[0][l],LPVP[0][l]);
                    //}
                    QFD = 0.0;
                    for (k = cusumVecLD[j]; k < cusumVecLD[j+1]; k++) 
                        QFD += pow(alphaPD[k],2)-pow(alpha[l][k],2);                    
                    detR = 0.0;                    
                    for (i = 0; i < m; i++)
                        detR += LPV[i][l] - LPVP[i][l];
                    detR *= 0.5;     	               	                
	                //probability of reverse direction 	                
	                postMeanVarEta2(p,m,LG,tol,ceta,Ngamma,YtildeP,sigma2ijP,X,gamma,RtCinv,&subMeanEta.vector,&subVarEta.matrix,U,mcm);
                    cSqRes2(p,m,LG,gamma,Ngamma,X,&subMeanEta.vector,Y,sqResP);                    
                    if (NCJ > 0){                
	                    subAlphaHat = gsl_vector_subvector(alphaHat,0,NCJ);
                        subD = gsl_matrix_submatrix(D,0,0,NCJ,NCJ);                        
                        if (isDz[j] == 0){ 
                            DeltaAlphaHatExp(m,p,l,tol,LPVP,sqResP,vecDelta[j],NCJ,cusumVecLD[j],cusumVecLD[j+1], 
                                             Z,sigma2zk[l],sigma2ijP,calpha[l],&subD.matrix,&subAlphaHat.vector,U,mcm);                        
                        }else{
                            gsl_vector_set(&subAlphaHat.vector,0,alphaPD[cusumVecLD[j]]);
                            gsl_matrix_set(&subD.matrix,0,0,1);
					    }  
                        gsl_matrix_scale(&subD.matrix,tuneAlpha[ND*l+j]);                        
                        move=0;                        
                        for (k = 0; k < vecLD[j]; k++)
                            if (vecDelta[j][k]==1) BaseSubAlpha[move++] = alpha[l][cusumVecLD[j]+k]; 
                        subAlphaP = gsl_vector_view_array(BaseSubAlpha,NCJ);
                        logMVNormC = logMVNormalpdf(NCJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);                                 
	                }else logMVNormC = 0.0;	            	                	                
                    //logAcp = detR+(-SPP+SPC)/2+logMVNormC-logMVNormP+QFD/(2*calpha[l]);
                    
                    Acp = exp(detR+(-SPP+SPC)/2+logMVNormC-logMVNormP-QFD/(2*calpha[l]))*
                              pow(2*M_PI*calpha[l],0.5*(NCJ-NPJ));                                                                                       
                    
                    unifRV = gsl_ran_flat(r,0.0,1.0);           	                	                
	                
	                //Rprintf("%s %i %i %i %i | %f %f | %f %f %f %f | %f %f %f %f \n","delta: ",
                    //sw,l,j,block,Acp,unifRV,SPC,SPP,logMVNormC,logMVNormP,QFD,calpha[l],detR,pow(2*M_PI*calpha[l],0.5*(NCJ-NPJ)));
	                
	                
	                if (Acp > unifRV){	      	   	    
                        if (NPJ > 0) acceptAlpha[l][j] += 1/((double)batchL);
                        for (k = 0; k < vecLD[j]; k++){
                            delta[l][cusumVecLD[j]+k] = vecDeltaP[j][k];
                            vecDelta[j][k] = vecDeltaP[j][k];
                            alpha[l][cusumVecLD[j]+k] = alphaPD[cusumVecLD[j]+k];                  
                        }            
                        for (i = 0; i < m; i++){
                             LPV[i][l] = LPVP[i][l];
                             sigma2ij[i][l] = sigma2ijP[i][l]; 
                             Ytilde[i*p+l] = YtildeP[i*p+l];
				        }
				        for (i = 0; i < (m*p); i++)                            
                            sqResC[i] = sqResP[i];			
                        SPC = SPP;
                        gsl_matrix_memcpy(St,StP);	     
                    }
                    else{
				    	if (j == (ND-1)){
				    	    for (i = 0; i < m; i++){
                                LPVP[i][l] = LPV[i][l];
                                YtildeP[i*p+l] = Ytilde[i*p+l];
                                sigma2ijP[i][l] = sigma2ij[i][l];
							}
						}
				    }                      
                }	    
	        }
	    }
	    //Rprintf("%s %i %i %i %f %f %f %f \n","after alpha: ",p,m,LG,ceta,sigma2zk[0]*SPC,LPV[0][0],Ytilde[0]*sqrt(sigma2zk[0]));
	    
        // - 4 - Update xi and psi
        //Rprintf("%i %s \n",sw,"xi & psi");
        for (l = 0; l < p; l++){
            for (j = 0; j < NC; j++){
		        if ((sw % batchL)==0 && WB > 0){
	                if (acceptPsi[l][j] > 0.25 && tunePsi[NC*l+j] < LUL) tunePsi[NC*l+j] += MIN(0.01,1/sqrt(WB)) * tunePsi[NC*l+j];
	                if (acceptPsi[l][j] <= 0.25 && tunePsi[NC*l+j] > LLL) tunePsi[NC*l+j] -= MIN(0.01,1/sqrt(WB)) * tunePsi[NC*l+j];
	                acceptPsi[l][j] = 0.0;
	                if (tunePsi[NC*l+j] < LLL) tunePsi[NC*l+j] = LLL;
	                if (tunePsi[NC*l+j] > LUL) tunePsi[NC*l+j] = LUL;
                }
		        for (k = 0; k < vecLC[j]; k++)
                    vecKsi[j][k] = ksi[l][cusumVecLC[j] + k];
		        gsl_ran_multinomial(r,maxBSC,1,blockSizeProbC,vecBSC);
                blockSize = 0;
                while(vecBSC[blockSize]==0) blockSize++;
                blockSize += 1;
                nBlocks = ceil((double)vecLC[j]/blockSize);
                gsl_ran_shuffle(r,indexX[j],vecLC[j],sizeof(int));
                for (block = 0; block < nBlocks; block++){
                    s = gsl_ran_flat(r,1.0,100000);
                    proposeBlockInd(s,vecKsi[j],vecLC[j],block,blockSize,indexX[j],comega[l][j],domega[l][j],vecKsiP[j]);
                    NPJ = 0;
                    NCJ = 0;
                    for (k = 0; k < vecLC[j]; k++){
                        NPJ += vecKsiP[j][k];
                        NCJ += vecKsi[j][k];
                    }
                    if (NPJ > 0){
                        subAlphaHat = gsl_vector_subvector(alphaHat,0,NPJ);
                        subD = gsl_matrix_submatrix(D,0,0,NPJ,NPJ);
                        DeltaPsiHat(m,p,l,LC,tol,omegaij,U,psi,cpsi[l],cusumVecLC[j],cusumVecLC[j+1],W,vecKsiP[j],NPJ, 
                                    &subAlphaHat.vector, &subD.matrix);
                        subAlphaP = gsl_vector_subvector(alphaP,0,NPJ);
                        gsl_matrix_scale(&subD.matrix,tunePsi[NC*l+j]);
                        s = gsl_ran_flat(r,1.0,100000);
                        sampleMN(s,NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);
                        logMVNormP = logMVNormalpdf(NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);
                    }else logMVNormP = 0.0;
                    move = 0;
                    for (k = 0; k < vecLC[j]; k++){
                        if (vecKsiP[j][k]==1) {psiPK[cusumVecLC[j]+k] = gsl_vector_get(&subAlphaP.vector,move++);						
					    }else{psiPK[cusumVecLC[j]+k] = 0;}
                        psiP[l][cusumVecLC[j]+k] = psiPK[cusumVecLC[j]+k];  
                    }
                    for (i = 0; i < m; i++){
                        LPSP[i][l] = LPS[i][l];
                        for (k = cusumVecLC[j]; k < cusumVecLC[j+1]; k++)
                            LPSP[i][l] += (psiPK[k]-psi[l][k])*W[m+k*m+i];						
                        omegaijP[i][l] = fi2[l]*exp(LPSP[i][l]);                        
                    }                    
                    QFD = 0.0;
                    for (k = cusumVecLC[j]; k < cusumVecLC[j+1]; k++)
                        QFD += pow(psiPK[k],2)-pow(psi[l][k],2);
                    detR = 1.0;
                    for (i = 0; i < m; i++)
                        detR *= gsl_ran_gamma_pdf(U[i][l],omegaijP[i][l]/2,2/omegaijP[i][l])/gsl_ran_gamma_pdf(U[i][l],omegaij[i][l]/2,2/omegaij[i][l]);
	                //probability of reverse direction
                    if (NCJ > 0){
                        subAlphaHat = gsl_vector_subvector(alphaHat,0,NCJ);
                        subD = gsl_matrix_submatrix(D,0,0,NCJ,NCJ);
                        DeltaPsiHat(m,p,l,LC,tol,omegaijP,U,psiP,cpsi[l],cusumVecLC[j],cusumVecLC[j+1],W,vecKsi[j],NCJ, 
                                    &subAlphaHat.vector, &subD.matrix);                                              
                        gsl_matrix_scale(&subD.matrix,tunePsi[NC*l+j]);
                        move=0;
                        for (k = 0; k < vecLC[j]; k++)
                            if (vecKsi[j][k]==1) BaseSubAlpha[move++] = psi[l][cusumVecLC[j]+k];
                        subAlphaP = gsl_vector_view_array(BaseSubAlpha,NCJ);
                        logMVNormC = logMVNormalpdf(NCJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);
	                }else logMVNormC = 0.0;
                    Acp = detR * pow(2*M_PI*cpsi[l],0.5*(NCJ-NPJ)) * exp(-QFD/(2*cpsi[l]) + logMVNormC - logMVNormP);
                    unifRV = gsl_ran_flat(r,0.0,1.0);
                    
                    //if (sw > 2000)for (k = cusumVecLC[j]; k < cusumVecLC[j+1]; k++)
                    //    fprintf(out_file17,"%s %i %i %i %i | %f %f \n","p vs c: ",sw,l,j,block,psiPK[k],psi[l][k]);
	                //if (sw > 2000) fprintf(out_file17,"%s %i %i %i %i | %f %f | %f %f %f %f %f %f \n",
                    //"Acp: ",sw,l,j,block,
                    //Acp,unifRV,
                    //detR,pow(2*M_PI*cpsi[l],0.5*(NCJ-NPJ)),QFD,cpsi[l],logMVNormC,logMVNormP);
                    
	                if (Acp > unifRV){
                        if (NPJ > 0) acceptPsi[l][j] += 1/((double)batchL);
                        for (k = 0; k < vecLC[j]; k++){                            
                            ksi[l][cusumVecLC[j]+k] = vecKsiP[j][k];
                            vecKsi[j][k] = vecKsiP[j][k];
                            psi[l][cusumVecLC[j]+k] = psiPK[cusumVecLC[j]+k];
                        }
                        for (i = 0; i < m; i++){
                             LPS[i][l] = LPSP[i][l];
                             omegaij[i][l] = omegaijP[i][l];                             
				        }				                                                        
                    }                    
                }
	        }
	    }	    
	    
	    // - 5 - sigma2_k, k=1,...,p
        //Rprintf("%i %s \n",sw,"sigma2_k");
		move = 0;
		for (k = 0; k < p; k++){		
            if ((sw % batchL)==0 && WB > 0){ 
	            if (acceptSZK[k] > 0.25 && tuneSigma[k] < SZKUL) tuneSigma[k] += MIN(0.01,1/sqrt(WB)) * tuneSigma[k]; 
	            if (acceptSZK[k] <= 0.25 && tuneSigma[k] > SZKLL) tuneSigma[k] -= MIN(0.01,1/sqrt(WB)) * tuneSigma[k];
	            acceptSZK[k] = 0.0;   
	            if (tuneSigma[k] < SZKLL) tuneSigma[k] = SZKLL;
	            if (tuneSigma[k] > SZKUL) tuneSigma[k] = SZKUL;
            }                        
            sigma2RP = sigma2zk[k] + gsl_ran_gaussian(r,sqrt(tuneSigma[k]));
            //sigma2RP = 0.3;
            while (sigma2RP <= 0) sigma2RP = sigma2zk[k] + gsl_ran_gaussian(r,sqrt(tuneSigma[k]));                        
            
            for (i = 0; i < m; i++){                
                sigma2ijP[i][k] = sigma2RP*exp(LPV[i][k]);
                YtildeP[i*p+k] = Y[i*p+k]*sqrt(U[i][k]/sigma2ijP[i][k]);
                //Rprintf("%i %i | %f %f %f | %i | %f %f %f \n",
                //k,0,Y[i*p+0],sigma2ijP[i][0],YtildeP[i*p+0],1,Y[i*p+1],sigma2ijP[i][1],YtildeP[i*p+1]);                
            }               
            computeStStar(YtildeP,time,m*p,0,p,StP);                          
            //print_matrix(StP);
            SPP = ScalcMult(p,m,LG,tol,ceta,Ngamma,YtildeP,sigma2ijP,X,gamma,RtCinv,StP,&Q2,U,mcm);                    
            //detR = 0.0;                    
            //for (i = 0; i < m; i++)
            //    detR += log(sigma2ij[i][k]/sigma2ijP[i][k]);
            detR = 0.5*m*log(sigma2zk[k]/sigma2RP);     	               	                            
            //Prior Ratio
            if (HNszk[k]==1){
                phi2G = szkParams[move++];
                priorLogR = 0.5*(sigma2zk[k]-sigma2RP)/phi2G;
			}else{
				alphaG = szkParams[move++];
				betaG = szkParams[move++];
                priorLogR = (alphaG+1)*log(sigma2zk[k]/sigma2RP) + betaG*(1/sigma2zk[k]-1/sigma2RP);    
			}                         
            Acp = exp(detR+(SPC-SPP)/2+priorLogR);
            unifRV = gsl_ran_flat(r,0.0,1.0);            
                                    	        
	        //fprintf(out_file16,"%s %i %i %f %f %f %f %f %f %f %f %f \n","sigma2: ",sw,k,Acp,sigma2zk[k],sigma2RP,
	        //detR+(SPC-SPP)/2+priorLogR,priorLogR,(SPC-SPP)/2,detR,SPC,SPP);            
            
            if (Acp > unifRV){
                sigma2zk[k] = sigma2RP;             
	            acceptSZK[k] += 1/((double)batchL);  
	            gsl_matrix_memcpy(St,StP);
	            for (i = 0; i < m; i++){
                    sigma2ij[i][k] = sigma2ijP[i][k];
                    Ytilde[i*p+k] = YtildeP[i*p+k];
                }
                SPC = SPP;                
            }
            else{
			    for (i = 0; i < m; i++){
                    sigma2ijP[i][k] = sigma2ij[i][k];                    
                    YtildeP[i*p+k] = Ytilde[i*p+k];
                }
			}
	    } 
	    
	    // - 6 - fi2_k, k=1,...,p
        //Rprintf("%i %s \n",sw,"phi2_k");
		move = 0;
		for (k = 0; k < p; k++){
            if ((sw % batchL)==0 && WB > 0){
	            if (acceptFi2[k] > 0.25 && tuneFi2[k] < Fi2UL) tuneFi2[k] += MIN(0.01,1/sqrt(WB)) * tuneFi2[k];
	            if (acceptFi2[k] <= 0.25 && tuneFi2[k] > Fi2LL) tuneFi2[k] -= MIN(0.01,1/sqrt(WB)) * tuneFi2[k];
	            acceptFi2[k] = 0.0;
	            if (tuneFi2[k] < Fi2LL) tuneFi2[k] = Fi2LL;
	            if (tuneFi2[k] > Fi2UL) tuneFi2[k] = Fi2UL;
            }
            sigma2RP = fi2[k] + gsl_ran_gaussian(r,sqrt(tuneFi2[k]));           
            while (sigma2RP <= 0) sigma2RP = fi2[k] + gsl_ran_gaussian(r,sqrt(tuneFi2[k]));
            for (i = 0; i < m; i++)
                omegaijP[i][k] = sigma2RP*exp(LPS[i][k]);
             detR = 1.0;
             for (i = 0; i < m; i++)
                 detR *= gsl_ran_gamma_pdf(U[i][k],omegaijP[i][k]/2,2/omegaijP[i][k])/gsl_ran_gamma_pdf(U[i][k],omegaij[i][k]/2,2/omegaij[i][k]);                           
            //Prior Ratio
            if (HNfi[k]==1){
                phi2G = fiParams[move++];
                priorLogR = 0.5*(fi2[k]-sigma2RP)/phi2G;
			}else{
				alphaG = fiParams[move++];
				betaG = fiParams[move++];
                priorLogR = (alphaG+1)*log(fi2[k]/sigma2RP) + betaG*(1/fi2[k]-1/sigma2RP);
			}
            Acp = detR*exp(priorLogR);
            unifRV = gsl_ran_flat(r,0.0,1.0);
	        //fprintf(out_file12,"%s %i %i %f %f %f %f %f %f \n","fi2: ",sw,k,Acp,fi2[k],sigma2RP,detR,priorLogR,detR);
            if (Acp > unifRV){
                fi2[k] = sigma2RP;
	            acceptFi2[k] += 1/((double)batchL);
	            for (i = 0; i < m; i++)
                    omegaij[i][k] = omegaijP[i][k];
            }
            else{//THIS IS PROBABLY NOT NEEDED
			    for (i = 0; i < m; i++)
                    omegaijP[i][k] = omegaij[i][k];                                   
			}
	    }   
	    
	    // - 7 - c_beta 
	    //Rprintf("%i %s \n",sw,"c_beta");
        if ((sw % batchL)==0 && WB > 0){ 
	        if (acceptCeta > 0.25 && tuneCb[0] < GUL) tuneCb[0] += MIN(0.01,1/sqrt(WB)) * tuneCb[0]; 
	        if (acceptCeta <= 0.25 && tuneCb[0] > GLL) tuneCb[0] -= MIN(0.01,1/sqrt(WB)) * tuneCb[0];
	        acceptCeta = 0.0;   
	        if (tuneCb[0] < GLL) tuneCb[0] = GLL;
	        if (tuneCb[0] > GUL) tuneCb[0] = GUL;
        }     
        cetahat = 1;        
        SPC = ScalcMult(p,m,LG,tol,ceta,Ngamma,Ytilde,sigma2ij,X,gamma,RtCinv,St,&Q2,U,mcm);                       
        elPrime = 99.9;        
        while(elPrime > 0.000000001 || -elPrime > 0.000000001){
            Sprime = -Q2/pow(cetahat+1,2);
            SDprime = 2*Q2/pow(cetahat+1,3);
            elPrime = -0.5*(Ngamma+p)/(cetahat+1)-Sprime/2-(alphaeta+1)/cetahat+betaeta/pow(cetahat,2);
            elDPrime = 0.5*(Ngamma+p)/(pow(cetahat+1,2))-SDprime/2+(alphaeta+1)/pow(cetahat,2)-2*betaeta/pow(cetahat,3);
            cetahat -= elPrime / elDPrime;
        }        
        
//Rprintf("%s %f %f %f %f %f %f %f \n","ceta 1: ",cetahat,elPrime,SPC*sigma2zk[0],Sprime*sigma2zk[0],SDprime*sigma2zk[0],Q2*sigma2zk[0],sigma2zk[0]);        
         
        cetaP = cetahat + gsl_ran_gaussian(r,sqrt(-tuneCb[0]/elDPrime));	    
        //cetaP = 545;
	    while(cetaP < 0) cetaP = cetahat + gsl_ran_gaussian(r,sqrt(-tuneCb[0]/elDPrime));
	    SPP = ScalcMult(p,m,LG,tol,cetaP,Ngamma,Ytilde,sigma2ij,X,gamma,RtCinv,St,&Q2,U,mcm);
	    Acp = exp(-0.5*(Ngamma+p)*(log(cetaP+1)-log(ceta+1))+(-SPP+SPC)/2-(alphaeta+1)*(log(cetaP)-log(ceta))
              + betaeta*(1/ceta-1/cetaP))*
                gsl_ran_gaussian_pdf(ceta-cetahat,sqrt(-tuneCb[0]/elDPrime))/
                gsl_ran_gaussian_pdf(cetaP-cetahat,sqrt(-tuneCb[0]/elDPrime));
	    unifRV = gsl_ran_flat(r,0.0,1.0);
        if (Acp > unifRV){
            ceta = cetaP;
            SPC = SPP;
	        acceptCeta += 1/((double)batchL);  
	    }
	    
//Rprintf("%s %i %f %f %f %f %f %f %f %f %f %f %i %f %f\n","ceta 2: ",
//sw,Acp,unifRV,sigma2zk[0]*SPC,sigma2zk[0]*SPP,ceta,cetahat,cetaP,sigma2zk[0]*Q2,alphaeta,betaeta,Ngamma,tuneCb[0],elDPrime);	                        

        // - 8 - c_alpha j
        //Rprintf("%i %s \n",sw,"c_alpha j");        
        move = 0;
        for (k = 0; k < p; k++){
            NPJ = 0;
            QFC = 0;
            for (l = 0; l < LD; l++){
                NPJ += delta[k][l];                        
                QFC += pow(alpha[k][l],2);
			}  
            if (HNca[k]==0){
				alphaG = calphaParams[move++]; 
				betaG = calphaParams[move++]; 
		        calpha[k] = 1/gsl_ran_gamma(r,alphaG+0.5*NPJ,1/(betaG+0.5*QFC));
		        //calpha[k] = 7;
		        //Rprintf("%s %i %f %f %i %f %f \n","calpha G: ",sw,alphaG,betaG,NPJ,QFC,calpha[k]);		    
		    }else if (HNca[k]==1){
                if ((sw % batchL)==0 && WB > 0){ 
	                if (acceptCa[k] > 0.25 && tuneCa[k] < FCaUL) tuneCa[k] += MIN(0.01,1/sqrt(WB)) * tuneCa[k]; 
	                if (acceptCa[k] <= 0.25 && tuneCa[k] > FCaLL) tuneCa[k] -= MIN(0.01,1/sqrt(WB)) * tuneCa[k];
	                acceptCa[k] = 0.0;   
	                if (tuneCa[k] < FCaLL) tuneCa[k] = FCaLL;
	                if (tuneCa[k] > FCaUL) tuneCa[k] = FCaUL;
                }            
                phi2G = calphaParams[move++];
                calphaP = calpha[k] + gsl_ran_gaussian(r,sqrt(tuneCa[k]));
                while (calphaP <= 0) calphaP = calpha[k] + gsl_ran_gaussian(r,sqrt(tuneCa[k]));                
                Acp = exp(-0.5*NPJ*log(calphaP/calpha[k]) + (QFC/2)*(1/calpha[k]-1/calphaP) + 
                      (calpha[k]-calphaP)/(2*phi2G));	        	            
	            unifRV = gsl_ran_flat(r,0.0,1.0);
	            //Rprintf("%s %i %f %i %f %f %f %f %f \n","calpha HN: ",sw,Acp,NPJ,calphaP,calpha[k],QFC,phi2G,tuneCa[k]);				
                if (Acp > unifRV){
                    calpha[k] = calphaP;
	                acceptCa[k] += 1/((double)batchL);                  
                }
	        } 
	    } 
	    
        // - 9 - c_psi k
        //Rprintf("%i %s \n",sw,"c_psi k");
        move = 0;
        for (k = 0; k < p; k++){
            NPJ = 0;
            QFC = 0;
            for (l = 0; l < LC; l++){
                NPJ += ksi[k][l];
                QFC += pow(psi[k][l],2);
			}
            if (HNcp[k]==0){
				alphaG = cpsiParams[move++];
				betaG = cpsiParams[move++];
		        cpsi[k] = 1/gsl_ran_gamma(r,alphaG+0.5*NPJ,1/(betaG+0.5*QFC));
		        //Rprintf("%s %i %f %f %i %f %f \n","cpsi : ",sw,alphaG,betaG,NPJ,QFC,cpsi[k]);
		    }else if (HNcp[k]==1){
                if ((sw % batchL)==0 && WB > 0){
	                if (acceptCpsi[k] > 0.25 && tuneCpsi[k] < CpUL) tuneCpsi[k] += MIN(0.01,1/sqrt(WB)) * tuneCpsi[k];
	                if (acceptCpsi[k] <= 0.25 && tuneCpsi[k] > CpLL) tuneCpsi[k] -= MIN(0.01,1/sqrt(WB)) * tuneCpsi[k];
	                acceptCpsi[k] = 0.0;
	                if (tuneCpsi[k] < CpLL) tuneCpsi[k] = CpLL;
	                if (tuneCpsi[k] > CpUL) tuneCpsi[k] = CpUL;
                }
                phi2G = cpsiParams[move++];
                calphaP = cpsi[k] + gsl_ran_gaussian(r,sqrt(tuneCpsi[k]));
                while (calphaP <= 0) calphaP = cpsi[k] + gsl_ran_gaussian(r,sqrt(tuneCpsi[k]));
                Acp = exp(-0.5*NPJ*log(calphaP/cpsi[k]) + (QFC/2)*(1/cpsi[k]-1/calphaP) +
                      (cpsi[k]-calphaP)/(2*phi2G));
	            unifRV = gsl_ran_flat(r,0.0,1.0);
	            //Rprintf("%s %i %f %i %f %f %f %f %f \n","calpha HN: ",sw,Acp,NPJ,calphaP,calpha[k],QFC,phi2G,tuneCa[k]);
                if (Acp > unifRV){
                    cpsi[k] = calphaP;
	                acceptCpsi[k] += 1/((double)batchL);
                }
	        }
	    }
	    
	    // - 10 - eta           
        //Rprintf("%i %s \n",sw,"eta");
        subMeanEta = gsl_vector_subvector(meanEta,0,Ngamma+p);
        subVarEta = gsl_matrix_submatrix(varEta,0,0,Ngamma+p,Ngamma+p);
        postMeanVarEta2(p,m,LG,tol,ceta,Ngamma,Ytilde,sigma2ij,X,gamma,RtCinv,&subMeanEta.vector,&subVarEta.matrix,U,mcm);                                   
        //print_matrix(&subVarEta.matrix);
        vecEta = gsl_vector_view_array(eta,Ngamma+p); 
        s = gsl_ran_flat(r,1.0,100000);
        sampleMN(s,Ngamma+p,&vecEta.vector,&subMeanEta.vector,&subVarEta.matrix,tol);            
        for (i = 0; i < m; i++){
            move = 0;
            for (k = 0; k < p; k++){
			    LPM[i][k] = 0;
                for (j = 0; j < (LG+1); j++)
                   if ((j==0) || (j>0 && gamma[k][j-1]==1)) LPM[i][k] += eta[move++]*X[i*(LG+1)+j];
                residuals[i][k] = Y[i*p+k]-LPM[i][k];
                if (p > 1) rtilde[i*p+k] = residuals[i][k]*sqrt(U[i][k]/sigma2ij[i][k]);
            }
        }        
        
 
        if (p > 1){   
			
			computeStStar(rtilde,time,m*p,0,p,St2);	         
            
            // - 11 - Update R 
            //Rprintf("%i %s \n",sw,"R");      
	        if ((sw % batchL)==0 && WB > 0){ 
	            //fprintf(out_file16, "%i %f %f \n", sw, acceptR, tuneR[0]);
	            if (acceptR > 0.25 && tuneR[0] > RtLL) tuneR[0] -= MIN(0.01,1/sqrt(WB)) * tuneR[0]; 
	            if (acceptR <= 0.25 && tuneR[0] < RtUL) tuneR[0] += MIN(0.01,1/sqrt(WB)) * tuneR[0];
		        acceptR = 0.0;   
	        } 
	        
	        //fprintf(out_file16, "%s %i %.2f %f \n", "sw:", sw, tuneR[0], acceptR); 	        	        	        
	        
	        gsl_matrix_memcpy(PSIt,EtC);            
            gsl_matrix_scale(PSIt,tuneR[0]-p-1);	    
            gsl_matrix_add(PSIt,St2);            
            gsl_matrix_memcpy(CopyPSIt,PSIt);
            
            //fprintf(out_file16, "%s %f %f %f \n", "RtC scaled:", gsl_matrix_get(PSIt,0,0),gsl_matrix_get(PSIt,0,1),gsl_matrix_get(PSIt,1,1));           
            
            ginv(p,tol,PSIt);
            
            //fprintf(out_file16, "%s %f %f %f \n", "inv:", gsl_matrix_get(PSIt,0,0),gsl_matrix_get(PSIt,0,1),gsl_matrix_get(PSIt,1,1));
            
            s = gsl_ran_flat(r,1.0,100000);
                      
            rwish(s,p,m+tuneR[0],PSIt,EtP);            
            
            //fprintf(out_file16, "%s %f %f %f \n", "rwish out:", gsl_matrix_get(EtP,0,0),gsl_matrix_get(EtP,0,1),gsl_matrix_get(EtP,1,1));
            
            ginv(p,tol,EtP);            
            
            //fprintf(out_file16, "%s %f %f %f \n", "EtP:", gsl_matrix_get(EtP,0,0),gsl_matrix_get(EtP,0,1),gsl_matrix_get(EtP,1,1));
            //fprintf(out_file16, "%s %f %f %f \n", "EtC:", gsl_matrix_get(EtC,0,0),gsl_matrix_get(EtC,0,1),gsl_matrix_get(EtC,1,1));
            
            decomposeEtoDS(p,0,EtP,DtP,RtP); 
            
            //fprintf(out_file16, "%s %f %f \n", "DtP:", gsl_matrix_get(DtP,0,0),gsl_matrix_get(DtP,1,1));
            //fprintf(out_file16, "%s %f \n", "RtP:", gsl_matrix_get(RtP,0,1));
            //fprintf(out_file16, "%s %f %f \n", "DtC:", gsl_matrix_get(DtC,0,0),gsl_matrix_get(DtC,1,1));
            //fprintf(out_file16, "%s %f \n", "RtC:", gsl_matrix_get(RtC,0,1));
            
            priorLogR = 0.0;
            move = 0;
           
            for (k = 0; k < (p-1); k++){
                for (l = (k+1); l < p; l++){
			        priorLogR += pow(FisherTr(gsl_matrix_get(RtP,k,l),FT[0])-theta[move],2) -
                                 pow(FisherTr(Rvec[move],FT[0])-theta[move],2);
                    move++;
                } 
            } 
            priorLogR *= -1/(2*tau2);    

            if (FT[0]==1){
                move = 0;
                for (k = 0; k < (p-1); k++){
                    for (l = (k+1); l < p; l++){
                        priorLogR += log(1+Rvec[move]) + log(1-Rvec[move]) - 
                                     log(1+gsl_matrix_get(RtP,k,l)) - log(1-gsl_matrix_get(RtP,k,l));
                        move++;
                    } 
                } 
		    }
		    
		    logLikP = logPropPdfDR(DtC,RtP,St2,DtC,p,0,m,0,1);
            logLikC = logPropPdfDR(DtC,RtC,St2,DtC,p,0,m,0,1);            
	 
	        gsl_matrix_memcpy(PSIt,EtC);            
            gsl_matrix_scale(PSIt,tuneR[0]-p-1);            	   
	        logPropDRP = logPropPdfDR(DtP,EtP,CopyPSIt,PSIt,p,p-1,m+tuneR[0]+p+1,tuneR[0],1);
            
	        gsl_matrix_memcpy(CopyPSIt,EtP);            
            gsl_matrix_scale(CopyPSIt,tuneR[0]-p-1);
	        gsl_matrix_add(CopyPSIt,St2);	    
	        
	        gsl_matrix_memcpy(PSIt,EtP);             
            gsl_matrix_scale(PSIt,tuneR[0]-p-1);	  	    
	        
	        logPropDRC = logPropPdfDR(DtC,EtC,CopyPSIt,PSIt,p,p-1,m+tuneR[0]+p+1,tuneR[0],1);

	        logAcp = logLikP - logLikC + priorLogR + logPropDRC - logPropDRP;		    		  
/*		                
            logPropDRP = logPropPdfDR(DtP,EtP,CopyPSIt,CopyPSIt,p,p-1,tuneR[0]+p+1,tuneR[0],1.0);                        
            
            gsl_matrix_memcpy(PSIt,RtP);
            gsl_matrix_scale(PSIt,tuneR[0]-p-1);  
            
            //fprintf(out_file16, "%s %f %f %f \n", "RtP scaled:", gsl_matrix_get(PSIt,0,0),gsl_matrix_get(PSIt,0,1),gsl_matrix_get(PSIt,1,1));              
            
            logPropDRC = logPropPdfDR(DtC,EtC,PSIt,PSIt,p,p-1,tuneR[0]+p+1,tuneR[0],1.0);    	        
	        
	        detR = logPropPdfDR(DtC,RtP,PSIt,RtC,p,0,m,m,0);    	  
	        gsl_matrix_memcpy(RtPinv,RtP);
	 
	        ginv(p,tol,RtPinv);

	        SPP = ScalcMult(p,m,LG,tol,ceta,Ngamma,Ytilde,sigma2ij,X,gamma,RtPinv,St,&Q2);
	        
	        logAcp = logPropDRC - logPropDRP + priorLogR + detR + (-SPP+SPC)/2;
*/
            Acp = exp(logAcp);
            unifRV = gsl_ran_flat(r,0.0,1.0); 
                       
if (1==0){  
fprintf(out_file14, "%i %.2f | %f %f %f | %f %f %f | %f %f \n%f %f %f %f %f %f %f\n",
sw, tuneR[0], 
gsl_matrix_get(EtC,0,0),gsl_matrix_get(EtC,0,1),gsl_matrix_get(EtC,1,1),
gsl_matrix_get(EtP,0,0),gsl_matrix_get(EtP,0,1),gsl_matrix_get(EtP,1,1),
gsl_matrix_get(RtP,0,1), gsl_matrix_get(RtC,0,1),
logPropDRC,logPropDRP,priorLogR,detR,SPP,SPC,Acp);	    
} 
                        
            if (Acp > unifRV){
                acceptR += 1/((double)batchL);  
	            //SPC = SPP; 
	            gsl_matrix_memcpy(RtC,RtP);
	            gsl_matrix_memcpy(DtC,DtP);
	            gsl_matrix_memcpy(EtC,EtP);	            
	            gsl_matrix_memcpy(RtCinv,RtP);	 
	            ginv(p,tol,RtCinv);	            
	            //gsl_matrix_memcpy(RtCinv,RtPinv);	            
	            move = 0;
                for (k = 0; k < (p-1); k++)
                    for (l = (k+1); l < p; l++)
                        Rvec[move++] = gsl_matrix_get(RtP,k,l);
            }
/*            
            //Rprintf("%i %s \n",sw,"R");      
	        if ((sw % batchL)==0 && WB > 0){ 
	            if (acceptR > 0.25 && tuneR[0] > RtLL) tuneR[0] -= MIN(0.01,1/sqrt(WB)) * tuneR[0]; 
	            if (acceptR <= 0.25 && tuneR[0] < RtUL) tuneR[0] += MIN(0.01,1/sqrt(WB)) * tuneR[0];
		        acceptR = 0.0;   
	        }
	        gsl_matrix_memcpy(PSIt,RtC);
            gsl_matrix_scale(PSIt,tuneR[0]-p-1);	    
            gsl_matrix_memcpy(CopyPSIt,PSIt);
            ginv(p,tol,PSIt);
            s = gsl_ran_flat(r,1.0,100000);
            rwish(s,p,tuneR[0],PSIt,EtP);
            ginv(p,tol,EtP);
            decomposeEtoDS(p,0,EtP,DtP,RtP); 
            priorLogR = 0.0;
            move = 0;
            for (k = 0; k < (p-1); k++){
                for (l = (k+1); l < p; l++){
			        priorLogR += pow(FisherTr(gsl_matrix_get(RtP,k,l),FT[0])-theta[move],2) -
                                 pow(FisherTr(Rvec[move],FT[0])-theta[move],2);
                    move++;
                } 
            } 
            priorLogR *= -1/(2*tau2);    
            if (FT[0]==1){
                move = 0;
                for (k = 0; k < (p-1); k++){
                    for (l = (k+1); l < p; l++){
                        priorLogR += log(1+Rvec[move]) + log(1-Rvec[move]) - 
                                     log(1+gsl_matrix_get(RtP,k,l)) - log(1-gsl_matrix_get(RtP,k,l));
                        move++;
                    } 
                } 
		    }
		    
if (1==0)		    fprintf(out_file16,"%s %i %i | %f %f %f %f %f %f %f \n",
"priorLogR: ",sw, FT[0], 
theta[0], gsl_matrix_get(RtC,0,1), Rvec[0], 
FisherTr(gsl_matrix_get(RtC,0,1),FT[0]),
gsl_matrix_get(RtP,0,1), 
FisherTr(gsl_matrix_get(RtP,0,1),FT[0]),
log(1+Rvec[0]) + log(1-Rvec[0]) - log(1+gsl_matrix_get(RtP,0,1)) - log(1-gsl_matrix_get(RtP,0,1)));
		    
            logPropDRP = logPropPdfDR(DtP,EtP,CopyPSIt,CopyPSIt,p,p-1,tuneR[0]+p+1,tuneR[0],1.0);                        
            
            gsl_matrix_memcpy(PSIt,RtP);
            gsl_matrix_scale(PSIt,tuneR[0]-p-1);            
            
            logPropDRC = logPropPdfDR(DtC,EtC,PSIt,PSIt,p,p-1,tuneR[0]+p+1,tuneR[0],1.0);    	        
	        
	        detR = logPropPdfDR(DtC,RtP,PSIt,RtC,p,0,m,m,0);    	  
	        gsl_matrix_memcpy(RtPinv,RtP);
	        ginv(p,tol,RtPinv);
	        
	        SPP = ScalcMult(p,m,LG,tol,ceta,Ngamma,Ytilde,sigma2ij,X,gamma,RtPinv,St,&Q2);
	        logAcp = logPropDRC - logPropDRP + priorLogR + detR + (-SPP+SPC)/2;
            Acp = exp(logAcp);
            unifRV = gsl_ran_flat(r,0.0,1.0);            
            if (1==1){  
                //print_matrix(RtC);
                fprintf(out_file16,"%s %i | %.2f %.2f %.2f | %.2f %.2f %.2f | %.3f %.3f | %.2f %.2f %.2f %.2f %.2f %.2f | %.1f\n",
"R: ",sw,
gsl_matrix_get(RtC,0,1),gsl_matrix_get(RtC,0,2),gsl_matrix_get(RtC,1,2),
gsl_matrix_get(RtP,0,1),gsl_matrix_get(RtP,0,2),gsl_matrix_get(RtP,1,2),
Acp,logAcp,logPropDRP,logPropDRC,priorLogR,detR,SPP,SPC,tuneR[0]);	    
            }                      
            if (Acp > unifRV){
                acceptR += 1/((double)batchL);  
	            SPC = SPP; 
	            gsl_matrix_memcpy(RtC,RtP);
	            gsl_matrix_memcpy(DtC,DtP);
	            gsl_matrix_memcpy(EtC,EtP);
	            gsl_matrix_memcpy(RtCinv,RtPinv);	            
	            move = 0;
                for (k = 0; k < (p-1); k++)
                    for (l = (k+1); l < p; l++)
                        Rvec[move++] = gsl_matrix_get(RtP,k,l);
            }
*/        
            // - 12 - Update Theta
            //Rprintf("%i %s \n",sw,"Theta");                         
            
            for (h1 = 0; h1 < H; h1++){			 
			    if (nmembers[h1]>0){                          
                    vecTheta = gsl_vector_view_array(vizTheta,nmembers[h1]);
                    move = 0;
                    for (k = 0; k < d; k++) 
                        if (compAlloc[k]==h1) grt[move++] = (FisherTr(Rvec[k],FT[0])/tau2 + muR[h1]/sigma2R) / (1/tau2 + 1/sigma2R);            
                    Ar = gsl_vector_view_array(grt,nmembers[h1]);                    
                    gsl_matrix_set_identity(A);                    
                    gsl_matrix_scale(A,1/(1/tau2 + 1/sigma2R));                            
                    subD = gsl_matrix_submatrix(A,0,0,nmembers[h1],nmembers[h1]);
                    s = gsl_ran_flat(r,1.0,100000); 
                    sampleMN(s,nmembers[h1],&vecTheta.vector,&Ar.vector,&subD.matrix,tol);        
                    move = 0;
                    for (k = 0; k < d; k++)                         
                        if (compAlloc[k]==h1) theta[k] = vizTheta[move++]; 
				}
			}                          
            //fprintf(out_file16,"%s %i %f %f %f \n","theta: ",sw,theta[0],sumtheta,grt[0]);				
                                    
            // - 13 - Update muR
            //Rprintf("%i %s \n",sw,"muR");            
            sstheta = 0.0;
            for (h1 = 0; h1 < H; h1++){
                sumtheta = 0.0;        
                for (k = 0; k < d; k++)
                    if (compAlloc[k]==h1) sumtheta += theta[k];
                temp = 1/(nmembers[h1]/sigma2R+1/varmuR);
                muR[h1] = temp*(sumtheta/sigma2R) + gsl_ran_gaussian(r,sqrt(temp));                             
                for (k = 0; k < d; k++)
                    if (compAlloc[k]==h1) sstheta += pow(theta[k]-muR[h1],2);
                //fprintf(out_file16,"%s %i %i %f %f %i \n","muR: ",sw,h1,temp,sstheta,nmembers[h1]);				
			}
            
            // - 14 - sigma2_R
            //Rprintf("%i %s \n",sw,"sigma2_R");
            if ((sw % batchL)==0 && WB > 0){ 
                if (acceptSigma2R > 0.25 && tuneSigma2R[0] < F1UL) tuneSigma2R[0] += MIN(0.01,1/sqrt(WB)) * tuneSigma2R[0]; 
	            if (acceptSigma2R <= 0.25 && tuneSigma2R[0] > F1LL) tuneSigma2R[0] -= MIN(0.01,1/sqrt(WB)) * tuneSigma2R[0];
	            acceptSigma2R = 0.0;   
	            if (tuneSigma2R[0] < F1LL) tuneSigma2R[0] = F1LL;
	            if (tuneSigma2R[0] > F1UL) tuneSigma2R[0] = F1UL;
            }            
            sigma2RP = sigma2R + gsl_ran_gaussian(r,sqrt(tuneSigma2R[0]));
            while (sigma2RP <= 0) sigma2RP = sigma2R + gsl_ran_gaussian(r,sqrt(tuneSigma2R[0]));            
            Acp = exp(-0.5*d*log(sigma2RP/sigma2R) + 0.5*sstheta*(1/sigma2R-1/sigma2RP) + 
                  0.5*(sigma2R-sigma2RP)/varsigma2R);
            //fprintf(out_file16,"%s %i %f %f %f %f %f \n","sigma2R: ",sw,Acp,sigma2RP,sigma2R,sstheta,varsigma2R);				      
	        unifRV = gsl_ran_flat(r,0.0,1.0);
            if (Acp > unifRV){
                sigma2R = sigma2RP;
	            acceptSigma2R += 1/((double)batchL);
            }
            
            // - 15 - Stick breaking Weights 
		    nmb = 0;
            for (h1 = 0; h1 < H-1; h1++){
                nmb += nmembers[h1];
                Vdp[h1] = gsl_ran_beta(r, (double) nmembers[h1]+1.0, (double) d-nmb+alphaDP);
                if (h1 == 0)
                    w[h1] = Vdp[h1];
                else 
                    w[h1] = Vdp[h1] * w[h1-1] * (1-Vdp[h1-1]) / Vdp[h1-1];
            }
            w[H-1] = w[H-2] * (1-Vdp[H-2]) / Vdp[H-2];	 
            
	        // - 16 - Allocation to  clusters          
            //puts("alloc");
            for (h1 = 0; h1 < H; h1++){
                for (i = 0; i < d; i++){
                    pdh[h1][i] = w[h1] * gsl_ran_gaussian_pdf(theta[i]-muR[h1],sqrt(sigma2R));                   
                    //if ((sw - burn) >= 0 && 1==0) fprintf(out_file16,"%i %i %i %f \n",sw,h1,i,pdh[h1][i]);                                    
                } 
            }

            for (i = 0; i < d; i++){
			    temp = 0.0;
                for (h1 = 0; h1 < H; h1++)
				    temp += pdh[h1][i];
				for (h1 = 0; h1 < H; h1++)
				    pdh[h1][i] /= temp;    
		    }		                 
	    
            s = gsl_ran_flat(r,1.0,100000);
            allocation(s,d,H,pdh,compAlloc,0);
            findNmembers(d,H,compAlloc,nmembers);
            //Rprintf("%i %i \n",nmembers[0],nmembers[1]);
        
            // - 17 - DP concentration parameter
            s = gsl_ran_flat(r,1.0,100000);
            newAlphaDP = updateAlpha(s,d,H,Alphaa,Alphab,TruncAlpha,nmembers,alphaDP);
            alphaDP = newAlphaDP;
	    }//if (p > 1)	    
        if (((sw - burn) >= 0) && (((sw - burn ) % thin) == 0)){
		    
            // - 18 - deviance = -2*LogLikelihood
            
            //SPC = ScalcMult(p,m,LG,tol,ceta,Ngamma,Ytilde,sigma2ij,X,gamma,RtCinv,St,&Q2,U,mcm);            
            //detR = 0.0;
            //for (i = 0; i < m; i++)
            //    for (j = 0; j < p; j++)
            //        detR += LPV[i][j] - U[i][j];
            //Rprintf("%s %f\n","sum lpv",detR);
            //for (j = 0; j < p; j++)
            //    detR += m*log(sigma2zk[j]);
            //Rprintf("%s %f\n","+m*log(s2)",detR);
            //detS = det(p,RtC);
            //detR += m*log(detS);    
            //Rprintf("%s %f\n","+m|R|",detR);
            //dev0 = SPC + (Ngamma+p)*log(ceta+1) + detR + m*p*log(2*M_PI); 
            //dev[0] += dev0;
             
	        //SPP = NormalQuadr(p,m,LG,Ngamma,Ytilde,sigma2ij,X,gamma,RtCinv,eta,U,mcm);
	        //dev1 = SPP + detR + m*p*log(2*M_PI);
            //dev[1] +=  dev1; 
            
	        //Rprintf("%s %f %f %f %f %f \n","deviance:",dev[0],m*p*log(2*M_PI),detR,
	        //(Ngamma+p)*log(ceta+1),SPC);
	        
	        if (compDev==1){                                   
                detR = m*log(det(p,RtC));            
                dev0 = mplog + detR;             
                for (j = 0; j < p; j++)
                    dev0 += m*log(sigma2zk[j]);
                for (i = 0; i < m; i++)
                    for (j = 0; j < p; j++)
                        dev0 += LPV[i][j];            
                for (i = 0; i < m; i++)
                    for (j = 0; j < p; j++)
                        dev0 -= 2*((omegaij[i][j]/2)*log(omegaij[i][j]/2) - gsl_sf_lngamma(omegaij[i][j]/2));
                for (i = 0; i < m; i++){  
				    for (j = 0; j < p; j++)
				        params[j] = Y[i*p+j] - LPM[i][j];				
				    for (j = 0; j < p; j++)
				        params[p+j] = RtCinv->data[j * RtCinv->tda + j] / sigma2ij[i][j];
				    l = 0;
				    for (k = 0; k < (p-1); k++)
                        for (j = k+1; j < p; j++)
                            params[2*p+(l++)] = RtCinv->data[k * RtCinv->tda + j] / 
                                                sqrt(sigma2ij[i][k]*sigma2ij[i][j]);
				    for (j = 0; j < p; j++)
				        params[p*(3+p)/2+j] = omegaij[i][j]/2; 
				    //if (p > 3)  p_hcubature(1,MultiNormalGammaPDF,params,p,lower,upper,0,0.0001,0,ERROR_INDIVIDUAL,&val,&err);
                    //if (p <= 3) p_pcubature(1,MultiNormalGammaPDF,params,p,lower,upper,0,0.0001,0,ERROR_INDIVIDUAL,&val,&err);     
				    if (p > 1) p_hcubature(1,MultiNormalGammaPDF,params,p,lower,upper,0,0.0001,0,ERROR_INDIVIDUAL,&val,&err);
                    if (p == 1) val = gsl_sf_gamma(params[2]+0.5)/pow(params[2]+params[0]*params[0]*params[1]/2,params[2]+0.5);
                    dev0 -= 2*log(val);		    	    
		            dev[0] += dev0;
			    }
			}	        	     
	        
	        // Write to files
            for (k = 0; k < p; k++)
                for (j = 0; j < LG; j++)
                    fprintf(out_file1, "%i ", gamma[k][j]);
            fprintf(out_file1, "\n");
            
            fprintf(out_file2, "%f \n", ceta);
            
            for (k = 0; k < p; k++)
                for (j = 0; j < LD; j++)
                    fprintf(out_file3, "%i ", delta[k][j]);
            fprintf(out_file3, "\n");
                        
            for (k = 0; k < p; k++)
                for (j = 0; j < LD; j++)
                    fprintf(out_file4, "%f ", alpha[k][j]);
            fprintf(out_file4, "\n");
                                    
            if (p > 1){
                for (k = 0; k < p; k++)
                    for (l = 0; l < p; l++)
                        fprintf(out_file5, "%f ", gsl_matrix_get(RtC,k,l)); 
                fprintf(out_file5, "\n");
                
                if (sw==LASTsw[0]){                
                    for (k = 0; k < p; k++)
                        for (l = 0; l < p; l++)
                            fprintf(out_file6, "%f ", gsl_matrix_get(DtC,k,l)); 
                    fprintf(out_file6, "\n");
                                    
                    for (k = 0; k < p; k++)
                        for (l = 0; l < p; l++)
                            fprintf(out_file6, "%f ", gsl_matrix_get(EtC,k,l)); 
                    fprintf(out_file6, "\n");			    
			    }

                //for (j = 0; j < d; j++) 
                //    fprintf(out_file6, "%f ", theta[j]); 
                //fprintf(out_file6, "\n");
             
                for (j = 0; j < H; j++) 
                    fprintf(out_file7, "%0.10f ", muR[j]);  
                fprintf(out_file7, "\n");
                                            
                fprintf(out_file8, "%0.10f \n", sigma2R);         
		    }
		    
            for (j = 0; j < p; j++)
                fprintf(out_file9, "%f ", calpha[j]);
            fprintf(out_file9, "\n");
            
            for (j = 0; j < p; j++)
                fprintf(out_file10, "%f ", sigma2zk[j]);
            fprintf(out_file10, "\n");  
            
            move = 0;
            for (k = 0; k < p; k++){                
                for (j = 0; j < (LG+1); j++)
                    if ((j==0) || (j>0 && gamma[k][j-1]==1)) fprintf(out_file11, "%f ", eta[move++]); 
                    else fprintf(out_file11, "%f ", 0.0);                
            }
            fprintf (out_file11, "\n");
            
            for (i = 0; i < d; i++)
                fprintf(out_file12, "%i ", compAlloc[i]);
            fprintf (out_file12, "\n");
            
            for (h1 = 0; h1 < H; h1++)
                fprintf(out_file13, "%i ", nmembers[h1]);
            fprintf(out_file13, "\n");
            
            //for (i = 0; i < d; i++){
            //    for (h1 = 0; h1 < H; h1++)
            //        fprintf(out_file14,"%f ",pdh[h1][i]);                                    
            //    fprintf(out_file14, "\n");
			//}
            fprintf(out_file14, "%f \n", dev0);            
			
			fprintf(out_file15, "%f \n", alphaDP); 
			
			for (k = 0; k < p; k++)
                for (j = 0; j < LC; j++)
                    fprintf(out_file16, "%f ", psi[k][j]);
            fprintf(out_file16, "\n"); 
            
            for (k = 0; k < p; k++)
                for (j = 0; j < LC; j++)
                    fprintf(out_file17, "%i ", ksi[k][j]);
            fprintf(out_file17, "\n");
            
            for (j = 0; j < p; j++)
                fprintf(out_file18, "%f ", fi2[j]);
            fprintf(out_file18, "\n");
            
            for (j = 0; j < p; j++)
                fprintf(out_file19, "%f ", cpsi[j]);
            fprintf(out_file19, "\n");                        
        }
        // If sw needs to be printed
        if ((sw==(sweeps-1)) && (!((sweeps % 1000)==0))) Rprintf("%i %s \n",sweeps, "posterior samples...");
    }//end of sw
    
    //Update LASTWB
    LASTWB[0] = WB;
    
    //Free up random number generator
    gsl_rng_free (r);

    //Close files
    fclose(out_file1); fclose(out_file2); fclose(out_file3);
    fclose(out_file4); fclose(out_file5); fclose(out_file6);
    fclose(out_file7); fclose(out_file8); fclose(out_file9);
    fclose(out_file10); fclose(out_file11); fclose(out_file12);
    fclose(out_file13); fclose(out_file14); fclose(out_file15);
    fclose(out_file16); fclose(out_file17); fclose(out_file18); 
    fclose(out_file19);

    //Free up gsl matrices
    gsl_matrix_free(PSIt); gsl_matrix_free(CopyPSIt);
    gsl_matrix_free(EtC); gsl_matrix_free(DtC); gsl_matrix_free(RtC); 
    gsl_matrix_free(EtP); gsl_matrix_free(DtP); gsl_matrix_free(RtP);
    gsl_matrix_free(RtCinv); gsl_matrix_free(RtPinv);
    gsl_matrix_free(A); gsl_matrix_free(D); gsl_matrix_free(varEta); 
    gsl_vector_free(meanEta);  gsl_vector_free(alphaHat);
    gsl_vector_free(alphaP); gsl_matrix_free(St); gsl_matrix_free(StP);
    gsl_matrix_free(St2);

    //Free jagged vectors
	for (j = 0; j < NG; j++) {free(indexG[j]); free(vecGamma[j]); free(vecGammaP[j]);}
	for (j = 0; j < ND; j++) {free(indexD[j]); free(vecDelta[j]); free(vecDeltaP[j]);}
    for (j = 0; j < NC; j++) {free(indexX[j]); free(vecKsi[j]); free(vecKsiP[j]);}
}
