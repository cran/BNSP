/* LongMultG.c Bayesian Nonparametric Long Multivariate Regression
 * Copyright (C) 2018  Georgios Papageorgiou, gpapageo@gmail.com 
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
/*
#include "matalg.h"
#include "pdfs.h"
#include "sampling.h"
#include "other.functions.h"
#include "mathm.h"
#include "spec.BCM.h"  
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

extern double ScalcMultLong(int m, int p, double tol, int LG, int Ngamma, int niMax, int *niVec, int *cusumniVec, int N, double sigma2ij[N][p], int dimLPC, double LPC[m][dimLPC][p][p], double *Y, double *X, int gamma[p][LG], gsl_matrix *RiAll, int *intime, double ceta, double *qf2);
extern double SPh(int T, int d, int H, int h, double tol, double *thetaTilde, int LG, int gamma[H][LG], int Ngamma, int *compAlloc, int nmembers, double ceta, double *AllBases, double *LPV, double *qf2);
extern void postMeanVarEtaLong(int m, int p, double tol, int LG, int Ngamma, int niMax, int *niVec, int *cusumniVec, int N, double sigma2ij[N][p], int dimLPC, double LPC[m][dimLPC][p][p], double *Y, double *X, int gamma[p][LG], gsl_matrix *RiAll, int *intime, double ceta, gsl_vector *MeanEta, gsl_matrix *varEta);
extern void cRes(int p, int m, int LG, int gamma[p][LG], int Ngamma, double *X, gsl_vector *MeanEta, double *Y, double *sqRes, double *BaseXg);
extern void cResCheck(int p, int m, int N, int niMax, int *niVec, int *cusumniVec, double sigma2ij[N][p], int dimLPC, double LPC[m][dimLPC][p][p], double *ResC, double *ResCheck, int SL);
extern void postMeanVarPSI(int m, int p, double tol, int *niVec,  int *cusumniVec, int LK, int NKsi, double *rsd, int ksi[p][p][LK], int *cusumC, double *C, int N, double sigma2ij[N][p], gsl_matrix *RiAll, int *intime, double *cpsi, gsl_vector *Mean, gsl_matrix *Var);
extern void MNCond(double tol, int start, int end, gsl_vector *mu, gsl_matrix *Sigma, double *W, gsl_vector *condMu, gsl_matrix *condSigma);
extern void iCovPostThetah(int T, int d, int H, int h, double tol, int LG, int gamma[H][LG], int *compAlloc, int nmembers, int *Ngamma, double ceta, double sigma2, double tau2, double *AllBases, double *LPV, gsl_matrix *A);
extern void postMeanVarEtaH(int T, int d, int H, int h, double tol, int LG, int gamma[H][LG], int *compAlloc, int *nmembers, int *Ngamma, double sigma2, double ceta, double *LPV, double *AllBases, double *thetaTilde, gsl_vector *MeanEta, gsl_matrix *varEta, int sw, int kk);
extern void setBaseZBSgh(int T, int d, int H, int h, int LG, int gamma[H][LG], int *compAlloc, double *AllBases, double *BaseZBSgh);
extern void cSqResh(int T, int d, int h, int *compAlloc, int *gamma, int *Ngamma, int LG, double *AllBases, gsl_vector *MeanEta, double *theta, double *sqRes);
extern void DeltaAlphaHat(int T, int d, double tol, double *LPV, double *sqRes, int *delta, int Ndelta, int start, int end, double *AllBases, double sigma2, double *sigma2t, double calpha, gsl_matrix *D, gsl_vector *alphaHat);

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX_PATH 300
 
void  longmultg(int *seed1, char **WorkingDir, 
                int *sbtmpN, 
                int *niVec, int *cusumniVec, int *intime2, int *intime, int *niMaxLUT, int *FUT,                
                int *cusumC, double *C, double *Y, double *X, double *Z, double *Xc, double *Zc, 
                int *LGDKc, int *NGDKc, int *vecLG, int *vecLD, int *vecLK, int *vecLGc, int *vecLDc,
                int *cusumVecLG, int *cusumVecLD, int *cusumVecLK, int *cusumVecLGc, int *cusumVecLDc,                  
                double *blockSizeProb, int *maxBSGDKc,                
                double *tuneCa, double *tuneSigma2, double *tuneCb, double *tuneAlpha, double *tuneSigma2R, 
                double *tuneR, double *tuneCpsi, double *tuneCbCor, double *tuneOmega, double *tuneComega, 
                double *pimu, double *cetaParams, double *pisigma, int *HNca, double *calphaParams,
                int *HNszk, double *szkParams, double *piphi, int *HNcpsi, double *cpsiParams, 
                double *cetaCorParams, int *HNco, double *comegaParams, double *pinufi,
                int *HNsigcor, double *sigmaCorParams,
                double *tau1, int *FT, double *dev, int *isDz,
                int *contParams, double *LASTAll,
                int *H1, double *DPparams)
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
         *out_file19, *out_file20, *out_file21, *out_file22, *out_file23, *out_file24, *out_file25;

    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.gamma.txt");
    out_file1 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.cbeta.txt");
    out_file2 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.delta.txt");
    out_file3 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.alpha.txt");
    out_file4 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.sigma2.txt");
    out_file5 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.calpha.txt");
    out_file6 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.beta.txt");
    out_file7 = fopen(path_name, "a"); 
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.psi.txt");
    out_file8 = fopen(path_name, "a");    
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.ksi.txt");
    out_file9 = fopen(path_name, "a");    
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.cpsi.txt"); 
    out_file10 = fopen(path_name, "a");    
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.R.txt");
    out_file11 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.nu.txt");
    out_file12 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.fi.txt");
    out_file13 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.omega.txt");
    out_file14 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.sigma2R.txt");
    out_file15 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.ceta.txt");
    out_file16 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.comega.txt");
    out_file17 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.eta.txt");
    out_file18 = fopen(path_name, "a");    
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.deviance.txt");
    out_file19 = fopen(path_name, "a");           
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.DE.txt");
    out_file20 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.compAlloc.txt");
    out_file21 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.nmembers.txt");
    out_file22 = fopen(path_name, "a");    
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.DPconc.txt");
    out_file23 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.probs.txt");
    out_file24 = fopen(path_name, "a");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.tune.txt");
    out_file25 = fopen(path_name, "a");

    // Sweeps, burn-in period, thin
    int sweeps = sbtmpN[0]; //number of posterior samples
    int burn = sbtmpN[1]; // integer burn-in period
    int thin = sbtmpN[2]; // integer thin

    // Dimensions
    int m = sbtmpN[3];
    int p = sbtmpN[4];
    int N = sbtmpN[5];
    int niMax = niMaxLUT[0];
    int LUT = niMaxLUT[1];    
    int LG = LGDKc[0];
    int LD = LGDKc[1];
    int LK = LGDKc[2];
    int LGc = LGDKc[3];
    int LDc = LGDKc[4];    
    int NG = NGDKc[0];
    int ND = NGDKc[1];
    int NK = NGDKc[2];
    int NGc = NGDKc[3];
    int NDc = NGDKc[4];
    int MVLD = sbtmpN[6];
    int d = p*(p-1)/2;
    int H = H1[0];//number of surfaces
    int mcm = sbtmpN[7];

    //Tolerance level
    double tol = 0.000000015; //tolerance level for eigen values when inverting matrices

    //Declare variables for loops:
    int sw, h, i, j, k, l, t, j2, k2, l2, move, move2;

    // Prior parameters
    //c.beta
    double alphaeta = cetaParams[0];
    double betaeta = cetaParams[1];
    //Rprintf("%s %f %f \n","c.beta: ",alphaeta,betaeta);

    //pi_mu & pi_sigma
    double cmu[p][NG];
    double dmu[p][NG];
    double csigma[p][ND];
    double dsigma[p][ND];
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
	//Rprintf("%s %f %f %f %f %f %f %f %f \n","c-d mu: ",cmu[0][0],dmu[0][0],cmu[0][1],dmu[0][1],cmu[1][0],dmu[1][0],cmu[1][1],dmu[1][1]);
	//Rprintf("%s %f %f %f %f %f %f %f %f \n","c-d sg",csigma[0][0],dsigma[0][0],csigma[0][1],dsigma[0][1],csigma[1][0],dsigma[1][0],csigma[1][1],dsigma[1][1]);		
	
	//pi_nu & pi_fi
    double cnu[NGc];
    double dnu[NGc];
    double cfi[NDc];
    double dfi[NDc];     
    move = 0;
    for (k = 0; k < NGc; k++){
        cnu[k] = pinufi[move++];
        dnu[k] = pinufi[move++];
	}    
    for (k = 0; k < NDc; k++){
        cfi[k] = pinufi[move++];
        dfi[k] = pinufi[move++];
	}

    //c.alpha & sigma^2_{zero k} & c.psi & sigma^2_cor
    double phi2G;  //for half normal prior
    double alphaG, betaG; //for IG prior
    
    //pi_phi
    double cphi[p][p][NK];
    double dphi[p][p][NK];
    move = 0;
	for (k = 0; k < p; k++){
		for (l = 0; l < p; l++){
	        for (t = 0; t < NK; t++){
                cphi[k][l][t] = piphi[move++];
                dphi[k][l][t] = piphi[move++];
                //Rprintf("%s %f %f \n","c-d phi: ",cphi[k][l][t],dphi[k][l][t]);
	        }
		}
    }
    
    // DP alpha ~ Gamma(a,b) I[alpha>c]
    double Alphaa = DPparams[0];  
    double Alphab = DPparams[1];
    double TruncAlpha = DPparams[2];

    //tau2
    double tau2 = tau1[0]*tau1[0];

    // Declare quantities of interest
    double eta[p*(LG+1)];
    int gamma[p][LG];
    int Ngamma;
    int *vecGamma[NG];
    for (j = 0; j < NG; j++)
        vecGamma[j] = malloc(sizeof(int) * vecLG[j]);  
    int delta[p][LD];
    int *vecDelta[ND];
    for (j = 0; j < ND; j++) 
        vecDelta[j] = malloc(sizeof(int) * vecLD[j]);
    double alpha[p][LD];
    double calpha[p];
    double ceta;
    double sigma2zk[p];
    int ksi[p][p][LK];  
    int Nksi;
    int *vecKsi[NK];  
    for (j = 0; j < NK; j++)
        vecKsi[j] = malloc(sizeof(int) * vecLK[j]);
    double psi[p][p][LK];
    double cPsi[p*p]; 
    double Rt[LUT][p*p];     
    double etaCor[LGc+1];
    int gammaCor[H][LGc];
    int NgammaCor[H];
    int NgammaT;   
    int *vecGammaCor[NGc];
    for (j = 0; j < NGc; j++)
        vecGammaCor[j] = malloc(sizeof(int) * vecLGc[j]);
    int deltaCor[LDc];
    int NdeltaCor;
    int *vecDeltaCor[NDc];
    for (j = 0; j < NDc; j++)
        vecDeltaCor[j] = malloc(sizeof(int) * vecLDc[j]);
    int vecNdeltaCor[NDc];
    double sigma2t[LUT];
    double omega[LDc];
    double sigma2cor;
    double comega;
    double cetaCor;
    double alphaDP; //concentration parameter
    int compAlloc[d];
    double pdh[H][d];
    int nmembers[H];      

    // Declare other quantities
    double LPV[N][p]; 
    double LPVP[N][p];
    double sigma2ij[N][p];
    double sigma2ijP[N][p];
    double alphaPD[LD];
    double BaseSubAlpha[MVLD];
    double sqResC[N*p];
    double sqResP[N*p];
    double ResC[N*p];
    double ResP[N*p];
    double ResCheck[N*p];
    double priorLogR, Acp, unifRV, QFC, QFP, QFD, detR, //detS, 
           logMVNormC, logMVNormP, SPC, SPP, temp, sigma2P, dev0, dev1,
           logLikP, logLikC, logPropDRP,logPropDRC, logAcp;
    int NPJ, NCJ, start;
    double cetahat, Sprime, SDprime, elPrime, elDPrime, Q2, cetaP, calphaP, cPsiP, omegaP;
    double W[p*p*LK];
    int dimLPC = niMax*(niMax-1)/2;
    double LPC[m][dimLPC][p][p];
    double LPCP[m][dimLPC][p][p];     
    int ksiP[p][p][LK]; 
    double psiPK[LK];
    double eta2[p*(LG+1)];
    double baseDt[p*p];
    double baseEt[p*p];
    double baseRt[p*p];
    double Dt[LUT][p*p];
    double Et[LUT][p*p];
    double rt[LUT*d];
    double theta[LUT*d];    
    double thetaTilde[LUT*d];
    double thetaTildeP[LUT*d];    
    double sigma2tP[LUT];
    double LPVcor[LUT];
    double LPVPcor[LUT];
    double omegaPD[LDc];
    double sqResCCor[LUT*d];
    double sqResPCor[LUT*d];
    double w[H]; //prior probabilities
    double Vdp[H-1]; //V[i] ~ Beta(1,alpha)
    double newAlphaDP; //store concentration parameter
    int nmb; //cumulative number of members
    double vizrt[LUT*d];
    double vizTheta[LUT*d];
    double vizZstar[(LUT*d)*(LGc+1)];
    double Zeta[LUT];
    double U[N][p];
    
    //Selecting block size: gamma_B, delta_B, ksi_B
    int block, blockSize;
    double nBlocks;    
    int maxBSG = maxBSGDKc[0];
    double blockSizeProbG[maxBSG];
    for (k = 0; k < maxBSG; k++)
        blockSizeProbG[k] = blockSizeProb[k];
    unsigned int vecBSG[maxBSG];    
    int maxBSD = maxBSGDKc[1];
    double blockSizeProbD[maxBSD];
    for (k = 0; k < maxBSD; k++)
        blockSizeProbD[k] = blockSizeProb[LG + k];
    unsigned int vecBSD[maxBSD];    
    int maxBSK = maxBSGDKc[2];
    double blockSizeProbK[maxBSK];
    for (k = 0; k < maxBSK; k++)
        blockSizeProbK[k] = blockSizeProb[LG + LD + k];
    unsigned int vecBSK[maxBSK];    
    int maxBSGc = maxBSGDKc[3];
    double blockSizeProbGc[maxBSGc];
    for (k = 0; k < maxBSGc; k++)
        blockSizeProbGc[k] = blockSizeProb[LG + LD + LK + k];
    unsigned int vecBSGc[maxBSGc];    
    int maxBSDc = maxBSGDKc[4];
    double blockSizeProbDc[maxBSDc];
    for (k = 0; k < maxBSDc; k++)
        blockSizeProbDc[k] = blockSizeProb[LG + LD + LK + LGc + k];
    unsigned int vecBSDc[maxBSDc];
    
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
    
    //Ksi indeces: for shuffling
    int *indexK[NK];
    for (j = 0; j < NK; j++)
        indexK[j] = malloc(sizeof(int) * vecLK[j]);
    for (j = 0; j < NK; j++)
        for (k = 0; k < vecLK[j]; k++)
            indexK[j][k] = k;
     
    //Gamma Cor indeces: for shuffling
    int *indexGc[NGc];
    for (j = 0; j < NGc; j++)
        indexGc[j] = malloc(sizeof(int) * vecLGc[j]);
    for (j = 0; j < NGc; j++)
        for (k = 0; k < vecLGc[j]; k++)
            indexGc[j][k] = k;

    //Delta Cor indeces: for shuffling
    int *indexDc[NDc];
    for (j = 0; j < NDc; j++)
        indexDc[j] = malloc(sizeof(int) * vecLDc[j]);
    for (j = 0; j < NDc; j++)
        for (k = 0; k < vecLDc[j]; k++)
            indexDc[j][k] = k;
    
    //Proposed gamma, delta, ksi
    int NgammaP;
    int *vecGammaP[NG];
    for (j = 0; j < NG; j++)
        vecGammaP[j] = malloc(sizeof(int) * vecLG[j]);
    int *vecDeltaP[ND];
    for (j = 0; j < ND; j++)
        vecDeltaP[j] = malloc(sizeof(int) * vecLD[j]);
    int NksiP;
    int *vecKsiP[NK];
    for (j = 0; j < NK; j++)
        vecKsiP[j] = malloc(sizeof(int) * vecLK[j]);        
    //int gammaPCor[LGc];
    int *vecGammaPCor[NGc];
    for (j = 0; j < NGc; j++)
        vecGammaPCor[j] = malloc(sizeof(int) * vecLGc[j]);
    int NgammaPCor;        
    int deltaPCor[LDc];
    int NdeltaPCor;
    int *vecDeltaPCor[NDc];
    for (j = 0; j < NDc; j++)
        vecDeltaPCor[j] = malloc(sizeof(int) * vecLDc[j]);
    int gammaMP[H][LGc]; //gamma matrix proposed
    int gammaC[LGc];
    int compAlloc2[d];
    for (j = 0; j < d; j++)
        compAlloc2[j] = H+2;    
    //
    double *BaseXg = (double*) malloc (N*p*p*(LG+1) * sizeof(double));
        
    // Declare and allocate gsl vectors and matrices
    gsl_matrix *St = gsl_matrix_alloc(p,p);
    gsl_matrix *PSIt = gsl_matrix_alloc(p,p);
    gsl_matrix *CopyPSIt = gsl_matrix_alloc(p,p);
    gsl_matrix *gEt = gsl_matrix_alloc(p,p);
    gsl_matrix *gDt = gsl_matrix_calloc(p,p);
    gsl_matrix *gRt = gsl_matrix_alloc(p,p);
    
    gsl_matrix *D = gsl_matrix_alloc(MVLD,MVLD);
    gsl_matrix *varEta = gsl_matrix_alloc(p*(LG+1),p*(LG+1));
    gsl_matrix *L = gsl_matrix_alloc(p*niMax,p*niMax);
    gsl_matrix *RiAll = gsl_matrix_alloc(p,p*LUT);
    gsl_matrix *varPsi = gsl_matrix_alloc(p*p*LK,p*p*LK);
    gsl_vector *meanEta = gsl_vector_alloc(p*(LG+1));
    gsl_vector *alphaHat = gsl_vector_alloc(MVLD);
    gsl_vector *alphaP = gsl_vector_alloc(MVLD);
    gsl_vector *meanPsi = gsl_vector_alloc(p*p*LK);
    
    gsl_matrix *A = gsl_matrix_alloc(LUT*d,LUT*d);
    gsl_matrix *varEta2 = gsl_matrix_alloc(LGc+1,LGc+1);
    gsl_vector *mutheta = gsl_vector_alloc(LUT*d);
    gsl_vector *meanEta2 = gsl_vector_alloc(LGc+1);
    
    gsl_vector *vecZero1 = gsl_vector_calloc(LUT);
    
    // GSL matrix and vector views
    gsl_matrix_view subD, subVarEta, subVarPsi, EtC, DtC, RtC, subA, Zstar;
    gsl_vector_view subAlphaHat, subAlphaP, subMeanEta, subMeanPsi, vecEta, vecTheta, Ar, subMutheta, vecMean;       
    
    //Make adaptive
    int batchL = 50; //batch length
    double WB; //which batch

    double acceptCeta = 0.0;
    double GLL = 0.01;
    double GUL = 1000;

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
    double SZKUL = 300;

    double acceptCa[p];
    for (k = 0; k < p; k++)
	    acceptCa[k] = 0.0;
    double FCaLL = 0.01;
    double FCaUL = 200;
    
    double acceptCpsi[p*p];
    for (k = 0; k < (p*p); k++)
	    acceptCpsi[k] = 0.0;
    double FCpsiLL = 0.01;
    double FCpsiUL = 200;
    
    double acceptRt[LUT]; 
    for (t = 0; t < LUT; t++)
        acceptRt[t] = 0.0;    
    double RtLL = p+2;
    double RtUL = 999999999999999;
    
    double acceptCetaCor = 0.0;
    double GLLcor = 0.001;
    double GULcor = 200;
    
    double acceptSigma2 = 0.0;
    double s2LL = 0.001;//0.00000000001;
    double s2UL = 200;
    
    double acceptOmega[NDc];
	for (j = 0; j < NDc; j++)
	    acceptOmega[j] = 0.0;	
    double omegaLL = 0.5;//0.00000000001;
    double omegaUL = 200;
    
    double acceptComega = 0.0;
    double coLL = 0.01;
    double coUL = 200;
    
    // Sampler initialization
    move2 = 0;
    
    // - 0 - R related: RiAll, Et, Dt, Rt, rt, theta
    //if (contParams[0]==1){
		for (t = 0; t < LUT; t++){
            for (j = 0; j < p*p; j++){            
                Rt[t][j] = LASTAll[move2++];            
                baseRt[j] = Rt[t][j];
			}                       
            subD = gsl_matrix_submatrix(RiAll,0,t*p,p,p);
            RtC = gsl_matrix_view_array(baseRt,p,p);
            gsl_matrix_memcpy(&subD.matrix,&RtC.matrix);            
            ginv(p,tol,&subD.matrix);
		}		
        for (t = 0; t < LUT; t++)
            for (j = 0; j < p*p; j++)            
                Dt[t][j] = LASTAll[move2++];                
        for (t = 0; t < LUT; t++)
            for (j = 0; j < p*p; j++)            
                Et[t][j] = LASTAll[move2++];        		
	//}else{        	                                
    //    for (t = 0; t < LUT; t++){
    //        for (j = 0; j < p*p; j++){            
    //            Et[t][j] = 0.0;
    //            Dt[t][j] = 0.0;            
    //            Rt[t][j] = 0.0;
	//	    }
	//	    for (j = 0; j < p; j++){
    //            Et[t][p*j+j] = 1;
    //            Dt[t][p*j+j] = 1;            
    //            Rt[t][p*j+j] = 1;
	//	    }
	//    }
	//    for (i = 0; i < LUT; i++){    
    //        subD = gsl_matrix_submatrix(RiAll,0,i*p,p,p);      
    //        gsl_matrix_set_identity(&subD.matrix);
	//    }
    //}
    
    for (t = 0; t < LUT; t++){
        move = 0;
        for (k = 0; k < (p-1); k++){
            for (l = (k+1); l < p; l++){
                rt[move+t*d] = Rt[t][k*p+l];
                move++;
            }
        }
    }
    
    for (t = 0; t < LUT; t++)
        for (i = 0; i < d; i++)
            theta[i+t*d] = FisherTr(rt[i+t*d],FT[0]) + gsl_ran_gaussian(r,tau1[0]);

    // - 1 - Gamma
    if (contParams[0]==1){
        for (j = 0; j < p; j++)
            for (k = 0; k < LG; k++)
                gamma[j][k] = LASTAll[move2++];
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
    if (contParams[0]==1){
        for (j = 0; j < p; j++)
            for (k = 0; k < LD; k++)
                delta[j][k] = LASTAll[move2++];
        for (j = 0; j < p; j++)
            for (k = 0; k < LD; k++)
                alpha[j][k] = LASTAll[move2++];
    }else{    
        for (j = 0; j < p; j++){
            for (k = 0; k < LD; k++){
                delta[j][k] = 0;
                alpha[j][k] = 0;
		    }
	    }
	}

	// - 3 - Sigma2zk
	if (contParams[0]==1){
        for (j = 0; j < p; j++)
            sigma2zk[j] = LASTAll[move2++];
    }else{        
        for (j = 0; j < p; j++)
            sigma2zk[j] = 1.0;
	}

    // - 4 - LPV, sigma2ij 
    for (i = 0; i < m; i++){
		for (t = 0; t < niVec[i]; t++){
		    for (j = 0; j < p; j++){
                LPV[cusumniVec[i]+t][j] = 0;
                for (k = 0; k < LD; k++)
                    LPV[cusumniVec[i]+t][j] += alpha[j][k]*Z[N+k*N+cusumniVec[i]+t];
                sigma2ij[cusumniVec[i]+t][j] = sigma2zk[j]*exp(LPV[cusumniVec[i]+t][j]);
                sigma2ijP[cusumniVec[i]+t][j] = sigma2ij[cusumniVec[i]+t][j];             
			}
        }
    }
           
    // - 5 - Ksi & Psi
    if (contParams[0]==1){
        for (j = 0; j < p; j++){
			for (k = 0; k < p; k++){
                for (l = 0; l < LK; l++){
                    ksi[j][k][l] = LASTAll[move2++];
                    ksiP[j][k][l] = ksi[j][k][l];                                     
		        }
	        }
	    }
	    for (j = 0; j < p; j++)
			for (k = 0; k < p; k++)
                for (l = 0; l < LK; l++)               
                    psi[j][k][l] = LASTAll[move2++];                    
    }else{    
        for (j = 0; j < p; j++){
			for (k = 0; k < p; k++){
                for (l = 0; l < LK; l++){
                    ksi[j][k][l] = 0;
                    ksiP[j][k][l] = ksi[j][k][l];
                    psi[j][k][l] = 0; //(j+k+l)/100;
		        }
	        }
	    }
    }       
    
    Nksi = 0;
    for (j = 0; j < p; j++)
        for (k = 0; k < p; k++)
            for (l = 0; l < LK; l++)
                Nksi += ksi[j][k][l];
	  
    // - 6 - LPC
    move = 0;
    for (i = 0; i < m; i++){	
		for (j = 0; j < (niVec[i]*(niVec[i]-1)/2); j++){			
            for (k = 0; k < p; k++){
				for (l = 0; l < p; l++){    
                    LPC[i][j][k][l] = 0;
                    for (k2 = 0; k2 < LK; k2++)
                        LPC[i][j][k][l] += C[k2+move*LK]*psi[k][l][k2];
                    LPCP[i][j][k][l] = LPC[i][j][k][l];
				}
	        }
	        move++;  
		}
    }
    
    
    
    // - 8 - c_eta    
    if (contParams[0]==0){ 
        ceta = 1; 
        cetahat = ceta;        
        SPC = ScalcMultLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPC,Y,X,gamma,RiAll,intime,ceta,&Q2);
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
        ceta = LASTAll[move2++];
	}
   	//Rprintf("%s %i %i %i %i %i %i %f %f \n"," - 8 - ",m,p,LG,Ngamma,niMax,N,Q2,SPC);
       
	// - 9 - calpha
    move = 0;
    if (contParams[0]==0){
        for (k = 0; k < p; k++){        
            if (HNca[k]==0){
                alphaG = calphaParams[move++]; betaG = calphaParams[move++];
                if (alphaG > 1.0) {calpha[k] = betaG/(alphaG-1);}else {calpha[k] = 110.0;}
            }else{
			    phi2G = calphaParams[move++];
			    calpha[k] = sqrt(2*phi2G/M_PI);
		    }
	    }
    }else{
		for (k = 0; k < p; k++)
		    calpha[k] = LASTAll[move2++];
	}
    
    // - 10 - cPsi
    move = 0;
    if (contParams[0]==0){
        for (k = 0; k < (p*p); k++){
            if (HNcpsi[k]==0){
                alphaG = cpsiParams[move++]; betaG = cpsiParams[move++];
                if (alphaG > 1.0) {cPsi[k] = betaG/(alphaG-1);}else{cPsi[k] = 110.0;}
            }else{
	            phi2G = cpsiParams[move++];
		        cPsi[k] = sqrt(2*phi2G/M_PI);		    
	        }
		}
	}else{
		for (k = 0; k < (p*p); k++)
			cPsi[k] = LASTAll[move2++];
	}
	//Rprintf("%s %i %f %f \n","cpsi main: ",k,cpsiParams[k],cPsi[k]); 
    
/*
    gsl_matrix *cs = gsl_matrix_alloc(4,4);
    gsl_vector *cm = gsl_vector_alloc(4);
        
          
    double base123[] = {91.22, 85.05, 90.46, 90.88, 96.73, 94.80, 92.08, 75.95,  88.34,  80.41,
                        85.05,165.93,175.82,175.74,191.82,171.13,172.78,155.75, 174.20, 170.89,
                        90.46,175.82,296.65,292.42,305.21,294.54,310.51,275.81, 305.72, 303.57,
                        90.88,175.74,292.42,396.46,401.97,387.13,426.62,368.39, 401.79, 396.14,
                        96.73,191.82,305.21,401.97,510.32,480.15,518.17,475.77, 524.75, 491.26,
                        94.80,171.13,294.54,387.13,480.15,568.26,600.55,552.54, 605.19, 568.11,
                        92.08,172.78,310.51,426.62,518.17,600.55,714.22,648.63, 733.22, 709.45,
                        75.95,155.75,275.81,368.39,475.77,552.54,648.63,703.70, 811.01, 751.21,
                        88.34,174.20,305.72,401.79,524.75,605.19,733.22,811.01,1081.90,1005.87,
                        80.41,170.89,303.57,396.14,491.26,568.11,709.45,751.21,1005.87,1028.68};
            
                                                                             
    double base321[] = { 9,8,7,6,5,4,3,2,1,0 };
    vecEta = gsl_vector_view_array(base321,10);    
    subD = gsl_matrix_view_array(base123,10,10);
    
    W[0] = 0; W[1] = 1; W[2] = 2; W[3] = 7;	W[4] = 8; W[5] = 9;   
	  
	MNCond(tol, 3, 6, &vecEta.vector, &subD.matrix, W, cm, cs);
	 
	//print_vector(cm);
*/
    // - 11 - Vec of gammas Cor
    NgammaT = 0;
    for (h = 0; h < H; h++){      
        for (k = 0; k < LGc; k++){
            if (contParams[0]==0){ 
                gammaCor[h][k] = 1;//gsl_ran_bernoulli(r,cnu/(cnu+dnu));
			}else{
				gammaCor[h][k] = LASTAll[move2++];
			}
            gammaMP[h][k] = gammaCor[h][k];
		} 
        NgammaCor[h] = 0;
        for (i = 0; i < LGc; i++)
            NgammaCor[h] += gammaCor[h][i];
        NgammaT += NgammaCor[h]; 
    } 
    
    // - 12 - Cluster allocation, alphaDP, Vdp & w
    if (contParams[0]==0){
        for (i = 0; i < d; i++)
            for (h = 0; h < H; h++)
                pdh[h][i] = 1/((double)H);    
        s = gsl_ran_flat(r,1.0,100000);
        allocation(s,d,H,pdh,compAlloc,0);
	}else{
	    for (i = 0; i < d; i++)
            compAlloc[i] = LASTAll[move2++];
	}    
    findNmembers(d,H,compAlloc,nmembers);
    if (contParams[0]==0){ 
        alphaDP = 0.5;
	}else{ 
        alphaDP = LASTAll[move2++];
	}        
    nmb = 0;
    for (h = 0; h < H-1; h++){
        nmb += nmembers[h];
        Vdp[h] = gsl_ran_beta(r, (double) nmembers[h]+1.0, (double) d-nmb+alphaDP);
        if (h == 0)
            w[h] = Vdp[h];
        else 
            w[h] = Vdp[h] * w[h-1] * (1-Vdp[h-1]) / Vdp[h-1];
    }
    w[H-1] = w[H-2] * (1-Vdp[H-2]) / Vdp[H-2]; 
        
    // - 13 - Vec of deltas Cor
    if (contParams[0]==0){
        for (j = 0; j < NDc; j++)
            for (k = 0; k < vecLDc[j]; k++)
                vecDeltaCor[j][k] = 0; //gsl_ran_bernoulli(r,cfi[j]/(cfi[j]+dfi[j]));
	}else{
	    for (j = 0; j < NDc; j++)
            for (k = 0; k < vecLDc[j]; k++)
                vecDeltaCor[j][k] = LASTAll[move2++];
	}
    for (j = 0; j < NDc; j++){
	    vecNdeltaCor[j] = 0;
        for (k = 0; k < vecLDc[j]; k++)
            vecNdeltaCor[j] += vecDeltaCor[j][k];
    }
    for (j = 0; j < NDc; j++)
	    for (k = 0; k < vecLDc[j]; k++)
	        deltaCor[cusumVecLDc[j]+k] = vecDeltaCor[j][k];        
    NdeltaCor = 0;
    for (j = 0; j < NDc; j++)
        NdeltaCor += vecNdeltaCor[j];
    for (k = 0; k < LDc; k++) 
        deltaPCor[k] = deltaCor[k];
    
    // - 14 - c_omega
    if (contParams[0]==0){
        if (HNco[0]==0){ 
	        alphaG = comegaParams[0]; betaG = comegaParams[1];
		    if (alphaG > 1.0) {comega = betaG/(alphaG-1);}else{comega = 110.0;}		    
	    }else{
	        phi2G = comegaParams[0];
            comega = sqrt(2*phi2G/M_PI); 
	    }
	}else{
	    comega = LASTAll[move2++];	
	} 
	
	// - 15 - sigma2cor, omega, LPVcor, sigma2t, QFC
    if (contParams[0]==0){ 
        sigma2cor = 1;
        for (j = 0; j < LDc; j++){
            if (deltaCor[j]==0){ 
                omega[j] = 0.0;
			}else{ 
				omega[j] = gsl_ran_gaussian(r,sqrt(comega));
			}
        }
	}else{
		sigma2cor = LASTAll[move2++];
		for (j = 0; j < LDc; j++)
            omega[j] = LASTAll[move2++];
	}    
    for (t = 0; t < LUT; t++){
        LPVcor[t] = 0.0;
        for (j = 0; j < LDc; j++)
            LPVcor[t] += omega[j]*Zc[LUT*d+j*LUT*d+t*d];
        sigma2t[t] = sigma2cor*exp(LPVcor[t]);
    }
    for (t = 0; t < LUT; t++)
        for (i = 0; i < d; i++)
            thetaTilde[i+t*d] = exp(-LPVcor[t]/2)*theta[i+t*d];
            
    // - 16 - c_eta (Cor)
    if (contParams[0]==0){
        cetaCor = 1;
        cetahat = cetaCor; //starting value is the current value
        Q2 = 0;
        SPC = 0;       
        for (h = 0; h < H; h++){          
            if (nmembers[h]>0){
                SPC += SPh(LUT,d,H,h,tol,thetaTilde,LGc,gammaCor,NgammaCor[h],compAlloc,nmembers[h],cetaCor,Xc,LPVcor,&temp);
                Q2 += temp;
            }
	    }    
        elPrime = 99.9;
        while(elPrime > 0.000000001 || -elPrime > 0.000000001){
            Sprime = -Q2/pow(cetahat+1,2);
            SDprime = 2*Q2/pow(cetahat+1,3);
            elPrime = -0.5*(NgammaT+H)/(cetahat+1)-Sprime/(2*sigma2cor)-(cetaCorParams[0]+1)/cetahat+cetaCorParams[1]/pow(cetahat,2);
            elDPrime = 0.5*(NgammaT+H)/(pow(cetahat+1,2))-SDprime/(2*sigma2cor)+(cetaCorParams[0]+1)/pow(cetahat,2)-2*cetaCorParams[1]/pow(cetahat,3);
            cetahat -= elPrime / elDPrime;
        }
        cetaCor = cetahat;
	}else{
        cetaCor = LASTAll[move2++];
	}
    
    //#############################################SAMPLER
    for (sw = 0; sw < sweeps; sw++){
        if (sw==0) Rprintf("%i %s \n",sw+1, "posterior sample...");
        if (((sw+1) % 1000)==0) Rprintf("%i %s \n",sw+1, "posterior samples...");

	    modf((sw+1)/batchL,&WB);
	    WB += contParams[2];
 
        // - 1 - gamma
        //Rprintf("%i %s \n",sw,"gamma");        
		SPC = ScalcMultLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPC,Y,X,gamma,RiAll,intime,ceta,&Q2); 		
		for (l = 0; l < p; l++){
		    for (j = 0; j < NG; j++){
		        for (k = 0; k < vecLG[j]; k++)
		            vecGamma[j][k] = gamma[l][cusumVecLG[j]+k];
		        gsl_ran_multinomial(r,maxBSG,1,blockSizeProbG,vecBSG);
                blockSize = 0;
	            while(vecBSG[blockSize]==0) blockSize++;
                blockSize += 1;
                nBlocks = ceil((double)vecLG[j]/blockSize);
                gsl_ran_shuffle(r,indexG[j],vecLG[j],sizeof(int));
	            for (block = 0; block < nBlocks; block++){
                    s = gsl_ran_flat(r,1.0,100000);
	                proposeBlockInd(s,vecGamma[j],vecLG[j],block,blockSize,indexG[j],cmu[l][j],dmu[l][j],vecGammaP[j]);
	                for (k = 0; k < vecLG[j]; k++)
	                    gamma[l][cusumVecLG[j]+k] = vecGammaP[j][k];
	                NPJ = 0;
                    NCJ = 0;
                    for (k = 0; k < vecLG[j]; k++){
                        NPJ += vecGammaP[j][k];
                        NCJ += vecGamma[j][k];
					}
					NgammaP = Ngamma - NCJ + NPJ;
	                SPP = ScalcMultLong(m,p,tol,LG,NgammaP,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPC,Y,X,gamma,RiAll,intime,ceta,&Q2); 	              
                    Acp = exp((SPC-SPP)/2)*pow(ceta+1,0.5*(NCJ-NPJ));
                    //Rprintf("%s %i %i %i %i %f %f %f %f %i %i\n","gamma: ",sw,l,j,block,Acp,SPC,SPP,ceta,NCJ,NPJ);
                    //fprintf(out_file24, "%s %i %f ","gamma: ", sw, Acp);
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
	                    
	    // - 2 - c_beta
	    //Rprintf("%i %s \n",sw,"c_beta");
        if ((sw % batchL)==0 && WB > 0){
	        if (acceptCeta > 0.25 && tuneCb[0] < GUL) tuneCb[0] += MIN(0.01,1/sqrt(WB)) * tuneCb[0];
	        if (acceptCeta <= 0.25 && tuneCb[0] > GLL) tuneCb[0] -= MIN(0.01,1/sqrt(WB)) * tuneCb[0];
	        acceptCeta = 0.0;
	        if (tuneCb[0] < GLL) tuneCb[0] = GLL;
	        if (tuneCb[0] > GUL) tuneCb[0] = GUL;
	        fprintf(out_file25,"%f ",tuneCb[0]);
        }
        cetahat = 1;
        SPC = ScalcMultLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPC,Y,X,gamma,RiAll,intime,ceta,&Q2); 
        elPrime = 99.9;
        while(elPrime > 0.000000001 || -elPrime > 0.000000001){
            Sprime = -Q2/pow(cetahat+1,2);
            SDprime = 2*Q2/pow(cetahat+1,3);
            elPrime = -0.5*(Ngamma+p)/(cetahat+1)-Sprime/2-(alphaeta+1)/cetahat+betaeta/pow(cetahat,2);
            elDPrime = 0.5*(Ngamma+p)/(pow(cetahat+1,2))-SDprime/2+(alphaeta+1)/pow(cetahat,2)-2*betaeta/pow(cetahat,3);
            cetahat -= elPrime / elDPrime;
        }
//Rprintf("%s %f %f %f %f %f %f %f \n","ceta 1: ",cetahat,elPrime,SPC,Sprime,SDprime,Q2,sigma2zk[0]);
        cetaP = cetahat + gsl_ran_gaussian(r,sqrt(-tuneCb[0]/elDPrime));
	    while(cetaP < 0) cetaP = cetahat + gsl_ran_gaussian(r,sqrt(-tuneCb[0]/elDPrime));
	    SPP = ScalcMultLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPC,Y,X,gamma,RiAll,intime,cetaP,&Q2); 	    
	    Acp = exp(-0.5*(Ngamma+p)*(log(cetaP+1)-log(ceta+1))+(-SPP+SPC)/2-(alphaeta+1)*(log(cetaP)-log(ceta))
              + betaeta*(1/ceta-1/cetaP))*
                gsl_ran_gaussian_pdf(ceta-cetahat,sqrt(-tuneCb[0]/elDPrime))/
                gsl_ran_gaussian_pdf(cetaP-cetahat,sqrt(-tuneCb[0]/elDPrime));
	    unifRV = gsl_ran_flat(r,0.0,1.0);
	    //fprintf(out_file24, "%s %i %f ","cbeta: ", sw, Acp);
        if (Acp > unifRV || sw < 5){
            ceta = cetaP;
            SPC = SPP;
	        acceptCeta += 1/((double)batchL);
	    }
//Rprintf("%s %i %f %f %f %f %f %f %f %f %f %f %i %f %f\n","ceta 2: ",
//sw,Acp,unifRV,SPC,SPP,ceta,cetahat,cetaP,Q2,alphaeta,betaeta,Ngamma,tuneCb[0],elDPrime);
	    
	    // - 3 - Update delta and alpha 
        //Rprintf("%i %s \n",sw,"delta & alpha");
        subMeanEta = gsl_vector_subvector(meanEta,0,Ngamma+p);
        subVarEta = gsl_matrix_submatrix(varEta,0,0,Ngamma+p,Ngamma+p);                    
        postMeanVarEtaLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPC,Y,X,gamma,
                           RiAll,intime,ceta,&subMeanEta.vector,&subVarEta.matrix);                            
        //cSqRes2(p,N,LG,gamma,Ngamma,X,&subMeanEta.vector,Y,sqResC);       
        cRes(p,N,LG,gamma,Ngamma,X,&subMeanEta.vector,Y,ResC,BaseXg);
	    cResCheck(p,m,N,niMax,niVec,cusumniVec,sigma2ij,dimLPC,LPC,ResC,ResCheck,0);               
        for (l2 = 0; l2 < (N*p); l2++) 
            sqResC[l2] = pow(ResCheck[l2],2);                
        for (l = 0; l < p; l++){
            for (j = 0; j < ND; j++){
		        if ((sw % batchL)==0 && WB > 0){
	                if (acceptAlpha[l][j] > 0.25 && tuneAlpha[ND*l+j] < HUL) tuneAlpha[ND*l+j] += MIN(0.01,1/sqrt(WB)) * tuneAlpha[ND*l+j];
	                if (acceptAlpha[l][j] <= 0.25 && tuneAlpha[ND*l+j] > HLL) tuneAlpha[ND*l+j] -= MIN(0.01,1/sqrt(WB)) * tuneAlpha[ND*l+j];
	                acceptAlpha[l][j] = 0.0;
	                if (tuneAlpha[ND*l+j] < HLL) tuneAlpha[ND*l+j] = HLL;
	                if (tuneAlpha[ND*l+j] > HUL) tuneAlpha[ND*l+j] = HUL;
	                fprintf(out_file25,"%f ",tuneAlpha[ND*l+j]);
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
                            DeltaAlphaHatExp(N,p,l,tol,LPV,sqResC,vecDeltaP[j],NPJ,cusumVecLD[j],cusumVecLD[j+1],
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
                        for (t = 0; t < niVec[i]; t++){
                            LPVP[cusumniVec[i]+t][l] = LPV[cusumniVec[i]+t][l];
                            for (k = cusumVecLD[j]; k < cusumVecLD[j+1]; k++){
                                LPVP[cusumniVec[i]+t][l] += (alphaPD[k]-alpha[l][k])*Z[N+k*N+cusumniVec[i]+t];
                                //if (1==0 && sw==0 && l==0 && j==0 && block==0)
                                //Rprintf("%s %i %i %i %i %f %f %f %f \n","new:", sw,l,j,block,alphaPD[k],alpha[l][k],Z[m+k*m+i],LPVP[i][l]);
						    }
                            sigma2ijP[cusumniVec[i]+t][l] = sigma2zk[l]*exp(LPVP[cusumniVec[i]+t][l]);
						}
                    }                   
                    SPP = ScalcMultLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ijP,dimLPC,LPC,Y,X,gamma,RiAll,intime,ceta,&Q2);                     
                    QFD = 0.0;
                    for (k = cusumVecLD[j]; k < cusumVecLD[j+1]; k++)
                        QFD += pow(alphaPD[k],2)-pow(alpha[l][k],2);
                    detR = 0.0;
                    for (i = 0; i < m; i++)
                        for (t = 0; t < niVec[i]; t++)
                            detR += LPV[cusumniVec[i]+t][l] - LPVP[cusumniVec[i]+t][l];
                    detR *= 0.5;                    
	                //probability of reverse direction
                    postMeanVarEtaLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ijP,dimLPC,LPC,Y,X,gamma,
                                       RiAll,intime,ceta,&subMeanEta.vector,&subVarEta.matrix);                                        
                    
                    //cSqRes2(p,N,LG,gamma,Ngamma,X,&subMeanEta.vector,Y,sqResP);                                                                    
                    cRes(p,N,LG,gamma,Ngamma,X,&subMeanEta.vector,Y,ResC,BaseXg);
	                cResCheck(p,m,N,niMax,niVec,cusumniVec,sigma2ij,dimLPC,LPC,ResC,ResCheck,0);               
                    for (l2 = 0; l2 < (N*p); l2++) 
                        sqResP[l2] = pow(ResCheck[l2],2);        
                    
                    if (NCJ > 0){
	                    subAlphaHat = gsl_vector_subvector(alphaHat,0,NCJ);
                        subD = gsl_matrix_submatrix(D,0,0,NCJ,NCJ);
                        if (isDz[j] == 0){ 
                            DeltaAlphaHatExp(N,p,l,tol,LPVP,sqResP,vecDelta[j],NCJ,cusumVecLD[j],cusumVecLD[j+1],
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
                    
                    Acp = exp(detR+(-SPP+SPC)/2+logMVNormC-logMVNormP-QFD/(2*calpha[l]))*
                              pow(2*M_PI*calpha[l],0.5*(NCJ-NPJ));
                    unifRV = gsl_ran_flat(r,0.0,1.0);
	                //fprintf(out_file24, "%s %i %i %i %i | %f %f | %f %f %f %f | %f %f %f %f \n","delta: ",
                    //sw,l,j,block,Acp,unifRV,SPC,SPP,logMVNormC,logMVNormP,QFD,calpha[l],detR,pow(2*M_PI*calpha[l],0.5*(NCJ-NPJ)));
	                //fprintf(out_file24, "%s %i %f ","delta: ", sw, Acp);
	                if (Acp > unifRV){
                        if (NPJ > 0) acceptAlpha[l][j] += 1/((double)batchL);
                        for (k = 0; k < vecLD[j]; k++){
                            delta[l][cusumVecLD[j]+k] = vecDeltaP[j][k];
                            vecDelta[j][k] = vecDeltaP[j][k];
                            alpha[l][cusumVecLD[j]+k] = alphaPD[cusumVecLD[j]+k];
                        }
                        for (i = 0; i < m; i++){
                            for (t = 0; t < niVec[i]; t++){
                                LPV[cusumniVec[i]+t][l] = LPVP[cusumniVec[i]+t][l];
                                sigma2ij[cusumniVec[i]+t][l] = sigma2ijP[cusumniVec[i]+t][l];
							}
						}
				        for (i = 0; i < (N*p); i++)
                            sqResC[i] = sqResP[i];
                        SPC = SPP;                        
                    }
                    else{
				    	if (j == (ND-1)){ //probably && block == (nBlocks-1) && LPVL not needed
				    	    for (i = 0; i < m; i++){
                                for (t = 0; t < niVec[i]; t++){
                                    LPVP[cusumniVec[i]+t][l] = LPV[cusumniVec[i]+t][l];
                                    sigma2ijP[cusumniVec[i]+t][l] = sigma2ij[cusumniVec[i]+t][l];
								}
							}
						}
				    }
                }
	        }
	    }
	    //Rprintf("%s %i %i %i %f %f %f %f \n","after alpha: ",p,m,LG,ceta,sigma2zk[0]*SPC,LPV[0][0],Ytilde[0]*sqrt(sigma2zk[0]));

	    // - 4 - sigma2_k, k=1,...,p
        //Rprintf("%i %s \n",sw,"sigma2_k");
		move = 0;
		for (k = 0; k < p; k++){
            if ((sw % batchL)==0 && WB > 0){
	            if (acceptSZK[k] > 0.25 && tuneSigma2[k] < SZKUL) tuneSigma2[k] += MIN(0.01,1/sqrt(WB)) * tuneSigma2[k];
	            if (acceptSZK[k] <= 0.25 && tuneSigma2[k] > SZKLL) tuneSigma2[k] -= MIN(0.01,1/sqrt(WB)) * tuneSigma2[k];
	            acceptSZK[k] = 0.0;
	            if (tuneSigma2[k] < SZKLL) tuneSigma2[k] = SZKLL;
	            if (tuneSigma2[k] > SZKUL) tuneSigma2[k] = SZKUL;
	            fprintf(out_file25,"%f ",tuneSigma2[k]);
            }
            sigma2P = sigma2zk[k] + gsl_ran_gaussian(r,sqrt(tuneSigma2[k]));
            while (sigma2P <= 0) sigma2P = sigma2zk[k] + gsl_ran_gaussian(r,sqrt(tuneSigma2[k]));
            for (i = 0; i < m; i++)
                for (t = 0; t < niVec[i]; t++)
                    sigma2ijP[cusumniVec[i]+t][k] = sigma2P*exp(LPV[cusumniVec[i]+t][k]);                                          
            SPP = ScalcMultLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ijP,dimLPC,LPC,Y,X,gamma,RiAll,intime,ceta,&Q2);             
            detR = 0.5*N*(log(sigma2zk[k])-log(sigma2P));
            //Prior Ratio
            if (HNszk[k]==1){
                phi2G = szkParams[move++];
                priorLogR = 0.5*(sigma2zk[k]-sigma2P)/phi2G;
			}else{
				alphaG = szkParams[move++];
				betaG = szkParams[move++];
                priorLogR = (alphaG+1)*(log(sigma2zk[k])-log(sigma2P)) + betaG*(1/sigma2zk[k]-1/sigma2P);
			}
            Acp = exp(detR+(SPC-SPP)/2+priorLogR);
            unifRV = gsl_ran_flat(r,0.0,1.0);
	        //Rprintf("%s %i %i %f %f %f %f %f %f %f %f %f %i %f %f %f \n","sigma2: ",sw,k,Acp,sigma2zk[k],sigma2P,
	        //                    detR+(SPC-SPP)/2+priorLogR,priorLogR,(SPC-SPP)/2,detR,SPC,SPP,HNszk[k],phi2G,alphaG,betaG);
	        //fprintf(out_file24, "%s %i %f ","sigma2: ", sw, Acp);
            if (Acp > unifRV){
                sigma2zk[k] = sigma2P;	            	            
	            for (i = 0; i < m; i++)
                    for (t = 0; t < niVec[i]; t++)
                        sigma2ij[cusumniVec[i]+t][k] = sigma2ijP[cusumniVec[i]+t][k];                       					                
                SPC = SPP;
                acceptSZK[k] += 1/((double)batchL);
            }
            else{
			    for (i = 0; i < m; i++)
                    for (t = 0; t < niVec[i]; t++)
                        sigma2ijP[cusumniVec[i]+t][k] = sigma2ij[cusumniVec[i]+t][k];                                        				
			}
	    }	    	   
	    
	    // - 5 - c_alpha 
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
		        //Rprintf("%s %i %i %f %f %i %f %f %f \n","calpha G: ",sw,k,alphaG,betaG,NPJ,QFC,calpha[k],alphaG+0.5*NPJ);
		    }else if (HNca[k]==1){
                if ((sw % batchL)==0 && WB > 0){
	                if (acceptCa[k] > 0.25 && tuneCa[k] < FCaUL) tuneCa[k] += MIN(0.01,1/sqrt(WB)) * tuneCa[k];
	                if (acceptCa[k] <= 0.25 && tuneCa[k] > FCaLL) tuneCa[k] -= MIN(0.01,1/sqrt(WB)) * tuneCa[k];
	                acceptCa[k] = 0.0;
	                if (tuneCa[k] < FCaLL) tuneCa[k] = FCaLL;
	                if (tuneCa[k] > FCaUL) tuneCa[k] = FCaUL;
	                fprintf(out_file25,"%f ",tuneCa[k]);
                }
                phi2G = calphaParams[move++];
                calphaP = calpha[k] + gsl_ran_gaussian(r,sqrt(tuneCa[k]));
                while (calphaP <= 0) calphaP = calpha[k] + gsl_ran_gaussian(r,sqrt(tuneCa[k]));
                Acp = exp(-0.5*NPJ*(log(calphaP)-log(calpha[k])) + (QFC/2)*(1/calpha[k]-1/calphaP) +
                      (calpha[k]-calphaP)/(2*phi2G));
                //fprintf(out_file24, "%s %i %f ","calpha: ", sw, Acp);      
	            unifRV = gsl_ran_flat(r,0.0,1.0);
	            //Rprintf("%s %i %f %i %f %f %f %f %f \n","calpha HN: ",sw,Acp,NPJ,calphaP,calpha[k],QFC,phi2G,tuneCa[k]);
                if (Acp > unifRV){
                    calpha[k] = calphaP;
	                acceptCa[k] += 1/((double)batchL);
                }
	        }
	    }
	    
	    // - 6 - eta 
        //Rprintf("%i %s \n",sw,"eta");
        subMeanEta = gsl_vector_subvector(meanEta,0,Ngamma+p);
        subVarEta = gsl_matrix_submatrix(varEta,0,0,Ngamma+p,Ngamma+p);
        postMeanVarEtaLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPC,Y,X,gamma,
                           RiAll,intime,ceta,&subMeanEta.vector,&subVarEta.matrix);            
        vecEta = gsl_vector_view_array(eta,Ngamma+p);
        s = gsl_ran_flat(r,1.0,100000);
        sampleMN(s,Ngamma+p,&vecEta.vector,&subMeanEta.vector,&subVarEta.matrix,tol);
	     
	    // - 7 - Psi & Ksi
        //Rprintf("%i %s \n",sw,"psi & ksi");        
        cRes(p,N,LG,gamma,Ngamma,X,&vecEta.vector,Y,ResC,BaseXg);                         
        start = 0;                
        for (k = 0; k < p; k++){
            for (l = 0; l < p; l++){
                for (j = 0; j < NK; j++){		        		            		            		            
		            move = 0;
                    for (k2 = 0; k2 < p; k2++){
                        for (l2 = 0; l2 < p; l2++){
                            for (j2 = 0; j2 < NK; j2++)
                                for (i = 0; i < vecLK[j2]; i++)
                                    if (ksi[k2][l2][cusumVecLK[j2]+i]==1 && !(k2==k && l2==l && j2==j)) 
                                        W[move++] = psi[k2][l2][cusumVecLK[j2]+i];			               		                    	                    		            
						}
					}		            
		            for (k2 = 0; k2 < vecLK[j]; k2++)
                        vecKsi[j][k2] = ksi[k][l][cusumVecLK[j] + k2];		            		            		            
		            gsl_ran_multinomial(r,maxBSK,1,blockSizeProbK,vecBSK);
                    blockSize = 0;
                    while(vecBSK[blockSize]==0) blockSize++;
                    blockSize += 1;
                    nBlocks = ceil((double)vecLK[j]/blockSize);
                    gsl_ran_shuffle(r,indexK[j],vecLK[j],sizeof(int));                    
                    for (block = 0; block < nBlocks; block++){                                            
                        s = gsl_ran_flat(r,1.0,100000);
                        proposeBlockInd(s,vecKsi[j],vecLK[j],block,blockSize,indexK[j],cphi[k][l][j],dphi[k][l][j],vecKsiP[j]);                                            
                        for (k2 = 0; k2 < vecLK[j]; k2++)
	                        ksiP[k][l][cusumVecLK[j] + k2] = vecKsiP[j][k2];                        
                        NPJ = 0;
                        NCJ = 0;                        
                        for (k2 = 0; k2 < vecLK[j]; k2++){
                            NPJ += vecKsiP[j][k2];
                            NCJ += vecKsi[j][k2];
                        }                                                 
                        NksiP = Nksi - NCJ + NPJ;                                                                                                 
                        /*
                        Rprintf("%s %i | %i %i %i | %i %i %i | %i %f %i | %i %i , %i %i | %i\n","Nksi's:", 
                                 sw, 
                                      k, l, p, 
                                             j, NK, vecLK[j], 
                                                             block, nBlocks, blockSize, 
                                                                                       Nksi, NCJ, NksiP, NPJ,
                                                                                                             start);                                                                                                                       
                        */
                        if (NPJ > 0){                                                                                                                         
                            subMeanPsi = gsl_vector_subvector(meanPsi,0,NksiP);
                            subVarPsi = gsl_matrix_submatrix(varPsi,0,0,NksiP,NksiP);                                                                 
                            postMeanVarPSI(m,p,tol,niVec,cusumniVec,LK,NksiP,ResC,ksiP,cusumC,C,N,sigma2ij,
                                           RiAll,intime,cPsi,&subMeanPsi.vector,&subVarPsi.matrix);                                                                                                                             
                            //puts("mean psi");
                            //print_vector(&subMeanPsi.vector);
                            //puts("var psi");
                            //print_matrix(&subVarPsi.matrix);                            
                            subAlphaHat = gsl_vector_subvector(alphaHat,0,NPJ); 
                            subD = gsl_matrix_submatrix(D,0,0,NPJ,NPJ);                                                                                                                                
                            //Rprintf("%i %s %i %i \n",sw, "MN start end", start, start+NPJ-1);   
                            MNCond(tol,start,start+NPJ-1,&subMeanPsi.vector,&subVarPsi.matrix,W,&subAlphaHat.vector,&subD.matrix);                                                                                
                            subAlphaP = gsl_vector_subvector(alphaP,0,NPJ);                            
                            s = gsl_ran_flat(r,1.0,100000);
                            //puts("cond mean");
                            //print_vector(&subAlphaHat.vector);
                            //puts("cond var");  
                            //print_matrix(&subD.matrix);    
                            sampleMN(s,NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);                        
                            logMVNormP = logMVNormalpdf(NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);                        
                        }else logMVNormP = 0.0;                                                                                                                                         
                        move = 0;
                        for (k2 = 0; k2 < vecLK[j]; k2++){
                            if (vecKsiP[j][k2]==1) {psiPK[cusumVecLK[j]+k2] = gsl_vector_get(&subAlphaP.vector,move++);						
					        }else{psiPK[cusumVecLK[j]+k2] = 0;}
                        }                        
                        move = 0;  
                        for (i = 0; i < m; i++){	
		                    for (j2 = 0; j2 < (niVec[i]*(niVec[i]-1)/2); j2++){			
                                LPCP[i][j2][k][l] = LPC[i][j2][k][l];
                                    for (k2 = cusumVecLK[j]; k2 < cusumVecLK[j+1]; k2++)//{
                                        LPCP[i][j2][k][l] += (psiPK[k2]-psi[k][l][k2])*C[k2+move*LK];
                                        //Rprintf("%i | %i %i %i | %i %i %i | %f %f %f \n",sw,k,l,j,i,j2,k2,psiPK[k2],psi[k][l][k2],C[k2+move*LK]);}				                   
	                            move++;  
		                    }
                        }                             
                        SPP = ScalcMultLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPCP,Y,X,gamma,RiAll,intime,ceta,&Q2);
                        QFD = 0.0;
                        for (k2 = cusumVecLK[j]; k2 < cusumVecLK[j+1]; k2++)
                            QFD += pow(psiPK[k2],2)-pow(psi[k][l][k2],2);	                    
	                    //probability of reverse direction                        
                        postMeanVarEtaLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPCP,Y,X,gamma,
                                           RiAll,intime,ceta,&subMeanEta.vector,&subVarEta.matrix);                    
                        vecEta = gsl_vector_view_array(eta2,Ngamma+p);
                        s = gsl_ran_flat(r,1.0,100000);
                        sampleMN(s,Ngamma+p,&vecEta.vector,&subMeanEta.vector,&subVarEta.matrix,tol);                         
                        //Rprintf("%s %i \n","eta2 ",sw);
                        //print_vector(&vecEta.vector);                                                                       
                        cRes(p,N,LG,gamma,Ngamma,X,&vecEta.vector,Y,ResP,BaseXg);                                                                                
                        //Rprintf("%i %s \n",sw,"3");
                        if (NCJ > 0){							
							subMeanPsi = gsl_vector_subvector(meanPsi,0,Nksi);
                            subVarPsi = gsl_matrix_submatrix(varPsi,0,0,Nksi,Nksi);        
                            postMeanVarPSI(m,p,tol,niVec,cusumniVec,LK,Nksi,ResP,ksi,cusumC,C,N,sigma2ij,
                                           RiAll,intime,cPsi,&subMeanPsi.vector,&subVarPsi.matrix);
							subAlphaHat = gsl_vector_subvector(alphaHat,0,NCJ); 
                            subD = gsl_matrix_submatrix(D,0,0,NCJ,NCJ);                                                 
                            MNCond(tol,start,start+NCJ-1,&subMeanPsi.vector,&subVarPsi.matrix,W,&subAlphaHat.vector,&subD.matrix);                            
                            move=0;
                            for (k2 = 0; k2 < vecLK[j]; k2++)
                                if (vecKsi[j][k2]==1) BaseSubAlpha[move++] = psi[k][l][cusumVecLK[j]+k2];
                            subAlphaP = gsl_vector_view_array(BaseSubAlpha,NCJ);                                                                                     
                            logMVNormC = logMVNormalpdf(NCJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);	                    	                    	                    	                                            
	                    }else logMVNormC = 0.0;	                    
                        //Rprintf("%i %s \n",sw,"4");
                        Acp = exp((-SPP+SPC)/2+logMVNormC-logMVNormP-QFD/(2*cPsi[k*p+l]))*pow(2*M_PI*cPsi[k*p+l],0.5*(NCJ-NPJ));
                        unifRV = gsl_ran_flat(r,0.0,1.0);
                        //fprintf(out_file24, "%s %i %f ","ksi & psi : ", sw, Acp);                                               
	                    if (Acp > unifRV){               					
							if (block == (nBlocks-1)) start += NPJ;									             
                            for (k2 = 0; k2 < vecLK[j]; k2++){                                
                                vecKsi[j][k2] = vecKsiP[j][k2];
                                psi[k][l][cusumVecLK[j]+k2] = psiPK[cusumVecLK[j]+k2];
                                ksi[k][l][cusumVecLK[j] + k2] = ksiP[k][l][cusumVecLK[j] + k2];
                            }                                                  
                            for (i = 0; i < m; i++)	
                                for (j2 = 0; j2 < (niVec[i]*(niVec[i]-1)/2); j2++)			
                                    LPC[i][j2][k][l] = LPCP[i][j2][k][l];                           
				            for (i = 0; i < (N*p); i++)
                                ResC[i] = ResP[i];
                            for (i = 0; i < (Ngamma+p); i++)
                                eta[i] = eta2[i];                                                                          
                            SPC = SPP;   
                            Nksi = NksiP;                                           			                                                                           
                        }
                        else{
							if (block == (nBlocks-1)){ 
							    start += NCJ;
							    for (k2 = 0; k2 < vecLK[j]; k2++)
						            ksiP[k][l][cusumVecLK[j] + k2] = ksi[k][l][cusumVecLK[j] + k2];
							}
				    	    if ((j == (NK-1)) && (block == (nBlocks-1))){				    	        
				    	        for (i = 0; i < m; i++)	
		                            for (j2 = 0; j2 < (niVec[i]*(niVec[i]-1)/2); j2++)		
                                        LPCP[i][j2][k][l] = LPC[i][j2][k][l];                          							    
						    }						    	                           
				        }
                    }
	            }
	        }
		}
	    		
	    //puts("end");
	    
	    //print_vector(&vecEta.vector);                                               		
	    		
		// - 8 - c_psi
        //Rprintf("%i %s \n",sw,"c_psi");
        move = 0;
        for (k = 0; k < p; k++){
            for (l = 0; l < p; l++){
                NPJ = 0;
                QFC = 0;                
                for (j = 0; j < LK; j++){
                    NPJ += ksi[k][l][j];
                    QFC += pow(psi[k][l][j],2);
                }                
                if (HNcpsi[k*p+l]==0){
		            alphaG = cpsiParams[move++]; 
		            betaG = cpsiParams[move++];		    
		            cPsi[k*p+l] = 1/gsl_ran_gamma(r,alphaG+0.5*NPJ,1/(betaG+0.5*QFC));		        
		        }else if (HNcpsi[k*p+l]==1){
                    if ((sw % batchL)==0 && WB > 0){
	                    if (acceptCpsi[k*p+l] > 0.25 && tuneCpsi[k*p+l] < FCpsiUL) tuneCpsi[k*p+l] += MIN(0.01,1/sqrt(WB)) * tuneCpsi[k*p+l];
	                    if (acceptCpsi[k*p+l] <= 0.25 && tuneCpsi[k*p+l] > FCpsiLL) tuneCpsi[k*p+l] -= MIN(0.01,1/sqrt(WB)) * tuneCpsi[k*p+l];
	                    acceptCpsi[k*p+l] = 0.0;
	                    if (tuneCpsi[k*p+l] < FCpsiLL) tuneCpsi[k*p+l] = FCpsiLL;
	                    if (tuneCpsi[k*p+l] > FCpsiUL) tuneCpsi[k*p+l] = FCpsiUL;
	                    fprintf(out_file25,"%f ",tuneCpsi[k*p+l]);
                    }
                    phi2G = cpsiParams[move++];
                    cPsiP = cPsi[k*p+l] + gsl_ran_gaussian(r,sqrt(tuneCpsi[k*p+l]));
                    while (cPsiP <= 0) cPsiP = cPsi[k*p+l] + gsl_ran_gaussian(r,sqrt(tuneCpsi[k*p+l]));
                    Acp = exp(-0.5*NPJ*(log(cPsiP)-log(cPsi[k*p+l])) + (QFC/2)*(1/cPsi[k*p+l]-1/cPsiP) + (cPsi[k*p+l]-cPsiP)/(2*phi2G));
	                unifRV = gsl_ran_flat(r,0.0,1.0);
                    //fprintf(out_file24, "%s %i %f ","cpsi: ", sw, Acp);
                    if (Acp > unifRV){
                        cPsi[k*p+l] = cPsiP;
	                    acceptCpsi[k*p+l] += 1/((double)batchL);
                    }
	            }
			}
		}
		
		// - 9 - R_t, t=1,2,...,T
        //Rprintf("%i %s \n",sw,"R_t");
        
        cResCheck(p,m,N,niMax,niVec,cusumniVec,sigma2ij,dimLPC,LPC,ResC,ResCheck,1); //update only when ResC is updated
	
	    //for (k = 0; k < (N*p); k++) Rprintf("%i %i %i %f %i \n",sw,t,k,ResCheck[k],intime2[k]);
	    //for (k = 0; k < (N*p); k++) Rprintf("%i %i \n", k, intime2[k]);
        
		for (t = 0; t < LUT; t++){	        
	        
	        if ((sw % batchL)==0 && WB > 0){ 
	            if (acceptRt[t] > 0.25 && tuneR[t] > RtLL) tuneR[t] -= MIN(0.01,1/sqrt(WB)) * tuneR[t]; 
	            if (acceptRt[t] <= 0.25 && tuneR[t] < RtUL) tuneR[t] += MIN(0.01,1/sqrt(WB)) * tuneR[t];	            
		        acceptRt[t] = 0.0;
		        if (tuneR[t] < RtLL) tuneR[t] = RtLL;
	            if (tuneR[t] > RtUL) tuneR[t] = RtUL;   	            	            
	            fprintf(out_file25,"%f ",tuneR[t]);
	        }
	        
	        
	        
	        if (t > 0 && (Acp > unifRV)){//otherwise ResC remains the same
	            postMeanVarEtaLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPC,Y,X,gamma,
                                   RiAll,intime,ceta,&subMeanEta.vector,&subVarEta.matrix);            
                vecEta = gsl_vector_view_array(eta,Ngamma+p);
                s = gsl_ran_flat(r,1.0,100000);
                sampleMN(s,Ngamma+p,&vecEta.vector,&subMeanEta.vector,&subVarEta.matrix,tol);	
	            cRes(p,N,LG,gamma,Ngamma,X,&vecEta.vector,Y,ResC,BaseXg);
	            cResCheck(p,m,N,niMax,niVec,cusumniVec,sigma2ij,dimLPC,LPC,ResC,ResCheck,1);
			}
			
			

			
            computeStStar(ResCheck,intime2,N*p,t,p,St);

            for (k = 0; k < p*p; k++){
                baseRt[k] = Rt[t][k];
                baseDt[k] = Dt[t][k];
                baseEt[k] = Et[t][k];
            }
            
            EtC = gsl_matrix_view_array(baseEt,p,p);
            DtC = gsl_matrix_view_array(baseDt,p,p);
            RtC = gsl_matrix_view_array(baseRt,p,p);           
	    	    
	        gsl_matrix_memcpy(PSIt,&EtC.matrix);
//Rprintf("%i %i %s \n",sw,t,"rwish - check");
//print_matrix(St);
//print_matrix(PSIt);
            gsl_matrix_scale(PSIt,tuneR[t]-p-1);      
	        gsl_matrix_add(PSIt,St);	    
            gsl_matrix_memcpy(CopyPSIt,PSIt);            
            ginv(p,tol,PSIt);            
            s = gsl_ran_flat(r,1.0,100000);
//print_matrix(PSIt);
//Rprintf("%s \n","before");            
            rwish(s,p,FUT[t]+tuneR[t],PSIt,gEt);           	    	                
//Rprintf("%s \n","after");
            ginv(p,tol,gEt);
            
            decomposeEtoDS(p,0,gEt,gDt,gRt);   
            priorLogR = 0.0;
            move = 0;
            for (k = 0; k < (p-1); k++){
                for (l = (k+1); l < p; l++){
					priorLogR += pow(FisherTr(gsl_matrix_get(gRt,k,l),FT[0])-theta[move+t*d],2) -
                                     pow(FisherTr(rt[move+t*d],FT[0])-theta[move+t*d],2);
                    move++;
                } 
            } 
            priorLogR *= -1/(2*tau2);
            
            if (FT[0]==1){
                move = 0;
                for (k = 0; k < (p-1); k++){
                    for (l = (k+1); l < p; l++){
                        priorLogR += log(1+rt[move+t*d]) + log(1-rt[move+t*d]) - 
                                     log(1+gsl_matrix_get(gRt,k,l)) - log(1-gsl_matrix_get(gRt,k,l));
                        move++;
                    } 
                } 
		    }
 
            logLikP = logPropPdfDR(gDt,gRt,St,gDt,p,0,FUT[t],0,1);
            logLikC = logPropPdfDR(gDt,&RtC.matrix,St,gDt,p,0,FUT[t],0,1);            
	 
	        gsl_matrix_memcpy(PSIt,&EtC.matrix);            
            gsl_matrix_scale(PSIt,tuneR[t]-p-1);	    
	        logPropDRP = logPropPdfDR(gDt,gEt,CopyPSIt,PSIt,p,p-1,FUT[t]+tuneR[t]+p+1,tuneR[t],1);
            
	        gsl_matrix_memcpy(CopyPSIt,gEt);            
            gsl_matrix_scale(CopyPSIt,tuneR[t]-p-1);
	        gsl_matrix_add(CopyPSIt,St);
	    
	        gsl_matrix_memcpy(PSIt,gEt);            
            gsl_matrix_scale(PSIt,tuneR[t]-p-1);	  
	    
	        logPropDRC = logPropPdfDR(&DtC.matrix,&EtC.matrix,CopyPSIt,PSIt,p,p-1,FUT[t]+tuneR[t]+p+1,tuneR[t],1);

	        logAcp = logLikP - logLikC + priorLogR + logPropDRC - logPropDRP;

            Acp = exp(logAcp);
            unifRV = gsl_ran_flat(r,0.0,1.0);                        
            //fprintf(out_file24, "%s %i %f ","Rt: ", sw, Acp);          
            if (Acp > unifRV || sw < 5){
                acceptRt[t] += 1/((double)batchL);  
	            for (k = 0; k < p; k++){
                    for (l = 0; l < p; l++){
                        Rt[t][k*p+l] = gsl_matrix_get(gRt,k,l);
                        Et[t][k*p+l] = gsl_matrix_get(gEt,k,l);
                        Dt[t][k*p+l] = gsl_matrix_get(gDt,k,l);
                    } 
                }
                move = 0;
                for (k = 0; k < (p-1); k++){
                    for (l = (k+1); l < p; l++){
                        rt[move+t*d] = Rt[t][k*p+l];
                        move++;
                    }
                }
                subD = gsl_matrix_submatrix(RiAll,0,t*p,p,p);
                gsl_matrix_memcpy(&subD.matrix,gRt);            
                ginv(p,tol,&subD.matrix);                      
            }
        } 
        
        // - 10 - Theta       
        //Rprintf("%i %s \n",sw,"Theta");                 
        for (h = 0; h < H; h++){			 
			if (nmembers[h]>0){                                
                vecTheta = gsl_vector_view_array(vizTheta,LUT*nmembers[h]);           
                subA = gsl_matrix_submatrix(A,0,0,LUT*nmembers[h],LUT*nmembers[h]);           
                iCovPostThetah(LUT,d,H,h,tol,LGc,gammaCor,compAlloc,nmembers[h],NgammaCor,cetaCor,sigma2cor,
                               tau2,Xc,LPVcor,&subA.matrix);                         
                move = 0; 
                for (t = 0; t < LUT; t++)
                    for (i = 0; i < d; i++)
                        if (compAlloc[i]==h) vizrt[move++] = FisherTr(rt[i+t*d],FT[0]);                                 
                Ar = gsl_vector_view_array(vizrt,LUT*nmembers[h]);            
                subMutheta = gsl_vector_subvector(mutheta,0,LUT*nmembers[h]);             
                gsl_blas_dgemv(CblasNoTrans,1.0,&subA.matrix,&Ar.vector,0.0,&subMutheta.vector);            
                gsl_vector_scale (&subMutheta.vector,1/tau2);            
                s = gsl_ran_flat(r,1.0,100000);             
                sampleMN(s,LUT*nmembers[h],&vecTheta.vector,&subMutheta.vector,&subA.matrix,tol);
                move2 = 0;
                for (t = 0; t < LUT; t++)
                    for (i = 0; i < d; i++) 
                        if (compAlloc[i]==h) theta[i+t*d] = vizTheta[move2++]; 
            }
        }
        
        for (t = 0; t < LUT; t++)
            for (i = 0; i < d; i++)
                thetaTilde[i+t*d] = exp(-LPVcor[t]/2)*theta[i+t*d];                
        
        // - 11 - Update Cor gamma 
        //Rprintf("%i %s \n",sw,"cor gamma");
        for (h = 0; h < H; h++){
		    for (j = 0; j < NGc; j++){		
		        gsl_ran_multinomial(r,maxBSGc,1,blockSizeProbGc,vecBSGc);
                blockSize = 0;
	            while(vecBSGc[blockSize]==0) blockSize++;
                blockSize += 1;
                nBlocks = ceil((double)vecLGc[j]/blockSize);
                gsl_ran_shuffle(r,indexGc[j],vecLGc[j],sizeof(int));
	            if (nmembers[h]>0) 
	                SPC = SPh(LUT,d,H,h,tol,thetaTilde,LGc,gammaCor,NgammaCor[h],compAlloc,nmembers[h],cetaCor,Xc,LPVcor,&Q2);
	            for (i = 0; i < vecLGc[j]; i++) 
	                vecGammaCor[j][i] = gammaCor[h][cusumVecLGc[j]+i];	 
	            for (block = 0; block < nBlocks; block++){
                    s = gsl_ran_flat(r,1.0,100000);
	                proposeBlockInd(s,vecGammaCor[j],vecLGc[j],block,blockSize,indexGc[j],cnu[j],dnu[j],vecGammaPCor[j]); 	                
	                for (i = 0; i < vecLGc[j]; i++) 
	                    gammaMP[h][cusumVecLGc[j]+i] = vecGammaPCor[j][i];	                    	            	               
	                NPJ = 0;
	                for (k = 0; k < vecLGc[j]; k++)
	                    NPJ += vecGammaCor[j][k] - vecGammaPCor[j][k];    
	                NgammaPCor = NgammaCor[h] - NPJ;                	                                        
	                if (nmembers[h]>0) 
	                    SPP = SPh(LUT,d,H,h,tol,thetaTilde,LGc,gammaMP,NgammaPCor,compAlloc,nmembers[h],cetaCor,Xc,LPVcor,&Q2);
                    if (nmembers[h]>0)
                        Acp = exp((-SPP+SPC)/(2*sigma2cor))*pow(cetaCor+1,0.5*((double)NPJ));
                    else Acp = 1.0;                                        
                    unifRV = gsl_ran_flat(r,0.0,1.0);            	                
	                if (Acp > unifRV){
                        for (k = 0; k < vecLGc[j]; k++){ 
	                        gammaCor[h][cusumVecLGc[j]+k] = gammaMP[h][cusumVecLGc[j]+k];	 
	                        vecGammaCor[j][k] = vecGammaPCor[j][k];
					    }
                        NgammaCor[h] = NgammaPCor; 
                        SPC = SPP;     
                    } 
                    else{  
	                    for (i = 0; i < vecLGc[j]; i++) 
	                        gammaMP[h][cusumVecLGc[j]+i] = gammaCor[h][cusumVecLGc[j]+i];
	                    NgammaPCor = NgammaCor[h];//IS THIS NEEDED? 
				    }				
                }  	    	   
	        }    
	    }
	    NgammaT = 0;
        for (h = 0; h < H; h++) 
            NgammaT += NgammaCor[h];
        
        // - 12 - Stick breaking Weights 
		nmb = 0;
        for (h = 0; h < H-1; h++){
            nmb += nmembers[h];
            Vdp[h] = gsl_ran_beta(r, (double) nmembers[h]+1.0, (double) d-nmb+alphaDP);
            if (h == 0)
                w[h] = Vdp[h];
            else 
                w[h] = Vdp[h] * w[h-1] * (1-Vdp[h-1]) / Vdp[h-1];
        }
        w[H-1] = w[H-2] * (1-Vdp[H-2]) / Vdp[H-2];	 
        
        // - 13 - Allocation to  clusters          
        //Rprintf("%i %s \n",sw "alloc");      
        for (h = 0; h < H; h++){                               
            if (nmembers[h]>0){
                subMeanEta = gsl_vector_subvector(meanEta2,0,NgammaCor[h]+1);
                subVarEta = gsl_matrix_submatrix(varEta2,0,0,NgammaCor[h]+1,NgammaCor[h]+1);                
                postMeanVarEtaH(LUT,d,H,h,tol,LGc,gammaCor,compAlloc,nmembers,NgammaCor,sigma2cor,cetaCor,LPVcor,Xc,thetaTilde,&subMeanEta.vector,&subVarEta.matrix,sw,0);                
                vecEta = gsl_vector_view_array(etaCor,NgammaCor[h]+1); 
                s = gsl_ran_flat(r,1.0,100000);
                sampleMN(s,NgammaCor[h]+1,&vecEta.vector,&subMeanEta.vector,&subVarEta.matrix,tol);
            }
            else if (nmembers[h] == 0) 
                for (j = 0; j < NgammaCor[h]+1; j++)
                    etaCor[j] = gsl_ran_gaussian(r,sqrt(sigma2cor*cetaCor));                                                                               
               
            if (((sw - burn) >= 0) && (((sw - burn ) % thin) == 0)){
                move = 0;
                fprintf(out_file18, "%f ", etaCor[move++]);
                for (j = 0; j < LGc; j++)
                    if (gammaCor[h][j]==1) fprintf(out_file18, "%f ", etaCor[move++]); else fprintf(out_file18, "%f ", 0.0);                
                if (h == (H-1)) fprintf (out_file18, "\n");                
            }
                
            vecEta = gsl_vector_view_array(etaCor,NgammaCor[h]+1);                                                            
            compAlloc2[0] = h;                 
            setBaseZBSgh(LUT,d,H,h,LGc,gammaCor,compAlloc2,Xc,vizZstar);                                                               	            
            Zstar = gsl_matrix_view_array(vizZstar,LUT,NgammaCor[h]+1);	                                                                    
            vecMean = gsl_vector_view_array(Zeta,LUT);
            gsl_blas_dgemv(CblasNoTrans,1.0,&Zstar.matrix,&vecEta.vector,0.0,&vecMean.vector);
            for (i = 0; i < d; i++){
                pdh[h][i] = w[h];
                for (t = 0; t < LUT; t++)
                    pdh[h][i] *= gsl_ran_gaussian_pdf(theta[i+t*d]-Zeta[t],sqrt(sigma2cor*exp(LPVcor[t])));//sqrt(sigma2t[t]) change                                                                       
            } 
        }
       
        for (i = 0; i < d; i++){
				temp = 0.0;
                for (h = 0; h < H; h++)
				    temp += pdh[h][i];
				for (h = 0; h < H; h++)
				    pdh[h][i] /= temp;    
		}		                 
	    
        s = gsl_ran_flat(r,1.0,100000);
        allocation(s,d,H,pdh,compAlloc,0);
        findNmembers(d,H,compAlloc,nmembers);
        
        // - 14 - DP concentration parameter
        s = gsl_ran_flat(r,1.0,100000); 
        newAlphaDP = updateAlpha(s,d,H,Alphaa,Alphab,TruncAlpha,nmembers,alphaDP);
        alphaDP = newAlphaDP;
       
        // - 15 - Update delta and omega
        //Rprintf("%i %s \n",sw,"delta & omega");        
        for (h = 0; h < H; h++){
            if (nmembers[h]>0){
                subMeanEta = gsl_vector_subvector(meanEta2,0,NgammaCor[h]+1);
                subVarEta = gsl_matrix_submatrix(varEta2,0,0,NgammaCor[h]+1,NgammaCor[h]+1);            
                postMeanVarEtaH(LUT,d,H,h,tol,LGc,gammaCor,compAlloc,nmembers,NgammaCor,sigma2cor,cetaCor,LPVcor,
                                Xc,thetaTilde,&subMeanEta.vector,&subVarEta.matrix,sw,0);                            
                for (i = 0; i < LGc; i++) 
	                gammaC[i] = gammaCor[h][i];                 
                cSqResh(LUT,d,h,compAlloc,gammaC,NgammaCor,LGc,Xc,&subMeanEta.vector,theta,sqResCCor);                
                
                if ((sw - burn) >= 0 && 1==0){
                    move = 0; 
                    fprintf(out_file24, "%f ", gsl_vector_get(&subMeanEta.vector,move++));
                    for (j = 0; j < LGc; j++)
                        if (gammaCor[h][j]==1) fprintf(out_file24, "%f ", gsl_vector_get(&subMeanEta.vector,move++)); else fprintf(out_file24, "%f ", 0.0);                
                    fprintf (out_file24, "\n");	       
	                for (j = 0; j < (LUT*d); j++)
                        fprintf(out_file24, "%f ", theta[j]);
                    fprintf (out_file24, "\n");
                    for (j = 0; j < (LUT*d); j++)
                        fprintf(out_file24, "%f ", sqResC[j]);
                    fprintf (out_file24, "\n");
                }
                                   
			}
        }
        SPC = 0;
        for (h = 0; h < H; h++){        
            if (nmembers[h]>0) 
                SPC += SPh(LUT,d,H,h,tol,thetaTilde,LGc,gammaCor,NgammaCor[h],compAlloc,nmembers[h],cetaCor,Xc,LPVcor,&Q2);	            
            if (1==0) fprintf(out_file24,"%i %i %s %i %i %f\n",sw,h,"SPC a & d:",nmembers[0],nmembers[1],SPC);            
        }
        QFC = 0.0;
        for (i = 0; i < LDc; i++)
            QFC += pow(omega[i],2);                 
        for (j = 0; j < NDc; j++){
		    if ((sw % batchL)==0 && WB > 0){ 
	            if (acceptOmega[j] > 0.25 && tuneOmega[j] < omegaUL) tuneOmega[j] += MIN(0.01,1/sqrt(WB)) * tuneOmega[j]; 
	            if (acceptOmega[j] <= 0.25 && tuneOmega[j] > omegaLL) tuneOmega[j] -= MIN(0.01,1/sqrt(WB)) * tuneOmega[j];
	            acceptOmega[j] = 0.0;   
	            if (tuneOmega[j] < omegaLL) tuneOmega[j] = omegaLL;
	            if (tuneOmega[j] > omegaUL) tuneOmega[j] = omegaUL;
	            fprintf(out_file25,"%f ",tuneOmega[j]);
            }            
		    gsl_ran_multinomial(r,maxBSDc,1,blockSizeProbDc,vecBSDc);
            blockSize = 0;
            while(vecBSDc[blockSize]==0) blockSize++;
            blockSize += 1;
            nBlocks = ceil((double)vecLDc[j]/blockSize);
            gsl_ran_shuffle(r,indexDc[j],vecLDc[j],sizeof(int)); 
            for (block = 0; block < nBlocks; block++){  
                s = gsl_ran_flat(r,1.0,100000);
                proposeBlockInd(s,vecDeltaCor[j],vecLDc[j],block,blockSize,indexDc[j],cfi[j],dfi[j],vecDeltaPCor[j]);            
                NPJ = 0;  
                for (k = 0; k < vecLDc[j]; k++)
                    NPJ += vecDeltaPCor[j][k];
                NdeltaPCor = NdeltaCor - vecNdeltaCor[j] + NPJ;
                for (k = 0; k < vecLDc[j]; k++)
	                deltaPCor[cusumVecLDc[j]+k] = vecDeltaPCor[j][k]; 	       	            	           	             	                            
                if (NPJ > 0){                
                    subAlphaHat = gsl_vector_subvector(alphaHat,0,NPJ);
                    subD = gsl_matrix_submatrix(D,0,0,NPJ,NPJ);
                    DeltaAlphaHat(LUT,d,tol,LPVcor,sqResCCor,deltaPCor,NPJ,cusumVecLDc[j],cusumVecLDc[j+1],
                                  Zc,sigma2cor,sigma2t,comega,&subD.matrix,&subAlphaHat.vector);             
                    subAlphaP = gsl_vector_subvector(alphaP,0,NPJ);
                    gsl_matrix_scale(&subD.matrix,tuneOmega[j]);
                    s = gsl_ran_flat(r,1.0,100000);
                    sampleMN(s,NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);
                    logMVNormP = logMVNormalpdf(NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol); 		                                                                       
                }else logMVNormP = 0.0;                                             
                move = 0;
                for (k = 0; k < vecLDc[j]; k++){
                    if (vecDeltaPCor[j][k]==1) omegaPD[cusumVecLDc[j]+k] = gsl_vector_get(&subAlphaP.vector,move++);                
                    else omegaPD[cusumVecLDc[j]+k] = 0;
                } 
                for (i = 0; i < LUT; i++){
                    LPVPcor[i] = LPVcor[i];
                    for (k = cusumVecLDc[j]; k < cusumVecLDc[j+1]; k++)
                        LPVPcor[i] += (omegaPD[k]-omega[k])*Zc[LUT*d+k*LUT*d+i*d]; 
                    sigma2tP[i] = sigma2cor*exp(LPVPcor[i]);
                } 
                for (t = 0; t < LUT; t++) 
                    for (i = 0; i < d; i++)
                        thetaTildeP[i+t*d] = exp(-LPVPcor[t]/2)*theta[i+t*d];  
                QFP = QFC;
                for (k = cusumVecLDc[j]; k < cusumVecLDc[j+1]; k++)
                    QFP += pow(omegaPD[k],2)-pow(omega[k],2);
                detR = 0.0;
                for (i = 0; i < LUT; i++)
                    detR += LPVcor[i] - LPVPcor[i];
                detR *= 0.5*d;
	            SPP = 0;
                for (h = 0; h < H; h++){        
                    if (nmembers[h]>0) 
                        SPP += SPh(LUT,d,H,h,tol,thetaTildeP,LGc,gammaCor,NgammaCor[h],compAlloc,nmembers[h],cetaCor,Xc,LPVPcor,&Q2);
                    if (1==0) fprintf(out_file24,"%i %i %s %i %i %f\n",sw,h,"SPP a & d:",nmembers[0],nmembers[1],SPP);
                }     	                 
	            //probability of reverse direction
	            for (h = 0; h < H; h++){
                    if (nmembers[h]>0){
                        subMeanEta = gsl_vector_subvector(meanEta2,0,NgammaCor[h]+1);
                        subVarEta = gsl_matrix_submatrix(varEta2,0,0,NgammaCor[h]+1,NgammaCor[h]+1);            
                        postMeanVarEtaH(LUT,d,H,h,tol,LGc,gammaCor,compAlloc,nmembers,NgammaCor,sigma2cor,cetaCor,LPVPcor,
                                        Xc,thetaTildeP,&subMeanEta.vector,&subVarEta.matrix,sw,0);                            
                        for (i = 0; i < LGc; i++) 
	                        gammaC[i] = gammaCor[h][i];                 
                        cSqResh(LUT,d,h,compAlloc,gammaC,NgammaCor,LGc,Xc,&subMeanEta.vector,theta,sqResPCor);                                   
			        }
                }	    	                                                                                             
                if (vecNdeltaCor[j] > 0){                
	                subAlphaHat = gsl_vector_subvector(alphaHat,0,vecNdeltaCor[j]);
                    subD = gsl_matrix_submatrix(D,0,0,vecNdeltaCor[j],vecNdeltaCor[j]);
                    DeltaAlphaHat(LUT,d,tol,LPVPcor,sqResPCor,deltaCor,vecNdeltaCor[j],cusumVecLDc[j],cusumVecLDc[j+1],
                                  Zc,sigma2cor,sigma2tP,comega,&subD.matrix,&subAlphaHat.vector);
                    gsl_matrix_scale(&subD.matrix,tuneOmega[j]);
                    move=0;
                    for (k = 0; k < vecLDc[j]; k++)
                        if (vecDeltaCor[j][k]==1) BaseSubAlpha[move++] = omega[cusumVecLDc[j]+k];
                    subAlphaP = gsl_vector_view_array(BaseSubAlpha,vecNdeltaCor[j]);
                    logMVNormC = logMVNormalpdf(vecNdeltaCor[j],&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);                                 
	            }else logMVNormC = 0.0;   
                Acp = exp(detR+(-SPP+SPC)/(2*sigma2cor)+logMVNormC-logMVNormP+(QFC-QFP)/(2*comega))*
                      pow(2*M_PI*comega,0.5*(NdeltaCor-NdeltaPCor));
                unifRV = gsl_ran_flat(r,0.0,1.0);    
                //fprintf(out_file24, "%s %i %f %f %f %f %f %f %f %f %f %f %i %i \n","d-o: ", sw, Acp,
                //detR,SPP,SPC,sigma2cor,logMVNormC,logMVNormP,QFC,QFP,comega,NdeltaCor,NdeltaPCor);       
	            if (Acp > unifRV){	      	   	    
                    if (NPJ > 0) acceptOmega[j] += 1/((double)batchL);
                    for (k = 0; k < vecLDc[j]; k++){
                        deltaCor[cusumVecLDc[j]+k] = deltaPCor[cusumVecLDc[j]+k];
                        vecDeltaCor[j][k] = vecDeltaPCor[j][k];
                        omega[cusumVecLDc[j]+k] = omegaPD[cusumVecLDc[j]+k];                  
                    }                
                    for (i = 0; i < LUT; i++){
                         LPVcor[i] = LPVPcor[i];
                         sigma2t[i] = sigma2tP[i];
				    }
				    for (i = 0; i < (LUT*d); i++){
                        thetaTilde[i] = thetaTildeP[i];
                        sqResCCor[i] = sqResPCor[i];
					}
                    SPC = SPP;
                    QFC = QFP;
                    NdeltaCor =  NdeltaPCor;
                    vecNdeltaCor[j] = NPJ;	      	     
                }
                else{
					for (k = 0; k < vecLDc[j]; k++){
                        deltaPCor[cusumVecLDc[j]+k] = deltaCor[cusumVecLDc[j]+k];
					    omegaPD[cusumVecLDc[j]+k] = omega[cusumVecLDc[j]+k];
					}
				}                      
            }	    
	    }
        
        // - 13 - sigma2cor
        //Rprintf("%i %s \n",sw,"sigma2cor");
		if (HNsigcor[0]==0){ 
            alphaG = sigmaCorParams[0]; 
            betaG = sigmaCorParams[1];
			sigma2cor = 1/gsl_ran_gamma(r,alphaG+0.5*LUT*d,1/(betaG+0.5*SPC));
			for (i = 0; i < LUT; i++)
                sigma2t[i] = sigma2cor*exp(LPVcor[i]);
	    }else if (HNsigcor[0]==1){
            if ((sw % batchL)==0 && WB > 0){ 
	            if (acceptSigma2 > 0.25 && tuneSigma2R[0] < s2UL) tuneSigma2R[0] += MIN(0.01,1/sqrt(WB)) * tuneSigma2R[0]; 
	            if (acceptSigma2 <= 0.25 && tuneSigma2R[0] > s2LL) tuneSigma2R[0] -= MIN(0.01,1/sqrt(WB)) * tuneSigma2R[0];
	            acceptSigma2 = 0.0;   
	            if (tuneSigma2R[0] < s2LL) tuneSigma2R[0] = s2LL;
	            if (tuneSigma2R[0] > s2UL) tuneSigma2R[0] = s2UL;
	            fprintf(out_file25,"%f ",tuneSigma2R[0]);
            }            
            phi2G = sigmaCorParams[0];
            sigma2P = sigma2cor + gsl_ran_gaussian(r,sqrt(tuneSigma2R[0]));
            while (sigma2P <= 0) sigma2P = sigma2cor + gsl_ran_gaussian(r,sqrt(tuneSigma2R[0]));            
            Acp = exp(-0.5*LUT*d*log(sigma2P/sigma2cor) + 0.5*SPC*(1/sigma2cor-1/sigma2P) + 
                  0.5*(sigma2cor-sigma2P)/phi2G);
	        unifRV = gsl_ran_flat(r,0.0,1.0);
	        //fprintf(out_file24, "%s %i %f ","sigma2cor: ", sw, Acp);
            if (Acp > unifRV){
                sigma2cor = sigma2P;
	            acceptSigma2 += 1/((double)batchL);  
	            for (i = 0; i < LUT; i++)
                    sigma2t[i] = sigma2cor*exp(LPVcor[i]);
            }
	    }
        
        // - 14 - c_eta 
	    //Rprintf("%i %s \n",sw,"c_eta");
        if ((sw % batchL)==0 && WB > 0){ 
	        if (acceptCetaCor > 0.25 && tuneCbCor[0] < GULcor) tuneCbCor[0] += MIN(0.01,1/sqrt(WB)) * tuneCbCor[0]; 
	        if (acceptCetaCor <= 0.25 && tuneCbCor[0] > GLLcor) tuneCbCor[0] -= MIN(0.01,1/sqrt(WB)) * tuneCbCor[0];
	        acceptCetaCor = 0.0;   
	        if (tuneCbCor[0] < GLLcor) tuneCbCor[0] = GLLcor;
	        if (tuneCbCor[0] > GULcor) tuneCbCor[0] = GULcor;
	        fprintf(out_file25,"%f ",tuneCbCor[0]);
        }     
        cetahat = 1;//starting value for NR
        SPC = 0;
        Q2 = 0; 
        for (h = 0; h < H; h++){          
            if (nmembers[h]>0){				   	            
                SPC += SPh(LUT,d,H,h,tol,thetaTilde,LGc,gammaCor,NgammaCor[h],compAlloc,nmembers[h],cetaCor,Xc,LPVcor,&temp);
                Q2 += temp;
            }
	    }
        elPrime = 99.9;
        while(elPrime > 0.000000001 || -elPrime > 0.000000001){
            Sprime = -Q2/pow(cetahat+1,2);
            SDprime = 2*Q2/pow(cetahat+1,3);
            elPrime = -0.5*(NgammaT+H)/(cetahat+1)-Sprime/(2*sigma2cor)-(cetaCorParams[0]+1)/cetahat+cetaCorParams[1]/pow(cetahat,2);
            elDPrime = 0.5*(NgammaT+H)/(pow(cetahat+1,2))-SDprime/(2*sigma2cor)+(cetaCorParams[0]+1)/pow(cetahat,2)-2*cetaCorParams[1]/pow(cetahat,3);
            cetahat -= elPrime / elDPrime;
        }
        cetaP = cetahat + gsl_ran_gaussian(r,sqrt(-tuneCbCor[0]/elDPrime));
	    while(cetaP < 0) cetaP = cetahat + gsl_ran_gaussian(r,sqrt(-tuneCbCor[0]/elDPrime));
	    SPP = 0;
        for (h = 0; h < H; h++){          
            if (nmembers[h]>0){
                SPP += SPh(LUT,d,H,h,tol,thetaTilde,LGc,gammaCor,NgammaCor[h],compAlloc,nmembers[h],cetaP,Xc,LPVcor,&temp);
            }
	    }	    
	    Acp = exp(-0.5*(NgammaT+H)*(log(cetaP+1)-log(cetaCor+1))+(-SPP+SPC)/(2*sigma2cor)-(cetaCorParams[0]+1)*(log(cetaP)-log(cetaCor))
              + cetaCorParams[1]*(1/cetaCor-1/cetaP))*
                gsl_ran_gaussian_pdf(cetaCor-cetahat,sqrt(-tuneCbCor[0]/elDPrime))/
                gsl_ran_gaussian_pdf(cetaP-cetahat,sqrt(-tuneCbCor[0]/elDPrime));
	    unifRV = gsl_ran_flat(r,0.0,1.0);
	    //fprintf(out_file24, "%s %i %f ","ceta cor: ", sw, Acp);
        if (Acp > unifRV){
            cetaCor = cetaP;
	        acceptCetaCor += 1/((double)batchL);  
	    }
        
	    // - 15 - c_omega
        //Rprintf("%i %s \n",sw,"c_omega");
        if (HNco[0]==0){ 
			alphaG = comegaParams[0]; betaG = comegaParams[1];
		    comega = 1/gsl_ran_gamma(r,alphaG+0.5*NdeltaCor,1/(betaG+0.5*QFC));		    
		    //fprintf(out_file24, "%s %i %f %f %i %f %f \n","comega: ", sw, alphaG, betaG, NdeltaCor, QFC, comega);
		    if ((sw % batchL)==0 && WB > 0) fprintf(out_file25,"\n");
		}else if (HNco[0]==1){
			phi2G=comegaParams[0];
            if ((sw % batchL)==0 && WB > 0){ 
	            if (acceptComega > 0.25 && tuneComega[0] < coUL) tuneComega[0] += MIN(0.01,1/sqrt(WB)) * tuneComega[0]; 
	            if (acceptComega <= 0.25 && tuneComega[0] > coLL) tuneComega[0] -= MIN(0.01,1/sqrt(WB)) * tuneComega[0];
	            acceptComega = 0.0;   
	            if (tuneComega[0] < coLL) tuneComega[0] = coLL;
	            if (tuneComega[0] > coUL) tuneComega[0] = coUL;
	            fprintf(out_file25,"%f \n",tuneComega[0]);
            }            
            omegaP = comega + gsl_ran_gaussian(r,sqrt(tuneComega[0]));
            while (omegaP <= 0) omegaP = comega + gsl_ran_gaussian(r,sqrt(tuneComega[0]));            
            Acp = exp(-0.5*NdeltaCor*log(omegaP/comega) + (QFC/2)*(1/comega-1/omegaP) + 
                  (comega-omegaP)/(2*phi2G));	        
	        unifRV = gsl_ran_flat(r,0.0,1.0);
	        //fprintf(out_file24, "%s %i %f %f \n","comega: ", sw, Acp, comega);
            if (Acp > unifRV){
                comega = omegaP;
	            acceptComega += 1/((double)batchL);                  
            }
	    }
        
	    //Post burn-in 	    
        if (((sw - burn) >= 0) && (((sw - burn ) % thin) == 0)){

            // - 16 - etaCor           
            //subMeanEta = gsl_vector_subvector(meanEta2,0,NgammaCor+1);//not needed
            //subVarEta = gsl_matrix_submatrix(varEta2,0,0,NgammaCor+1,NgammaCor+1);//when all code is running
            //postMeanVarEta(LUT,d,tol,gammaCor,NgammaCor,LGc,sigma2cor,cetaCor,LPVcor,Xc,thetaTilde,&subMeanEta.vector,&subVarEta.matrix,sw);
            //gsl_matrix_scale(&subVarEta.matrix,sigma2cor);
            //vecEta = gsl_vector_view_array(etaCor,NgammaCor+1); 
            //s = gsl_ran_flat(r,1.0,100000);
            //sampleMN(s,NgammaCor+1,&vecEta.vector,&subMeanEta.vector,&subVarEta.matrix,tol);
            
            // - 17 - deviance = -2*LogLikelihood
            SPC = ScalcMultLong(m,p,tol,LG,Ngamma,niMax,niVec,cusumniVec,N,sigma2ij,dimLPC,LPC,Y,X,gamma,RiAll,intime,ceta,&Q2);
            
            detR = 0.0;
            
            for (i = 0; i < N; i++)
                for (j = 0; j < p; j++)
                    detR += LPV[i][j];
           
            for (j = 0; j < p; j++)
                detR += N*log(sigma2zk[j]);
            
            //FIX
            //for (j = 0; j < LUT; j++){
            //    //specify RtC
            //    detS = det(p,RtC);
            //    detR += FUT[j]*log(detS); //where FUT is the vector of frequencies of unique times 
			//}
            
            dev0 = SPC + (Ngamma+p)*log(ceta+1) + detR + N*p*log(2*M_PI);
            dev[0] += dev0;

            
            SPP = 0; // FIX NormalQuadr(p,m,LG,Ngamma,Ytilde,sigma2ij,X,gamma,RtCinv,eta);
            dev1 = SPP + detR + N*p*log(2*M_PI);
            dev[1] += dev1;

	        //Rprintf("%s %f %f %f %f %f \n","deviance:",dev[0],m*p*log(2*M_PI),detR,
	        //(Ngamma+p)*log(ceta+1),SPC);
	        
	        // - 18 - Write to files
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
            
            for (j = 0; j < p; j++)
                fprintf(out_file5, "%f ", sigma2zk[j]);
            fprintf(out_file5, "\n");
            
            for (j = 0; j < p; j++)
                fprintf(out_file6, "%f ", calpha[j]);
            fprintf(out_file6, "\n");
            
            move = 0;
            for (k = 0; k < p; k++){
                for (j = 0; j < (LG+1); j++)
                    if ((j==0) || (j>0 && gamma[k][j-1]==1)) fprintf(out_file7, "%f ", eta[move++]);
                    else fprintf(out_file7, "%f ", 0.0);
            }
            fprintf (out_file7, "\n");
            
            for (k = 0; k < p; k++)
				for (l = 0; l < p; l++)
                    for (j = 0; j < LK; j++)
                        fprintf(out_file8, "%f ", psi[k][l][j]);
            fprintf (out_file8, "\n");
            
            for (k = 0; k < p; k++)
				for (l = 0; l < p; l++)
                    for (j = 0; j < LK; j++)
                        fprintf(out_file9, "%i ", ksi[k][l][j]);
            fprintf (out_file9, "\n");
            
            for (j = 0; j < (p*p); j++)
                fprintf(out_file10, "%f ", cPsi[j]);
            fprintf(out_file10, "\n");
            
            for (t = 0; t < LUT; t++)
                for (j = 0; j < p*p; j++)
                    fprintf(out_file11, "%f ", Rt[t][j]);
            fprintf(out_file11, "\n");
            
            for (h = 0; h < H; h++){
                for (j = 0; j < LGc; j++)
                    fprintf(out_file12, "%i ", gammaCor[h][j]);
                if (h == (H-1)) fprintf(out_file12, "\n");
            }
            
            for (j = 0; j < LDc; j++)
                fprintf(out_file13, "%i ", deltaCor[j]);
            fprintf(out_file13, "\n");

            for (j = 0; j < LDc; j++) 
                fprintf(out_file14, "%f ", omega[j]);
            fprintf (out_file14, "\n");

            fprintf(out_file15, "%0.10f \n", sigma2cor);            

            fprintf(out_file16, "%f \n", cetaCor);

            fprintf(out_file17, "%f \n", comega);                       
            
            fprintf(out_file19, "%f %f \n", dev0, dev1);

            if (sw==contParams[1]){                    
                for (t = 0; t < LUT; t++)
                    for (j = 0; j < p*p; j++)
                        fprintf(out_file20, "%f ", Dt[t][j]);
                fprintf(out_file20, "\n");
                                    
                for (t = 0; t < LUT; t++)
                    for (j = 0; j < p*p; j++)
                        fprintf(out_file20, "%f ", Et[t][j]);
                fprintf(out_file20, "\n");                
			}
			
			for (i = 0; i < d; i++)
                fprintf(out_file21, "%i ", compAlloc[i]);
            fprintf (out_file21, "\n");
            
            for (h = 0; h < H; h++)
                fprintf(out_file22, "%i ", nmembers[h]);
            fprintf(out_file22, "\n");
			
			fprintf(out_file23, "%f \n", alphaDP);
			
			for (i = 0; i < d; i++){
                for (h = 0; h < H; h++)
                    fprintf(out_file24,"%f ",pdh[h][i]);                                    
                fprintf(out_file24, "\n");
			}      
        }
        // If sw needs to be printed  
        if ((sw==(sweeps-1)) && (!((sweeps % 1000)==0))) Rprintf("%i %s \n",sweeps, "posterior samples...");
    }//end of sw
   
    //Update LASTWB
    contParams[2] = WB;

    //Free up random number generator
    gsl_rng_free (r);

    //Close files 
    fclose(out_file1); fclose(out_file2); fclose(out_file3);
    fclose(out_file4); fclose(out_file5); fclose(out_file6);
    fclose(out_file7); fclose(out_file8); fclose(out_file9);
    fclose(out_file10); fclose(out_file11); fclose(out_file12);    
    fclose(out_file13); fclose(out_file14); fclose(out_file15);
    fclose(out_file16); fclose(out_file17); fclose(out_file18);
    fclose(out_file19); fclose(out_file20); fclose(out_file21); 
    fclose(out_file22); fclose(out_file23); fclose(out_file24);fclose(out_file25); 

    //Free up gsl matrices
    gsl_matrix_free(St); gsl_matrix_free(PSIt); gsl_matrix_free(CopyPSIt);   
    gsl_matrix_free(gEt); gsl_matrix_free(gDt); gsl_matrix_free(gRt);
    gsl_matrix_free(D); gsl_matrix_free(varEta);
    gsl_vector_free(meanEta);  gsl_vector_free(alphaHat);
    gsl_vector_free(alphaP); 
    gsl_matrix_free(L); 
    gsl_matrix_free(RiAll);
    gsl_matrix_free(varPsi); gsl_vector_free(meanPsi);
    gsl_matrix_free(A);
    gsl_matrix_free(varEta2);
    gsl_vector_free(mutheta);
    gsl_vector_free(meanEta2);
    gsl_vector_free(vecZero1);

    //
    free(BaseXg);

    //Free jagged vectors
	for (j = 0; j < NG; j++) {free(indexG[j]); free(vecGamma[j]); free(vecGammaP[j]);}
	for (j = 0; j < ND; j++) {free(indexD[j]); free(vecDelta[j]); free(vecDeltaP[j]);}
	for (j = 0; j < NK; j++) {free(indexK[j]); free(vecKsi[j]); free(vecKsiP[j]);}
    for (j = 0; j < NGc; j++) {free(indexGc[j]); free(vecGammaCor[j]); free(vecGammaPCor[j]);}
	for (j = 0; j < NDc; j++) {free(indexDc[j]); free(vecDeltaCor[j]); free(vecDeltaPCor[j]);}	
}
