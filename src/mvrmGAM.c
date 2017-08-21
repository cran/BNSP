/* mvrmGAM.c Bayesian semiparametric regression models for mean and variance functions 
 * Copyright (C) 2017  Georgios Papageorgiou, gpapageo@gmail.com
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
//#include "matalg.h"
//#include "pdfs.h"
//#include "sampling.h"
//#include "other.functions.h"
//#include "mathm.h"
#include "spec.BCM.h" 
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

void mvrmC(int *seed1, char **WorkingDir, int *WF1,
         int *sweeps1, int *burn1, int *thin1,
         double *Y, double *X, double *Z, int *n1, int *LG1, int *LD1,
         double *blockSizeProb1, int *maxBSG1, double *blockSizeProb2, int *maxBSD1, 
         double *f1, double *g1, double *h,
         int *NG1, int *ND1, int *vecLG, int *vecLD, int *cusumVecLG, int *cusumVecLD, int *MVLD,
         double *cetaParams, double *calphaParams, double *pimu, double *pisigma,
         int *HN1, double *sigmaParams)
{
    gsl_set_error_handler_off();

    // Random number generator initialization
    int seed = seed1[0]; 
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);
    unsigned long int s; 

    // Specify directory
    int WF = WF1[0]; // indicator: 1 = write files, 0 = no files.
    const int MAX_PATH = 300;
    char path_name[MAX_PATH + 1];

    // Open files
    FILE *out_file1, *out_file2, *out_file3, *out_file4, *out_file5, *out_file6, *out_file7;

    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.gamma.txt");
    out_file1 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.delta.txt");
    out_file2 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.alpha.txt");
    out_file3 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.sigma2.txt");
    out_file4 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.cbeta.txt");
    out_file5 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.calpha.txt");
    out_file6 = fopen(path_name, "wt");
    snprintf(path_name, MAX_PATH, "%s%s", *WorkingDir, "BNSP.beta.txt");
    out_file7 = fopen(path_name, "wt");
    
    // Sweeps, burn-in period, thin
    int sweeps = sweeps1[0]; 
    int burn = burn1[0]; 
    int thin = thin1[0]; 

    // Dimensions
    int n = n1[0]; 
    int LG = LG1[0];
    int LD = LD1[0];
    int NG = NG1[0];
    int ND = ND1[0];

    //Tolerance level 
    double tol = 0.01; 
    
    //Declare variables for loops:
    int sw, i, j, k;

    // Prior parameters
      //c.eta
    double alphaeta = cetaParams[0];
    double betaeta = cetaParams[1];
      //c.alpha
    double alphaalpha = calphaParams[0];
    double betaalpha = calphaParams[1];
      //pi_mu & pi_sigma
    double cmu[NG];
    double dmu[NG];
    double csigma[ND];
    double dsigma[ND];     
    for (k = 0; k < NG; k++){
        cmu[k] = pimu[k];
        dmu[k] = pimu[k+NG];
	}    
    for (k = 0; k < ND; k++){
        csigma[k] = pisigma[k];
        dsigma[k] = pisigma[k+ND];
	}    
      //sigma2
    int HN=HN1[0];
    double phi2sigma;  //for half normal prior
    double alphasigma; //for IG prior
    double betasigma;
    if (HN==1) phi2sigma=sigmaParams[0];
    else if (HN==0) {alphasigma=sigmaParams[0]; betasigma=sigmaParams[1];} 

    //Tuning papameters
    double f = f1[0];
    double g = g1[0];

    // Declare quantities of interest
    double eta[LG+1];
    int gamma[LG];
    int Ngamma;
    int *vecGamma[NG];
    for (j = 0; j < NG; j++)
        vecGamma[j] = malloc(sizeof(int) * vecLG[j]);
    int vecNgamma[NG];
    int delta[LD];
    int Ndelta;
    int *vecDelta[ND];
    for (j = 0; j < ND; j++)
        vecDelta[j] = malloc(sizeof(int) * vecLD[j]);
    int vecNdelta[ND];
    double sigma2z[n];
    double alpha[LD];
    double sigma2;
    double calpha;
    double ceta;

    // Declare other quantities
    int move;
    double yTilde[n];
    double yTildeP[n];
    double sigma2zP[n];  
    double LPV[n];
    double LPVP[n]; 
    double alphaPD[LD];
    double BaseSubAlpha[MVLD[0]];
    double sqResC[n];
    double sqResP[n];
    double Acp, unifRV, QFC, QFP, detR,
           logMVNormC, logMVNormP, SPC, SPP, sigma2P;
    double cetahat, Sprime, SDprime, elPrime, elDPrime, Q2, cetaP;
    
    //For selecting block size gamma_B or delta_B
    int block, blockSize;
    double nBlocks;
    int maxBSG = maxBSG1[0];
    double blockSizeProbG[maxBSG];
    for (k = 0; k < maxBSG; k++)
        blockSizeProbG[k] = blockSizeProb1[k];
    unsigned int vecBSG[maxBSG];
    int maxBSD = maxBSD1[0]; 
    double blockSizeProbD[maxBSD];
    for (k = 0; k < maxBSD; k++)
        blockSizeProbD[k] = blockSizeProb2[k];
    unsigned int vecBSD[maxBSD];

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

    //Proposed gamma & delta
    int gammaP[LG];
    int NgammaP;
    int *vecGammaP[NG];
    for (j = 0; j < NG; j++)
        vecGammaP[j] = malloc(sizeof(int) * vecLG[j]);        
    int deltaP[LD];
    int NdeltaP;
    int *vecDeltaP[ND];
    for (j = 0; j < ND; j++)
        vecDeltaP[j] = malloc(sizeof(int) * vecLD[j]);    
    int NPJ; //nonzero in jth proposed

    // Declare and allocate gsl vectors and matrices
    gsl_matrix *D = gsl_matrix_alloc(MVLD[0],MVLD[0]);
    gsl_matrix *varEta = gsl_matrix_alloc(LG+1,LG+1);
    gsl_vector *meanEta = gsl_vector_alloc(LG+1);
    gsl_vector *alphaHat = gsl_vector_alloc(MVLD[0]);
    gsl_vector *alphaP = gsl_vector_alloc(MVLD[0]);

    // GSL Matrix and vector views
    gsl_matrix_view subD, subVarEta;
    gsl_vector_view subAlphaHat,subAlphaP, subMeanEta, vecEta;
    
    //Make adaptive
    int batchL = 50; //batch length
    double WB; //which batch 
    double acceptCeta = 0.0;
    double GLL = 0.01;
    double GUL = 200;
    double acceptSigma2 = 0.0;
    double FLL = 0.01;
    double FUL = 200;
    double acceptAlpha[ND];
	for (j = 0; j < ND; j++)
	    acceptAlpha[j] = 0.0;	
    double HLL = 3;
    double HUL = 200;

    // Sampler initialization
    // - 1 - Vec of gammas
    for (j = 0; j < NG; j++)
        for (k = 0; k < vecLG[j]; k++)
            vecGamma[j][k] = 1; //gsl_ran_bernoulli(r,cmu/(cmu+dmu)); 
    for (j = 0; j < NG; j++){
	    vecNgamma[j] = 0;
        for (k = 0; k < vecLG[j]; k++)
            vecNgamma[j] += vecGamma[j][k];
    }
	for (j = 0; j < NG; j++)
	    for (k = 0; k < vecLG[j]; k++)
	        gamma[cusumVecLG[j]+k] = vecGamma[j][k];
	Ngamma = 0;
	for (j = 0; j < NG; j++)
	    Ngamma += vecNgamma[j]; 
    for (k = 0; k < LG; k++) 
        gammaP[k] = gamma[k];
    NgammaP = Ngamma;
        
    // - 2 - Vec of delta
    for (j = 0; j < ND; j++)
        for (k = 0; k < vecLD[j]; k++)
            vecDelta[j][k] = 0; //gsl_ran_bernoulli(r,csigma[j]/(csigma[j]+dsigma[j]));
    for (j = 0; j < ND; j++){
	    vecNdelta[j] = 0;
        for (k = 0; k < vecLD[j]; k++)
            vecNdelta[j] += vecDelta[j][k];
    }
    for (j = 0; j < ND; j++)
	    for (k = 0; k < vecLD[j]; k++)
	        delta[cusumVecLD[j]+k] = vecDelta[j][k];        
    Ndelta = 0;
    for (j = 0; j < ND; j++)
        Ndelta += vecNdelta[j];
    for (k = 0; k < LD; k++) 
        deltaP[k] = delta[k];
        
    // - 3 - c_alpha
    calpha = 110.0;//1/gsl_ran_gamma(r,alphaalpha,1/betaalpha);
    
    // - 4 - sigma2, alpha, LPV, sigma2t, QFC
    sigma2 = 1.0;
    for (k = 0; k < LD; k++){
        if (delta[k]==0) alpha[k] = 0.0;
        else alpha[k] = gsl_ran_gaussian(r,sqrt(calpha)); 
    }
    for (i = 0; i < n; i++){
        LPV[i] = 0.0;
        for (k = 0; k < LD; k++)
            if (delta[k]==1) LPV[i] += alpha[k]*Z[n+k*n+i];
        sigma2z[i] = sigma2*exp(LPV[i]);
        yTilde[i] = exp(-LPV[i]/2)*Y[i];
    }                
    QFC = 0.0;
    for (k = 0; k < LD; k++)
        QFC += pow(alpha[k],2);
    for (k = 0; k < LD; k++)
        alphaPD[k] = alpha[k];
    
    // - 5 - c_eta
    ceta = 1;
    cetahat = ceta; //starting value is the current value
    SPC = SPcalc(n,1,tol,yTilde,gamma,Ngamma,LG,ceta,X,LPV,&Q2);
    elPrime = 99.9;
    while(elPrime > 0.000000001 || -elPrime > 0.000000001){
        Sprime = -Q2/pow(cetahat+1,2);
        SDprime = 2*Q2/pow(cetahat+1,3);
        elPrime = -0.5*(Ngamma+1)/(cetahat+1)-Sprime/(2*sigma2)-(alphaeta+1)/cetahat+betaeta/pow(cetahat,2);
        elDPrime = 0.5*(Ngamma+1)/(pow(cetahat+1,2))-SDprime/(2*sigma2)+(alphaeta+1)/pow(cetahat,2)-2*betaeta/pow(cetahat,3);
        cetahat -= elPrime / elDPrime;
    }
    ceta = cetahat;

    //#############################################SAMPLER
    for (sw = 0; sw < sweeps; sw++){
        if (sw==0) Rprintf("%i %s \n",sw+1, "posterior sample...");
        if (((sw+1) % 1000)==0) Rprintf("%i %s \n",sw+1, "posterior samples...");
        
	    modf(sw/batchL,&WB);
       	    
        // - 1 - Update gamma 
        //Rprintf("%i %s \n",sw,"gamma");
		for (j = 0; j < NG; j++){		
		    gsl_ran_multinomial(r,maxBSG,1,blockSizeProbG,vecBSG);
            blockSize = 0;
	        while(vecBSG[blockSize]==0) blockSize++;
            blockSize += 1;
            nBlocks = ceil((double)vecLG[j]/blockSize);
            gsl_ran_shuffle(r,indexG[j],vecLG[j],sizeof(int));
	        SPC = SPcalc(n,1,tol,yTilde,gamma,Ngamma,LG,ceta,X,LPV,&Q2);
	        for (block = 0; block < nBlocks; block++){
                s = gsl_ran_flat(r,1.0,100000);
	            proposeBlockInd(s,vecGamma[j],vecLG[j],block,blockSize,indexG[j],cmu[j],dmu[j],vecGammaP[j]); 
	            for (k = 0; k < vecLG[j]; k++)
	                gammaP[cusumVecLG[j]+k] = vecGammaP[j][k];	            	            
	            NPJ = 0;
                for (k = 0; k < vecLG[j]; k++)
                    NPJ += vecGammaP[j][k];
	            NgammaP = Ngamma - vecNgamma[j] + NPJ; 
	            SPP = SPcalc(n,1,tol,yTilde,gammaP,NgammaP,LG,ceta,X,LPV,&Q2);
                Acp = exp((-SPP+SPC)/(2*sigma2))*pow(ceta+1,0.5*(vecNgamma[j]-NPJ));
                unifRV = gsl_ran_flat(r,0.0,1.0);            
	            if (Acp > unifRV){
                    for (k = 0; k < vecLG[j]; k++){
	                    gamma[cusumVecLG[j]+k] = gammaP[cusumVecLG[j]+k];
	                    vecGamma[j][k] = vecGammaP[j][k];
					}
                    Ngamma = NgammaP; 
                    vecNgamma[j] = NPJ;
                    SPC = SPP;
                } 
                else{  
                    for (k = 0; k < vecLG[j]; k++)
	                    gammaP[cusumVecLG[j]+k] = gamma[cusumVecLG[j]+k];
	                 NgammaP = Ngamma;
				}				
            }  	    	   
	    }    
        
	    // - 2 - Update delta and alpha
        //Rprintf("%i %s \n",sw,"delta & alpha");
        subMeanEta = gsl_vector_subvector(meanEta,0,Ngamma+1);
        subVarEta = gsl_matrix_submatrix(varEta,0,0,Ngamma+1,Ngamma+1);
        postMeanVarEta(n,1,tol,gamma,Ngamma,LG,sigma2,ceta,LPV,X,yTilde,&subMeanEta.vector,&subVarEta.matrix,sw);
        cSqRes(n,1,gamma,Ngamma,LG,X,&subMeanEta.vector,Y,sqResC);       
        for (j = 0; j < ND; j++){
		    if ((sw % batchL)==0 && WB > 0){ 
	            if (acceptAlpha[j] > 0.25 && h[j] < HUL) h[j] += MIN(0.01,1/sqrt(WB)) * h[j]; 
	            if (acceptAlpha[j] <= 0.25 && h[j] > HLL) h[j] -= MIN(0.01,1/sqrt(WB)) * h[j];
	            acceptAlpha[j] = 0.0;   
	            if (h[j] < HLL) h[j] = HLL;
	            if (h[j] > HUL) h[j] = HUL;
            }            
		    gsl_ran_multinomial(r,maxBSD,1,blockSizeProbD,vecBSD);
            blockSize = 0;
            while(vecBSD[blockSize]==0) blockSize++;
            blockSize += 1;
            nBlocks = ceil((double)vecLD[j]/blockSize);
            gsl_ran_shuffle(r,indexD[j],vecLD[j],sizeof(int)); 
            for (block = 0; block < nBlocks; block++){  
                s = gsl_ran_flat(r,1.0,100000);
                proposeBlockInd(s,vecDelta[j],vecLD[j],block,blockSize,indexD[j],csigma[j],dsigma[j],vecDeltaP[j]);            
                NPJ = 0;  
                for (k = 0; k < vecLD[j]; k++)
                    NPJ += vecDeltaP[j][k];
                NdeltaP = Ndelta - vecNdelta[j] + NPJ;
                for (k = 0; k < vecLD[j]; k++)
	                deltaP[cusumVecLD[j]+k] = vecDeltaP[j][k]; 	       	            	           	             	                            
                if (NPJ > 0){                
                    subAlphaHat = gsl_vector_subvector(alphaHat,0,NPJ);
                    subD = gsl_matrix_submatrix(D,0,0,NPJ,NPJ);
                    DeltaAlphaHat(n,1,tol,LPV,sqResC,deltaP,NPJ,cusumVecLD[j],cusumVecLD[j+1],
                                  Z,sigma2,sigma2z,calpha,&subD.matrix,&subAlphaHat.vector);              
                    subAlphaP = gsl_vector_subvector(alphaP,0,NPJ);
                    gsl_matrix_scale(&subD.matrix,h[j]);
                    s = gsl_ran_flat(r,1.0,100000);
                    sampleMN(s,NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);
                    logMVNormP = logMVNormalpdf(NPJ,&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol); 		                                                                       
                }else logMVNormP = 0.0;                                             
                move = 0;
                for (k = 0; k < vecLD[j]; k++){
                    if (vecDeltaP[j][k]==1) alphaPD[cusumVecLD[j]+k] = gsl_vector_get(&subAlphaP.vector,move++);                
                    else alphaPD[cusumVecLD[j]+k] = 0;
                } 
                for (i = 0; i < n; i++){
                    LPVP[i] = LPV[i];
                    for (k = cusumVecLD[j]; k < cusumVecLD[j+1]; k++)
                        LPVP[i] += (alphaPD[k]-alpha[k])*Z[n+k*n+i];
                    sigma2zP[i] = sigma2*exp(LPVP[i]);
                    yTildeP[i] = exp(-LPVP[i]/2)*Y[i];
                }   
                QFP = QFC;
                for (k = cusumVecLD[j]; k < cusumVecLD[j+1]; k++)
                    QFP += pow(alphaPD[k],2)-pow(alpha[k],2);
                detR = 0.0;
                for (i = 0; i < n; i++)
                    detR += LPV[i] - LPVP[i];
                detR *= 0.5;
	            SPP = SPcalc(n,1,tol,yTildeP,gamma,Ngamma,LG,ceta,X,LPVP,&Q2);     
	            //probability of reverse direction
	            postMeanVarEta(n,1,tol,gamma,Ngamma,LG,sigma2,ceta,LPVP,X,yTildeP,&subMeanEta.vector,&subVarEta.matrix,sw);
                cSqRes(n,1,gamma,Ngamma,LG,X,&subMeanEta.vector,Y,sqResP);	    	                             
                if (vecNdelta[j] > 0){                
	                subAlphaHat = gsl_vector_subvector(alphaHat,0,vecNdelta[j]);
                    subD = gsl_matrix_submatrix(D,0,0,vecNdelta[j],vecNdelta[j]);
                    DeltaAlphaHat(n,1,tol,LPVP,sqResP,delta,vecNdelta[j],cusumVecLD[j],cusumVecLD[j+1],
                                  Z,sigma2,sigma2zP,calpha,&subD.matrix,&subAlphaHat.vector);              
                    gsl_matrix_scale(&subD.matrix,h[j]);
                    move=0;
                    for (k = 0; k < vecLD[j]; k++)
                        if (vecDelta[j][k]==1) BaseSubAlpha[move++] = alpha[cusumVecLD[j]+k];
                    subAlphaP = gsl_vector_view_array(BaseSubAlpha,vecNdelta[j]);
                    logMVNormC = logMVNormalpdf(vecNdelta[j],&subAlphaP.vector,&subAlphaHat.vector,&subD.matrix,tol);                                 
	            }else logMVNormC = 0.0;   
                Acp = exp(detR+(-SPP+SPC)/(2*sigma2)+logMVNormC-logMVNormP+(QFC-QFP)/(2*calpha))*
                      pow(2*M_PI*calpha,0.5*(Ndelta-NdeltaP));
                unifRV = gsl_ran_flat(r,0.0,1.0);           
	            if (Acp > unifRV){	      	   	    
                    if (NPJ > 0) acceptAlpha[j] += 1/((double)batchL);
                    for (k = 0; k < vecLD[j]; k++){
                        delta[cusumVecLD[j]+k] = deltaP[cusumVecLD[j]+k];
                        vecDelta[j][k] = vecDeltaP[j][k];
                        alpha[cusumVecLD[j]+k] = alphaPD[cusumVecLD[j]+k];                  
                    }                
                    for (i = 0; i < n; i++){
                         LPV[i] = LPVP[i];
                         sigma2z[i] = sigma2zP[i];
                         yTilde[i] = yTildeP[i];
                         sqResC[i] = sqResP[i];				
				    }
                    SPC = SPP;
                    QFC = QFP;
                    Ndelta =  NdeltaP;
                    vecNdelta[j] = NPJ;	      	     
                }
                else{
					for (k = 0; k < vecLD[j]; k++){
                        deltaP[cusumVecLD[j]+k] = delta[cusumVecLD[j]+k];
					    alphaPD[cusumVecLD[j]+k] = alpha[cusumVecLD[j]+k];
					}
				}                      
            }	    
	    }
        
        // - 3 - sigma2
        //Rprintf("%i %s \n",sw,"sigma2");
		if (HN==0) sigma2 = 1/gsl_ran_gamma(r,alphasigma+0.5*n,1/(betasigma+0.5*SPC));
        else if (HN==1){
            if ((sw % batchL)==0 && WB > 0){ 
	            if (acceptSigma2 > 0.25 && f < FUL) f += MIN(0.01,1/sqrt(WB)) * f; 
	            if (acceptSigma2 <= 0.25 && f > FLL) f -= MIN(0.01,1/sqrt(WB)) * f;
	            acceptSigma2 = 0.0;   
	            if (f < FLL) f = FLL;
	            if (f > FUL) f = FUL;
            }            
            sigma2P = sigma2 + gsl_ran_gaussian(r,sqrt(f));
            while (sigma2P < 0) sigma2P = sigma2 + gsl_ran_gaussian(r,sqrt(f));
            Acp = exp(-0.5*n*log(sigma2P/sigma2) + (SPC/2)*(1/sigma2-1/sigma2P) + 
                  (sigma2-sigma2P)/(2*phi2sigma));
	        unifRV = gsl_ran_flat(r,0.0,1.0);
            if (Acp > unifRV){
                sigma2 = sigma2P;
	            acceptSigma2 += 1/((double)batchL);  
                for (i = 0; i < n; i++)
                    sigma2z[i] = sigma2*exp(LPV[i]);
            }
	    }
	    
        // - 4 - c_eta 
	    //Rprintf("%i %s \n",sw,"c_eta");
        if ((sw % batchL)==0 && WB > 0){ 
	        if (acceptCeta > 0.25 && g < GUL) g += MIN(0.01,1/sqrt(WB)) * g; 
	        if (acceptCeta <= 0.25 && g > GLL) g -= MIN(0.01,1/sqrt(WB)) * g;
	        acceptCeta = 0.0;   
	        if (g < GLL) g = GLL;
	        if (g > GUL) g = GUL;
        }     
        cetahat = 1;//starting value for NR
        SPC = SPcalc(n,1,tol,yTilde,gamma,Ngamma,LG,ceta,X,LPV,&Q2);
        elPrime = 99.9;
        while(elPrime > 0.000000001 || -elPrime > 0.000000001){
            Sprime = -Q2/pow(cetahat+1,2);
            SDprime = 2*Q2/pow(cetahat+1,3);
            elPrime = -0.5*(Ngamma+1)/(cetahat+1)-Sprime/(2*sigma2)-(alphaeta+1)/cetahat+betaeta/pow(cetahat,2);
            elDPrime = 0.5*(Ngamma+1)/(pow(cetahat+1,2))-SDprime/(2*sigma2)+(alphaeta+1)/pow(cetahat,2)-2*betaeta/pow(cetahat,3);
            cetahat -= elPrime / elDPrime;
        }
        cetaP = cetahat + gsl_ran_gaussian(r,sqrt(-g/elDPrime));
	    while(cetaP < 0) cetaP = cetahat + gsl_ran_gaussian(r,sqrt(-g/elDPrime));
	    SPP = SPcalc(n,1,tol,yTilde,gamma,Ngamma,LG,cetaP,X,LPV,&Q2);
	    Acp = exp(-0.5*(Ngamma+1)*(log(cetaP+1)-log(ceta+1))+(-SPP+SPC)/(2*sigma2)-(alphaeta+1)*(log(cetaP)-log(ceta))
              + betaeta*(1/ceta-1/cetaP))*
                gsl_ran_gaussian_pdf(ceta-cetahat,sqrt(-g/elDPrime))/
                gsl_ran_gaussian_pdf(cetaP-cetahat,sqrt(-g/elDPrime));
	    unifRV = gsl_ran_flat(r,0.0,1.0);
        if (Acp > unifRV){
            ceta = cetaP;
	        acceptCeta += 1/((double)batchL);  
	    }
	    
	    // - 5 - c_alpha
	    //Rprintf("%i %s \n",sw,"c_alpha");
        calpha = 1/gsl_ran_gamma(r,alphaalpha+0.5*Ndelta,1/(betaalpha+0.5*QFC));
        
        if (((sw - burn) >= 0) && (((sw - burn ) % thin) == 0) && (WF == 1)){
            // - 6 - eta              
            subMeanEta = gsl_vector_subvector(meanEta,0,Ngamma+1);//not needed
            subVarEta = gsl_matrix_submatrix(varEta,0,0,Ngamma+1,Ngamma+1);//when all code is running
            postMeanVarEta(n,1,tol,gamma,Ngamma,LG,sigma2,ceta,LPV,X,yTilde,&subMeanEta.vector,&subVarEta.matrix,sw);
            vecEta = gsl_vector_view_array(eta,Ngamma+1); 
            s = gsl_ran_flat(r,1.0,100000);
            sampleMN(s,Ngamma+1,&vecEta.vector,&subMeanEta.vector,&subVarEta.matrix,tol);
            // write to files
            for (k = 0; k < LG; k++)
                fprintf(out_file1, "%i ", gamma[k]);
            fprintf(out_file1, "\n");
            for (k = 0; k < LD; k++)
                fprintf(out_file2, "%i ", delta[k]);
            fprintf(out_file2, "\n");
            for (k = 0; k < LD; k++) 
                fprintf(out_file3, "%f ", alpha[k]);
            fprintf (out_file3, "\n");
            fprintf(out_file4, "%f \n", sigma2);
            fprintf(out_file5, "%f \n", ceta);
            fprintf(out_file6, "%f \n", calpha);
            move = 0;
            fprintf(out_file7, "%f ", eta[move++]);
            for (k = 0; k < LG; k++)
                if (gamma[k]==1) fprintf(out_file7, "%f ", eta[move++]); else fprintf(out_file7, "%f ", 0.0);                
            fprintf (out_file7, "\n");
        }
        // If sw needs to be printed
        if (sw==(sweeps-1) && (!(sweeps % 1000)==0)) Rprintf("%i %s \n",sweeps, "posterior samples...");
    }//end of sw
    //Update adaptive parameters
    f1[0] = f;
    g1[0] = g;

    //Free up random number generator
    gsl_rng_free (r);

    //Close files
    fclose(out_file1); fclose(out_file2); fclose(out_file3);
    fclose(out_file4); fclose(out_file5); fclose(out_file6);
    fclose(out_file7);

    //Free up gsl matrices
    gsl_matrix_free(varEta); gsl_vector_free(meanEta); gsl_matrix_free(D); gsl_vector_free(alphaHat); 
    gsl_vector_free(alphaP);
	
	//Free jagged vectors
	for (j = 0; j < NG; j++) {free(indexG[j]); free(vecGamma[j]); free(vecGammaP[j]);}
	for (j = 0; j < ND; j++) {free(indexD[j]); free(vecDelta[j]); free(vecDeltaP[j]);}
}
