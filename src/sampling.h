/* sampling.h Includes sampling methods
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

#include <math.h>

//Sample from a univariate truncated normal distribution. Arguments: seed, subject index,
//mean and standard deviation, lower and upper bounds, store latent
void sampleTUN(unsigned long int s, int i, double tnmean, double tnsd, double lower, double upper, double *latentx){
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    if ((lower - upper) == 0)
        latentx[i] = lower;
    else{
        latentx[i] = lower - 1.0;
        if ((lower < tnmean) && (upper < tnmean)){
            while ((trunc(pow(10,10)*(latentx[i] - lower)) < 0) || (trunc(pow(10,10)*(latentx[i] - upper)) > 0))
                latentx[i] = tnmean - gsl_ran_gaussian_tail(r,-upper+tnmean,tnsd);
        }
        else if ((lower > tnmean) && (upper > tnmean)){
            while ((trunc(pow(10,10)*(latentx[i] - lower)) < 0) || (trunc(pow(10,10)*(latentx[i] - upper)) > 0))
                latentx[i] = tnmean + gsl_ran_gaussian_tail(r,lower-tnmean,tnsd);
        }
        else{
            while ((trunc(pow(10,10)*(latentx[i] - lower)) < 0) || (trunc(pow(10,10)*(latentx[i] - upper)) > 0))
                latentx[i] = tnmean+gsl_ran_gaussian(r,tnsd);
        }
    }
    gsl_rng_free(r);
}

//Sample from a multivariate truncated normal distribution. Arguments: seed, dim, lower bounds,
//upper bounds, sample from previous iteration, mean, covariance, tolerance.
//Samples are stored in gsl vector y. Uses Robert (2009) Simulation of truncated normal variables, 
//and formula (3.2) for inverting the matrices. 
void sampleTMN(unsigned long int s, int p, double *L, double *U, gsl_vector *y, gsl_vector *m,
               gsl_matrix *Sigma, double tol){
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    int i;
    double res1, res2, res3, temp;
    double cm[2*p];
    gsl_vector *V1, *V2;
    gsl_matrix *cSigma, *V;
    V1 = gsl_vector_alloc(p);
    V2 = gsl_vector_alloc(p);
    cSigma = gsl_matrix_alloc(p,p);
    V = gsl_matrix_alloc(p,p);
    gsl_matrix_memcpy(V,Sigma);
    ginv(p,tol,V);
    //Inverse(p,V);
    gsl_vector_view Sigmai;
    gsl_vector_view Vi;
    for (i=0; i < p; i++){
        gsl_vector_memcpy(V1,y);
        gsl_matrix_memcpy(cSigma,Sigma);
        gsl_vector_sub(V1,m);
        gsl_vector_set(V1,i,0);
        gsl_blas_dgemv(CblasNoTrans,1.0,V,V1,0.0,V2);
        gsl_vector_set(V2,i,0);
        Sigmai = gsl_matrix_row(cSigma,i);
        gsl_vector_set(&Sigmai.vector,i,0);
        gsl_blas_ddot(&Sigmai.vector,V2,&res1);
        Vi = gsl_matrix_row(V,i);
        gsl_blas_ddot(&Sigmai.vector,&Vi.vector,&res2);
        gsl_blas_ddot(&Vi.vector,V1,&res3);
        cm[i] = gsl_vector_get(m,i) + res1 - res2*res3/gsl_matrix_get(V,i,i);
        gsl_blas_dgemv(CblasNoTrans,1.0,V,&Sigmai.vector,0.0,V1);
        gsl_blas_ddot(&Sigmai.vector,V1,&res1);
        cm[i+p] = gsl_matrix_get(Sigma,i,i) - res1 + res2*res2/gsl_matrix_get(V,i,i);
        temp = L[i]-0.01; // temp replaces y[i] to avoid using gsl get and set
        if (fabs(L[i] - U[i]) < 0.0001) temp = L[i]/2 + U[i]/2;
        else {
            if ((L[i] < cm[i]) && (U[i] < cm[i])){
                while ((trunc(pow(10,10)*(temp - L[i])) < 0) || (trunc(pow(10,10)*(temp - U[i])) > 0))
                    temp = cm[i]-gsl_ran_gaussian_tail(r,-U[i]+cm[i],sqrt(cm[i+p]));
            }
            else if ((L[i] > cm[i]) && (U[i] > cm[i])){
                while ((trunc(pow(10,10)*(temp - L[i])) < 0) || (trunc(pow(10,10)*(temp - U[i])) > 0))
                    temp = cm[i]+gsl_ran_gaussian_tail(r,L[i]-cm[i],sqrt(cm[i+p]));
            }
            else{
                while ((trunc(pow(10,10)*(temp - L[i])) < 0) || (trunc(pow(10,10)*(temp - U[i])) > 0))
                    temp = cm[i]+gsl_ran_gaussian(r,sqrt(cm[i+p]));
            }
        }
        gsl_vector_set(y,i,temp);
    }
    gsl_vector_free(V1);
    gsl_vector_free(V2);
    gsl_matrix_free(cSigma);
    gsl_matrix_free(V);
    gsl_rng_free(r);
}

//Sample from a multivariate truncated normal distribution, with lower bound fixed at zero and no upper bound. 
//Arguments: seed, dimension, gsl vector for storing samples, mean, covariance, tolerance.
//Uses Robert (2009) Simulation of truncated normal variables, and formula (3.2) for inverting the matrices.
void sampleTMN2(unsigned long int s, int p, gsl_vector *y, gsl_vector *m, gsl_matrix *Sigma, double tol){
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    int i, k;
    double res1, res2, res3, temp;
    double cm[2*p];
    gsl_vector *V1, *V2;
    gsl_matrix *cSigma, *V;
    V1 = gsl_vector_alloc(p);
    V2 = gsl_vector_alloc(p);
    cSigma = gsl_matrix_alloc(p,p);
    V = gsl_matrix_alloc(p,p);
    gsl_matrix_memcpy(V,Sigma);
    ginv(p,tol,V);    
    gsl_vector_view Sigmai;
    gsl_vector_view Vi;
    //initialize vector y to be the mean or if mean is out of bounds, set it to zero
    for (i=0; i < p; i++){
		temp = gsl_vector_get(m,i); 
		if (temp < 0) temp = 0;
		gsl_vector_set(y,i,temp);
	}
    for (k=0; k < 3; k++){ //repeat 3 times, the first 2 are burn-in
        for (i=0; i < p; i++){
            gsl_vector_memcpy(V1,y);
            gsl_matrix_memcpy(cSigma,Sigma);
            gsl_vector_sub(V1,m);
            gsl_vector_set(V1,i,0);
            gsl_blas_dgemv(CblasNoTrans,1.0,V,V1,0.0,V2);
            gsl_vector_set(V2,i,0);
            Sigmai = gsl_matrix_row(cSigma,i);
            gsl_vector_set(&Sigmai.vector,i,0);
            gsl_blas_ddot(&Sigmai.vector,V2,&res1);
            Vi = gsl_matrix_row(V,i);
            gsl_blas_ddot(&Sigmai.vector,&Vi.vector,&res2);
            gsl_blas_ddot(&Vi.vector,V1,&res3);
            cm[i] = gsl_vector_get(m,i) + res1 - res2*res3/gsl_matrix_get(V,i,i);
            gsl_blas_dgemv(CblasNoTrans,1.0,V,&Sigmai.vector,0.0,V1);
            gsl_blas_ddot(&Sigmai.vector,V1,&res1);
            cm[i+p] = gsl_matrix_get(Sigma,i,i) - res1 + res2*res2/gsl_matrix_get(V,i,i);
            if (cm[i] > 0){
                temp = cm[i]+gsl_ran_gaussian(r,sqrt(cm[i+p]));
                while (temp < 0)
                    temp = cm[i]+gsl_ran_gaussian(r,sqrt(cm[i+p]));	
            }
            else{
			    temp = cm[i]+gsl_ran_gaussian_tail(r,-cm[i],sqrt(cm[i+p]));
                while (temp < 0)
                    temp = cm[i]+gsl_ran_gaussian_tail(r,-cm[i],sqrt(cm[i+p]));
            }
            gsl_vector_set(y,i,temp);
        }
	}
    gsl_vector_free(V1);
    gsl_vector_free(V2);
    gsl_matrix_free(cSigma);
    gsl_matrix_free(V);
    gsl_rng_free(r);
}
 
//Sample from a multivariate normal distribution. Arguments: seed, dim, vector to store sample, mean, covariance, tolerance.
//Samples are stored in gsl vector y.
void sampleMN(unsigned long int s, int p, gsl_vector *y, gsl_vector *mu, gsl_matrix *Sigma, double tol){
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    int i;
    gsl_matrix *SigmaHalf = gsl_matrix_alloc(p,p);
    gsl_vector *temp = gsl_vector_alloc(p);
    gsl_matrix_memcpy(SigmaHalf,Sigma);
    matHalf(p,tol,SigmaHalf);
    for (i=0; i < p; i++) // N(0,1)
	    gsl_vector_set(temp,i,gsl_ran_gaussian(r,1));
    gsl_blas_dgemv(CblasNoTrans,1.0,SigmaHalf,temp,0.0,y);
    gsl_vector_add(y,mu);
    gsl_matrix_free(SigmaHalf);
    gsl_vector_free(temp);
    gsl_rng_free(r);
}

//Sample from a Wishart distribution. Takes as input seed s, dimension p, degrees of freedom n, scale matrix scale, random
//variate rw. The random variate is stored in rw.
void rwish(unsigned long int s, int p, double n, gsl_matrix *scale, gsl_matrix *rw){
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    int i,j;

    // Declare matrices
    gsl_matrix *CC, *Z, *Prod1;
    CC = gsl_matrix_alloc(p,p);
    Z = gsl_matrix_calloc(p,p);
    Prod1 = gsl_matrix_alloc(p,p);

    // Set matrix CC equal to the scale matrix
    gsl_matrix_memcpy(CC,scale);

    // Cholesky decomposition
    //gsl_set_error_handler_off();
    gsl_linalg_cholesky_decomp(CC);
    //gsl_set_error_handler(NULL);

    // Make CC lower triangular (gsl returns full matrix)
    for (i = 1; i < p; i++)
         for (j = 0; j < i; j++)
           gsl_matrix_set(CC,i,j,0.0);

    // Continue with matrix Z which require Chi-sq random variates
    for (i = 0; i < p; i++)
       gsl_matrix_set(Z, i, i, sqrt(gsl_ran_chisq(r,(n-i))));

    for (i = 0; i < (p-1); i++)
       for (j = (i+1); j < p; j++)
          gsl_matrix_set(Z,i,j,gsl_ran_ugaussian(r));

    // Finally, the random variate is given by crossprod(Z %*% CC) = t(Z%*%CC) %*% Z%*%CC
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Z,CC,0.0,Prod1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Prod1,Prod1,0.0,rw);

    // Free up gsl objects
    gsl_matrix_free(CC);
    gsl_matrix_free(Z);
    gsl_matrix_free(Prod1);
    gsl_rng_free(r);
}

//Given seed, sample size, number of components, (a,b) prior parameters, number of subjects per cluster, and current value
//of the concentration parameter, the function returns a new value for the concentration paramter
double updateAlpha(unsigned long int s, int n, int ncomp, double a, double b, double TruncAlpha, int *nmembers, double alpha){
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    int h;
    double eta, pieta, unifRV, newalpha, K;
    newalpha = 0.0;
    eta = gsl_ran_beta(r,alpha+1.0,(double) n);
    K = 0.0;
    for (h = 0; h < ncomp; h++)
        if (nmembers[h] > 0) K++;
    pieta = (a+K-1.0)/(a+K-1.0+n*(b-log(eta)));
    unifRV = gsl_ran_flat(r,0.0,1.0);
    newalpha = 0.0;
    while (newalpha < TruncAlpha){
        if (unifRV < pieta)
            newalpha = gsl_ran_gamma(r,a+K,(1.0/(b-log(eta))));
        else
            newalpha = gsl_ran_gamma(r,a+K-1.0,(1.0/(b-log(eta))));
    }
    gsl_rng_free(r);
    return(newalpha);
}

//Given component  allocation and eta that plays the role of the mean, update z: for spatial only
void updatez(unsigned long int s, int n, int ncomp, int *compAlloc, double eta[ncomp][n], double z[ncomp][n]){
    int j,h;
    double muZ;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    for (j = 0; j < n; j++){
        for (h = 0; h < ncomp; h++){
            muZ = eta[h][j];
            if ((h < compAlloc[j]) && (muZ > 0)) z[h][j] = muZ - gsl_ran_gaussian_tail(r, muZ, 1.0);
            if ((h < compAlloc[j]) && (muZ < 0)){
                z[h][j] = 10.0;
                while (z[h][j] > 0) {z[h][j] = muZ + gsl_ran_gaussian(r,1);}
            }
            if ((h == compAlloc[j]) && (muZ > 0)){
                z[h][j] = -10.0;
                while (z[h][j] < 0) {z[h][j] = muZ + gsl_ran_gaussian(r,1);}
            }
            if ((h == compAlloc[j]) && (muZ < 0)) z[h][j] = muZ + gsl_ran_gaussian_tail(r, -muZ, 1.0);;
            if (h > compAlloc[j]) z[h][j] = muZ + gsl_ran_gaussian(r,1);
        }
    }
    gsl_rng_free(r);
}

//Impute GMRF
void imputeGMRF(unsigned long int s, int n, int ncomp, double alphau, double phiu, double lu, double *eigenvl, gsl_matrix *qij,
                double z[ncomp][n], double eta[ncomp][n]){

    int k, h, i;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    gsl_matrix *dgnl = gsl_matrix_calloc(n,n);
    gsl_matrix *WM = gsl_matrix_alloc(n,n);
    gsl_matrix *WM2 = gsl_matrix_alloc(n,n);
    gsl_vector *W = gsl_vector_alloc(n);
    gsl_vector *W1 = gsl_vector_alloc(n);
    gsl_vector *V = gsl_vector_alloc(n);

    for (k = 0; k < n; k++)
        gsl_matrix_set(dgnl,k,k,1/sqrt(phiu*phiu*lu*eigenvl[k]+phiu*phiu+1.0));

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,qij,dgnl,0.0,WM2);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,WM2,qij,0.0,WM);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,WM,WM,0.0,WM2);

    //Impute Gaussian Markov random fields
    for (h = 0; h < ncomp; h++){
        for (i = 0; i < n; i++) // N(0,1)
	        gsl_vector_set(W,i,gsl_ran_gaussian(r,1));
        for (i = 0; i < n; i++)
            gsl_vector_set(V,i,alphau*phiu*phiu+z[h][i]);
        gsl_blas_dgemv(CblasNoTrans,1.0,WM,W,0.0,W1);
        gsl_blas_dgemv(CblasNoTrans,1.0,WM2,V,1.0,W1);
        for (i = 0; i < n; i++)
            eta[h][i] = gsl_vector_get(W1,i);
    }
    gsl_matrix_free(dgnl); gsl_matrix_free(WM); gsl_matrix_free(WM2);
    gsl_vector_free(W); gsl_vector_free(W1); gsl_vector_free(V);
    gsl_rng_free(r);
}

double updatespatialalpha(unsigned long int s, int n, int ncomp, int *nmembers, double phiu, double mualpha, double sigmalpha,
                          double eta[ncomp][n], double nClusters){
    int h, i;
    double sumeta, alphau;

    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    sumeta = 0.0;

    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0){
            for (i = 0; i < n; i++)
                sumeta += eta[h][i];
        }
    }

    alphau = (phiu*phiu*sumeta+mualpha/(sigmalpha*sigmalpha))/(n*nClusters*phiu*phiu+1/(sigmalpha*sigmalpha)) +
             gsl_ran_gaussian(r,1/sqrt(n*nClusters*phiu*phiu+1/(sigmalpha*sigmalpha)));

    gsl_rng_free(r);
    return(alphau);
}

double updatespatialphiu(unsigned long int s, int n, double lu, double alphaphi,
                         double betaphi, double nClusters, double *BS){
    double BS1, BS2, phisq, phiu;

    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    BS1 = BS[0]; BS2 = BS[1];

    phisq = gsl_ran_gamma(r, alphaphi+n*nClusters/2.0, 1.0/(betaphi+0.5*(BS1+lu*BS2)));
    phiu = sqrt(phisq);

    gsl_rng_free(r);
    return(phiu);
}

double updatespatiallu(unsigned long int s, int n, double *eigenvl, double lu, double steplu, double phiu,
                       double Mlu, double nClusters, double *BS){
    int i;
    double proplu, Detr, logRt, Rt, ludecision, BS2;

    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    BS2 = BS[1];

    proplu = lu + gsl_ran_gaussian(r,steplu);

    if ((proplu < Mlu) && (proplu > 0.0)){
        Detr = 0.0;
        for (i = 0; i < n; i++)
            Detr += log(eigenvl[i]*proplu+1)-log(eigenvl[i]*lu+1);
        logRt = Detr*nClusters/2-(phiu*phiu/2)*BS2*(proplu-lu);
        Rt = exp(logRt);
        if (Rt > 1.0) Rt = 1.0;
        ludecision = gsl_ran_flat(r, 0.0, 1.0);
        if (ludecision < Rt) lu = proplu;
    }
    gsl_rng_free(r);
    return(lu);
}

//obtained proposed values for the One Res Ltnt model.
void propose(unsigned long int s, double *XiC, double *XiP, int nRespPars, int j, double *prec, int family){
    int k;
    double alpha, beta, min;
    double L = 9999.99;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    if (family == 1)
        XiP[0] = gsl_ran_gamma(r,prec[0]*XiC[0]*XiC[0],1/(prec[0]*XiC[0]));
    else if (family == 2){
        beta = XiC[0] - 1 + XiC[0]*(1-XiC[0])*(1-XiC[0])*prec[0];
        if (beta < 0.001) beta = 0.001;
        alpha = beta * XiC[0]/(1-XiC[0]);
        XiP[0] = gsl_ran_beta(r,alpha,beta);
    }
    else if (family == 3 || family == 4)
        for (k = 0; k < nRespPars; k++)
            XiP[k] = gsl_ran_gamma(r,prec[k]*XiC[k]*XiC[k],1/(prec[k]*XiC[k]));
    else if (family == 5){
        //Rprintf("%s %f %f \n","FIC",XiC[0],XiC[1]);
        if (j == 0){ 
            XiP[0] = gsl_ran_gamma(r,prec[0]*XiC[0]*XiC[0],1/(prec[0]*XiC[0]));
            XiP[1] = XiC[1];
		}
        else if (j == 1){
            XiP[0] = XiC[0]; 
            min = 0.05;
            //if (1-XiC[0]/4 > min) min = 1-XiC[0]/4; 
            //XiP[1] = XiC[1] + gsl_ran_gaussian(r,1/sqrt(prec[1])); gsl_ran_gaussian_tail(r,,1/sqrt(prec[1]))
            //while(XiP[1] < 0.5 || XiP[1] < (1-XiP[0]/4)) XiP[1] = XiC[1] + gsl_ran_gaussian(r,);
            sampleTUN(s,1,XiC[1],1/sqrt(prec[1]),min,L,XiP);        
        }
        //Rprintf("%s %f %f \n","FIP",XiP[0],XiP[1]);
    }
    else if (family == 6){
        for (k = 0; k < nRespPars; k++)
            XiP[k] = gsl_ran_gamma(r,prec[k]*XiC[k]*XiC[k],1/(prec[k]*XiC[k]));
        while((XiP[1] < 0.3) || (XiP[1] < 1/(1+2*XiP[0]))) XiP[1] = gsl_ran_gamma(r,prec[1]*XiC[1]*XiC[1],1/(prec[1]*XiC[1]));
    }
    else if (family == 7){
        for (k = 0; k < nRespPars; k++)
            XiP[k] = gsl_ran_gamma(r,prec[k]*XiC[k]*XiC[k],1/(prec[k]*XiC[k]));
        while(XiP[1] < 0.1) XiP[1] = gsl_ran_gamma(r,prec[1]*XiC[1]*XiC[1],1/(prec[1]*XiC[1]));
    }
    else if (family == 8){
        do {
            for (k = 0; k < (nRespPars-1); k++)
                XiP[k] = gsl_ran_gamma(r,prec[k]*XiC[k]*XiC[k],1/(prec[k]*XiC[k]));
            XiP[2] = XiC[2] + gsl_ran_gaussian(r,1/prec[2]);
            if (XiC[2] <= (XiP[1]/2-1)) while(XiP[2] > (XiP[1]/2-1)) XiP[2] = XiC[2] + gsl_ran_gaussian(r,1/prec[2]);
            else if (XiC[2] > (XiP[1]/2-1)) XiP[2] = XiC[2] - gsl_ran_gaussian_tail(r,-(XiP[1]/2-1)+XiC[2],1/prec[2]);
        } while ((XiP[0]*(XiP[1]-2*XiP[2]-1)-XiP[2]*XiP[2]) < 0);
    }
    for (k = 0; k < nRespPars; k++) 
        if (XiP[k] < 0.00001) XiP[k] = 0.00001;
    gsl_rng_free(r);
}

//obtained proposed values for the One Res Ltnt model. prec is the precision around the current estimate
void propose2(unsigned long int s, double *XiC, double *XiP, int nRespPars, double *prec, int family){
    int k;
    double alpha, beta;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    if (family == 1)
        XiP[0] = gsl_ran_gamma(r,prec[0]*XiC[0]*XiC[0],1/(prec[0]*XiC[0])); //has mean XiC[0] and variance=1/prec
    else if (family == 2){
        beta = XiC[0] - 1 + XiC[0]*(1-XiC[0])*(1-XiC[0])*prec[0];
        if (beta < 0.001) beta = 0.001;
        alpha = beta * XiC[0]/(1-XiC[0]);
        XiP[0] = gsl_ran_beta(r,alpha,beta); // has mean XiC[0] and variance=1/prec
    }
    else if (family == 3 || family == 4)
        for (k = 0; k < nRespPars; k++)
            XiP[k] = gsl_ran_gamma(r,prec[k]*XiC[k]*XiC[k],1/(prec[k]*XiC[k]));
    else if (family == 5){
        XiP[0] = XiC[0] + gsl_ran_gaussian(r,1/prec[0]);
        XiP[1] = XiC[1] + gsl_ran_gaussian(r,1/prec[1]);
        XiP[2] = XiC[2] + gsl_ran_gaussian(r,1/prec[2]);
        while(XiP[2] < 0.5) XiP[2] = XiC[2] + gsl_ran_gaussian(r,1/prec[2]);
    }
    else if (family == 6){
        for (k = 0; k < nRespPars; k++)
            XiP[k] = gsl_ran_gamma(r,prec[k]*XiC[k]*XiC[k],1/(prec[k]*XiC[k]));
        while(XiP[1] < 0.3) XiP[1] = gsl_ran_gamma(r,prec[1]*XiC[1]*XiC[1],1/(prec[1]*XiC[1]));
    }
    else if (family == 7){
        for (k = 0; k < nRespPars; k++)
            XiP[k] = gsl_ran_gamma(r,prec[k]*XiC[k]*XiC[k],1/(prec[k]*XiC[k]));
        while(XiP[1] < 0.1) XiP[1] = gsl_ran_gamma(r,prec[1]*XiC[1]*XiC[1],1/(prec[1]*XiC[1]));
    }
    else if (family == 8){
        //Rprintf("%s %f %f %f ","cur:",XiC[0],XiC[1],XiC[2]);
        for (k = 0; k < (nRespPars-1); k++)
            XiP[k] = gsl_ran_gamma(r,prec[k]*XiC[k]*XiC[k],1/(prec[k]*XiC[k]));
        XiP[2] = XiC[2] + gsl_ran_gaussian(r,1/prec[2]);
        //Rprintf("%s %f %f %f ","prp:",XiP[0],XiP[1],XiP[2]);
        while(XiP[2] > (XiP[1]/2-1)) XiP[2] = XiC[2] + gsl_ran_gaussian(r,1/prec[2]);
        //Rprintf("%s %f ","prp2:",XiP[2]);
    }
    gsl_rng_free(r);
}

// from (4) of Wood (1994): dimension m = 2, lambda, ksi, sample
void rvMF(unsigned long int s, int m, double lambda, double *mode, double *out)
{
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    double d = m - 1;
    double b = d / (sqrt(4 * lambda * lambda + d * d) + 2 * lambda); 
    double x = (1 - b) / (1 + b);
    double c = lambda * x + d * log(1 - x * x);
    double z, w, u;
    double test = 0;
    int rs;
    double y[2];
	while(!test){
	    z = gsl_ran_beta(r, d / 2, d / 2); 
	    w = (1 - (1 + b) * z) / (1 - (1 - b) * z);
	    u = gsl_ran_flat(r, 0.0, 1.0);
	    if(lambda * w + d * log(1 - x * w) - c >= log(u)) 
		    test = 1;	   
	}
	rs = 2 * gsl_ran_bernoulli(r, 0.5) - 1;
	y[0] = sqrt(1 - w * w) * rs;
	y[1] = w;
	out[0] = mode[1] * y[0] + mode[0] * y[1];
	out[1] = -mode[0] * y[0] + mode[1] * y[1];
	gsl_rng_free(r);
}
