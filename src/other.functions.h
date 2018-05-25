/* other.functions.h Includes miscellaneous functions
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

//Allocation of observations to clusters with arguments (seed, number of sampling units, number of components,
//probability matrix ncomp x n, a vector that will store the integers indicating cluster allocation.
void allocation(unsigned long int s, int n, int ncomp, double Prob[ncomp][n], int *compAlloc, int sw){
    int i, h, komp;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double compProbi[ncomp];
    unsigned int vecAlloc[ncomp];

    for (i = 0; i < n; i++){
        for (h = 0; h < ncomp; h++)
            compProbi[h] = Prob[h][i];

        //if (sw>3 && i==3){
        //    Rprintf("%s %i %i \n","comp prob: ", i,sw);
        //    for (h = 0; h < ncomp; h++) Rprintf("%f ",compProbi[h]);
        //    Rprintf("\n");
        //}

        gsl_ran_multinomial(r,ncomp,1,compProbi,vecAlloc);
        komp = 0;
        while(vecAlloc[komp]==0) komp++;
        compAlloc[i] = komp;

        //double temp;
        //if (komp >= ncomp) {
            //temp = 0.0;
            //Rprintf("%i\n", i);
            //for (h = 0; h < ncomp; h++) Rprintf("%f ",compProbi[h]);
            //for (h = 0; h < ncomp; h++) temp += compProbi[h];
            //Rprintf("%s %f ","total: ",temp*100);
            //if (temp==0.0) Rprintf("%s %i %i","zero sw i: ",sw,i);
            //Rprintf("\n");
        //}

    }
    gsl_rng_free(r);
}

//Find the number of sampling units per cluster. Arguments: number of sampling units, number of components,
//vector that will contain the number of subjects per cluster, and component allocations
void findNmembers(int n, int ncomp, int *compAlloc, int *nmembers){
    int i, h;
    for (h = 0; h < ncomp; h++) nmembers[h] = 0;
    for (h = 0; h < ncomp; h++)
        for (i = 0; i < n; i++)
            if (compAlloc[i]==h) nmembers[h] += 1;
}

//Logit function
double logit(double p){
    return log(p/(1-p));
}

//Inverse logit function
double invlogit(double p){
    return exp(p)/(1+exp(p));
}

//Initialize Xi coefficients. For Poisson and Binomial mixtures they are set equal to the mean of the cluster;
//for empty clusters, they are set randomly to mean + error.
void initGLMOneResLtnt1(unsigned long int s, int *Y, double *H, int n, int ncomp, int nRespPars,
                        int *nmembers, int *compAlloc, double Xi[ncomp][nRespPars], double Ymean, int family){
    double sumhY1, sumhH;
    int i, h;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0){
            sumhY1 = 0.0;
            sumhH = 0.0;
            for (i = 0; i < n; i++){
                if (compAlloc[i]==h){
                    sumhY1 += (double) Y[i];
                    sumhH += H[i];
                }
            }
            if (sumhY1 > 0)
                Xi[h][0] = sumhY1/sumhH;
            else
                Xi[h][0] = 1.0/sumhH;
        }
        else if (nmembers[h] == 0){
                 if (family == 1 || family > 2) Xi[h][0] = exp(log(Ymean) + gsl_ran_gaussian(r,1.0));
                 if (family == 2) Xi[h][0] = invlogit(logit(Ymean) + gsl_ran_gaussian(r,1.0));
        }
        if (family==8) Xi[h][0] = sqrt(2*Xi[h][0]);
        if (family == 5 || family == 6 || family == 7) Xi[h][1] = 1.0;
        if (family == 8) {Xi[h][1] = 3.0; Xi[h][2] = 0.0;}
    }
    gsl_rng_free(r);
}

//Initialize Xi coefficients for negative binomial and beta binomial mixtures
void initGLMOneResLtnt2(unsigned long int s, int *Y, double *H, int n, int ncomp, int nRespPars,
                         int *nmembers, int *compAlloc, double Xi[ncomp][nRespPars], int family){
    double sumhY1, sumhY2, sumhH;
    double Ybarh, Hbarh, Yvarh;
    double Ymean, Hmean, Yvar;
    int i, h;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    sumhY1 = 0.0;
    sumhY2 = 0.0;
    sumhH = 0.0;
    for (i = 0; i < n; i++){
        sumhY1 += (double) Y[i];
        sumhY2 += (double) Y[i]*Y[i];
        sumhH += H[i];
    }
    Ymean = sumhY1 / n;
    Hmean = sumhH / n;
    Yvar = (sumhY2 - n*Ymean*Ymean)/(n-1);
    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 1){ //problem with variance when n=1
            sumhY1 = 0.0;
            sumhY2 = 0.0;
            sumhH = 0.0;
            for (i = 0; i < n; i++){
                if (compAlloc[i]==h){
                    sumhY1 += (double) Y[i];
                    sumhY2 += (double) Y[i]*Y[i];
                    sumhH += H[i];
                }
            }
            Ybarh = sumhY1 / nmembers[h];
            if (Ybarh == 0.0) Ybarh = 1.0 / nmembers[h];
            Hbarh = sumhH / nmembers[h];
            Yvarh = (sumhY2 - nmembers[h]*Ybarh*Ybarh)/(nmembers[h]-1);
            if (family == 3){
                if (Yvarh==Ybarh) Yvarh += 1.0;
                Xi[h][1] = Hbarh * Ybarh / (Yvarh - Ybarh);
                if (Xi[h][1] < 0.1) Xi[h][1] = 0.1;
                Xi[h][0] = Xi[h][1] * Ybarh / Hbarh;
            }
            else if (family == 4){
                Xi[h][0] = 2;//(-Ybarh*Yvarh+Ybarh*Ybarh*(Hbarh-Ybarh))/(Hbarh*Yvarh-Ybarh*(Hbarh-Ybarh));
                if (Xi[h][0] < 0.1) Xi[h][0] = 0.1;
                Xi[h][1] = Xi[h][0] * (Hbarh/Ybarh-1);
            }
        }
        else if (nmembers[h] < 2){
            if (family == 3){
                Xi[h][1] = Hmean * Ymean / (Yvar - Ymean) * exp(gsl_ran_gaussian(r,1.0));
                if (Xi[h][1] < 0.1) Xi[h][1] = 0.1;
                Xi[h][0] = Xi[h][1] * Ymean / Hmean * exp(gsl_ran_gaussian(r,1.0));
            }
            else if (family == 4){
                Xi[h][0] = (Ymean*Yvar-Ymean*Ymean*(Hmean-Ymean))/(Hmean*Yvar-Ymean*(Hmean-Ymean)) * exp(gsl_ran_gaussian(r,1.0));
                if (Xi[h][0] < 0.1) Xi[h][0] = 0.1;
                Xi[h][1] = Xi[h][0] * (Hmean/Ymean-1) * exp(gsl_ran_gaussian(r,1.0));
            }
        }
        //if (h==36) Rprintf("%s %i %f %f %f %f %f %f %f \n","here new: ",nmembers[36],Xi[36][0],Xi[36][1],Xi[h][0],Xi[h][1],Hbarh,Ybarh,Yvarh);
    }
    gsl_rng_free(r);
}

//Initialize Poisson regression coefs: intercepts are set equal to the mean SMR of the cluster;
//for empty clusters they are set randomly to mean SMR + error. Slopes are set to zero.
void initPoissonRegCoef(unsigned long int s, int *Y, double *E, int n, int nreg, int nres, int ncomp,
                        int *nmembers, int *compAlloc, double *smrbar, double beta[ncomp][nreg*nres]){
    double sumhN[nres];
    double sumhD[nres];
    int i, h, k;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    for (h = 0; h < ncomp; h++)
        for (k = 0; k < (nreg*nres); k++)
            beta[h][k] = 0.0;
    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0){
            for (k = 0; k < nres; k++){
                sumhN[k] = 0.0;
                sumhD[k] = 0.0;
            }
            for (i = 0; i < n; i++){
                if (compAlloc[i]==h){
                    for (k = 0; k < nres; k++){
                        sumhN[k] += (double) Y[i+k*n];
                        sumhD[k] += E[i+k*n];
                    }
                }
            }
            for (k = 0; k < nres; k++){
                if (sumhN[k] > 0)
                    beta[h][k*nreg] = log(sumhN[k]/sumhD[k]);
                else
                    beta[h][k*nreg] = log((sumhN[k]+1.0)/sumhD[k]);
            }
        }
        if (nmembers[h] == 0){
            for (k = 0; k < nres; k++) beta[h][k*nreg] = log(smrbar[k]) + gsl_ran_gaussian(r,1.0);
        }
    }
    gsl_rng_free(r);
}

//Initialize regression coefs for mixed responses: count, binomial, continuous. Intercepts are set equal to the mean
//of the cluster; for empty clusters they are set randomly to mean + error. Slopes are set to zero.
void initPoisBinomConRegCoef(unsigned long int s, double *Y, double *E, int n, int *cumnreg, int nres, int totNreg,
                             int ncomp, int *nmembers, int *compAlloc, double *ybar, double beta[ncomp][totNreg]){
    double sumhN[nres];
    double sumhD[nres];
    int i, h, k;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    for (h = 0; h < ncomp; h++)
        for (k = 0; k < totNreg; k++)
            beta[h][k] = 0.0;
    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0){
            for (k = 0; k < nres; k++){
                sumhN[k] = 0.0;
                sumhD[k] = 0.0;
            }
            for (i = 0; i < n; i++){
                if (compAlloc[i]==h){
                    for (k = 0; k < nres; k++){
                        sumhN[k] += Y[i+k*n];
                        sumhD[k] += E[i+k*n];
                    }
                }
            }
            if (sumhN[0] > 0.0)
                beta[h][cumnreg[0]] = log(sumhN[0]/sumhD[0]);
            else
                beta[h][cumnreg[0]] = log((sumhN[0]+1.0)/sumhD[0]);
            beta[h][cumnreg[1]] = logit((sumhN[1]+1.0)/(sumhD[1]+2.0));
            beta[h][cumnreg[2]] = sumhN[2]/sumhD[2];
        }
        if (nmembers[h] == 0){
            beta[h][cumnreg[0]] = log(ybar[0]) + gsl_ran_gaussian(r,1.0);
            beta[h][cumnreg[1]] = logit(ybar[1]) + gsl_ran_gaussian(r,1.0);
            beta[h][cumnreg[2]] = ybar[2] + gsl_ran_gaussian(r,1.0);
        }
    }
    gsl_rng_free(r);
}

//Initialize cluster means of the confounders (the clustering variables that is): for non-empty clusters
//mu_h = mean of the members, for empty clusters mu_h = overall mean + error.
void initMuh(unsigned long int s, double *W, int n, int nconf, int ncomp, double *wbar, double *wsd, int *compAlloc,
             int *nmembers, double muh[ncomp][nconf]){
    double sumh[nconf];
    int i, h, k;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0){
            for (k = 0; k < nconf; k++)
                sumh[k] = 0.0;
            for (i = 0; i < n; i++){
                if (compAlloc[i]==h){
                    for (k = 0; k < nconf; k++)
                        sumh[k] += W[k*n+i];
                }
            }
            for (k = 0; k < nconf; k++)
                muh[h][k] = sumh[k]/nmembers[h];
        }
        if (nmembers[h]==0){
            for (k = 0; k < nconf; k++)
                muh[h][k] = wbar[k] + gsl_ran_gaussian(r,wsd[k]);
        }
    }
    gsl_rng_free(r);
}

//Initialize cluster means of the confounders (the clustering variables that is): for non-empty clusters
//mu_h = mean of the members, for empty clusters mu_h = overall mean + error.
void initMuhFx(unsigned long int s, double *X, double *latentx, int n, int p, int ncomp, double *xbar, double *xsd,
               double ybar, double ysd, int *compAlloc, int *nmembers, double muh[ncomp][p]){
    double sumh[p];
    int i, h, k;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0){
            for (k = 0; k < p; k++)
                sumh[k] = 0.0;
            for (i = 0; i < n; i++){
                if (compAlloc[i]==h){
                    sumh[0] += latentx[i];
                    for (k = 0; k < p-1; k++)
                        sumh[k+1] += X[k*n+i];
                }
            }
            for (k = 0; k < p; k++)
                muh[h][k] = sumh[k]/nmembers[h];
        }
        if (nmembers[h]==0){
                muh[h][0] = ybar + gsl_ran_gaussian(r,ysd);
            for (k = 0; k < p-1; k++)
                muh[h][k+1] = xbar[k] + gsl_ran_gaussian(r,xsd[k]);
        }
    }
    gsl_rng_free(r);
}

//Returns c_{y-1}(E*gamma) (lower) and c_{y}(E*gamma) (upper) limits for 1 sampling unit, for nres count responses,
//given design matrix X, responses Y, offsets E, sample size, number of regressors, number of responses,
//sampling unit, regression coefficients, and two vectors where lower and upper limits are stored.
void calcLimits(double *X, int *Y, double *E, int n, int nreg, int nres, int i, double *beta, double *lower, double *upper){
    int k, j;
    double loggamma;
    for (k = 0; k < nres; k++){
        loggamma = 0;
        for (j = 0; j < nreg; j++)
            loggamma += beta[nreg*k+j] * X[j*n+i];
        if (Y[k*n+i]==0)
            lower[k]=-999.99;
        else
            lower[k] = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[k*n+i]-1,E[k*n+i]*exp(loggamma)));
        if (lower[k] < -999.99) lower[k] = -999.99;
        if (lower[k] > 999.99) lower[k] = 999.99;
        upper[k] = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[k*n+i],E[k*n+i]*exp(loggamma)));
        if (upper[k] < -999.99) upper[k] = -999.99;
        if (upper[k] > 999.99) upper[k] = 999.99;
    }
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) limits for ith sampling unit,
//given responses Y, offsets H, sampling unit, rates, and two doubles where lower and upper limits are stored.
void calcGLMLimits(int Y, double H, double *Xi, double *lower, double *upper, int family){
    double lmt = 9999.99;
    double normConst, saveCdf, mu, f, B; 
    saveCdf = 0;
    if (Y==0)
        *lower = -lmt;
    else{
        if (family == 1) saveCdf = gsl_cdf_poisson_P(Y-1,H*Xi[0]);            
        else if (family == 2) saveCdf = gsl_cdf_binomial_P(Y-1,Xi[0],(int) H);
        else if (family == 3) saveCdf = gsl_cdf_negative_binomial_P(Y-1,Xi[1]/(H+Xi[1]),Xi[0]);
        else if (family == 4) saveCdf = cdf_beta_binomial_P2(H,Y-1,Xi[0],Xi[1]);
        else if (family == 5) saveCdf = cdf_generalized_poisson_P3(Y-1,H*Xi[0],Xi[1],&normConst);//cdf_generalized_poisson_P3(Y-1,H*Xi[0],Xi[1]/sqrt(H),&normConst);
        else if (family == 6) saveCdf = cdf_com_poisson_P2(Y-1,H*Xi[0],Xi[1]);
        else if (family == 7) saveCdf = cdf_hyper_poisson_P2(Y-1,H*Xi[0],Xi[1]);
        else if (family == 8) saveCdf = cdf_tri_parametric_P2(Y-1,H*Xi[0],Xi[1],Xi[2]);
        *lower = gsl_cdf_ugaussian_Pinv(saveCdf);
    }
    if (*lower < -lmt) *lower = -lmt; 
    if (*lower > lmt) *lower = lmt;
    if (family == 1) saveCdf += gsl_ran_poisson_pdf(Y,H*Xi[0]);
    else if (family == 2) saveCdf += gsl_ran_binomial_pdf(Y,Xi[0],(int) H);
    else if (family == 3) saveCdf += gsl_ran_negative_binomial_pdf(Y,Xi[1]/(H+Xi[1]),Xi[0]);
    else if (family == 4) saveCdf += exp(gsl_sf_lnchoose(H,Y) + gsl_sf_lnbeta(Y+Xi[0],H-Y+Xi[1]) - gsl_sf_lnbeta(Xi[0],Xi[1]));
    else if (family == 5){		
		mu = H*Xi[0]; f = Xi[1]; B = -mu/(f-1);
		if ((f>=1)||(f < 1 && Y < B))
		    saveCdf += exp(log(mu)+(Y-1)*log(mu+(f-1)*Y)-Y*log(f)-(mu+(f-1)*Y)/f-gsl_sf_lnfact(Y))/normConst;		
	}
    else if (family == 6) *upper = gsl_cdf_ugaussian_Pinv(cdf_com_poisson_P2(Y,H*Xi[0],Xi[1]));
    else if (family == 7) *upper = gsl_cdf_ugaussian_Pinv(cdf_hyper_poisson_P2(Y,H*Xi[0],Xi[1]));
    else if (family == 8) *upper = gsl_cdf_ugaussian_Pinv(cdf_tri_parametric_P2(Y,H*Xi[0],Xi[1],Xi[2]));
    if (saveCdf > 1.0) saveCdf = 1.0;
    *upper = gsl_cdf_ugaussian_Pinv(saveCdf);
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) limits for ith sampling unit for both Y* and X*,
//given responses Y, offsets H, covariates X, sampling unit, rates, and two doubles where lower and upper limits are stored.
void calcGLMLimitsYX(int *Y, double *H, double *X, int i, int NDC, int n, double *Xi, double *lower, double *upper, int family){
    int j;
    double lmt = 9999.99;
    calcGLMLimits(Y[i],H[i],Xi,&lower[0],&upper[0],family);
/*    
    double saveCdf, mu, f, B, normConst;
    if (Y[i]==0)
        lower[0] = -lmt;
    else{
        if (family == 1) lower[0] = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[i]-1,H[i]*Xi[0]));
        else if (family == 2) lower[0] = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(Y[i]-1,Xi[0],(int) H[i]));
        else if (family == 3) lower[0] = gsl_cdf_ugaussian_Pinv(gsl_cdf_negative_binomial_P(Y[i]-1,Xi[1]/(H[i]+Xi[1]),Xi[0]));
        else if (family == 4) lower[0] = gsl_cdf_ugaussian_Pinv(cdf_beta_binomial_P(H[i],Y[i]-1,Xi[0],Xi[1]));
        else if (family == 5){
			 saveCdf = cdf_generalized_poisson_P3(Y[i]-1,H[i]*Xi[0],Xi[1],&normConst);
			 lower[0] = gsl_cdf_ugaussian_Pinv(saveCdf);
		 }
        else if (family == 6) lower[0] = gsl_cdf_ugaussian_Pinv(cdf_com_poisson_P2(Y[i]-1,H[i]*Xi[0],Xi[1]));
        else if (family == 7) lower[0] = gsl_cdf_ugaussian_Pinv(cdf_hyper_poisson_P2(Y[i]-1,H[i]*Xi[0],Xi[1]));
        else if (family == 8) lower[0] = gsl_cdf_ugaussian_Pinv(cdf_tri_parametric_P2(Y[i]-1,H[i]*Xi[0],Xi[1],Xi[2]));
    }
    if (lower[0] < -lmt) lower[0] = -lmt;
    if (lower[0] > lmt) lower[0] = lmt;
    if (family == 1) upper[0] = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[i],H[i]*Xi[0]));
    else if (family == 2) upper[0] = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(Y[i],Xi[0],(int) H[i]));
    else if (family == 3) upper[0] = gsl_cdf_ugaussian_Pinv(gsl_cdf_negative_binomial_P(Y[i],Xi[1]/(H[i]+Xi[1]),Xi[0]));
    else if (family == 4) upper[0] = gsl_cdf_ugaussian_Pinv(cdf_beta_binomial_P(H[i],Y[i],Xi[0],Xi[1]));
    else if (family == 5){ 
        mu = H[i]*Xi[0]; f = Xi[1]; B = -mu/(f-1);
        if ((f>=1)||(f < 1 && Y[i] < B))
		    saveCdf += exp(log(mu)+(Y[i]-1)*log(mu+(f-1)*Y[i])-Y[i]*log(f)-(mu+(f-1)*Y[i])/f-gsl_sf_lnfact(Y[i]))/normConst;
		if (saveCdf > 1.0) saveCdf = 1.0;
        upper[0] = gsl_cdf_ugaussian_Pinv(saveCdf);
	}
    else if (family == 6) upper[0] = gsl_cdf_ugaussian_Pinv(cdf_com_poisson_P2(Y[i],H[i]*Xi[0],Xi[1]));
    else if (family == 7) upper[0] = gsl_cdf_ugaussian_Pinv(cdf_hyper_poisson_P2(Y[i],H[i]*Xi[0],Xi[1]));
    else if (family == 8) upper[0] = gsl_cdf_ugaussian_Pinv(cdf_tri_parametric_P2(Y[i],H[i]*Xi[0],Xi[1],Xi[2]));
    if (upper[0] < -lmt) upper[0] = -lmt;
    if (upper[0] > lmt) upper[0] = lmt;
*/
    for (j = 0; j < NDC; j++){
        if (X[j*n+i]==0) {lower[j+1] = -lmt; upper[j+1] = 0;}
        else if (X[j*n+i]==1) {lower[j+1] = 0; upper[j+1] = lmt;}
	}
}

void calcGLMLimits2(int *Y, double *H, int i, double *Xi, double *lower, double *upper, int family, int *BFE){
    double normConst;
    double lmt = 9999.99;
    double lp = Xi[0] + Xi[1] * BFE[i];
    if (Y[i]==0)
        *lower = -lmt;
    else{
        if (family == 1) *lower = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[i]-1,H[i]*Xi[0]));
        else if (family == 2) *lower = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(Y[i]-1,Xi[0],(int) H[i]));
        else if (family == 3) *lower = gsl_cdf_ugaussian_Pinv(gsl_cdf_negative_binomial_P(Y[i]-1,Xi[1]/(H[i]+Xi[1]),Xi[0]));
        else if (family == 4) *lower = gsl_cdf_ugaussian_Pinv(cdf_beta_binomial_P(H[i],Y[i]-1,Xi[0],Xi[1]));
        else if (family == 5) *lower = gsl_cdf_ugaussian_Pinv(cdf_generalized_poisson_P3(Y[i]-1,H[i]*exp(lp),Xi[2],&normConst));
        else if (family == 6) *lower = gsl_cdf_ugaussian_Pinv(cdf_com_poisson_P(Y[i]-1,H[i]*Xi[0],Xi[1]));
        else if (family == 7) *lower = gsl_cdf_ugaussian_Pinv(cdf_hyper_poisson_P2(Y[i]-1,H[i]*Xi[0],Xi[1]));
        else if (family == 8) *lower = gsl_cdf_ugaussian_Pinv(cdf_tri_parametric_P(Y[i]-1,Xi[0],Xi[1],Xi[2]));
    }
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
    if (family == 1) *upper = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[i],H[i]*Xi[0]));
    else if (family == 2) *upper = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(Y[i],Xi[0],(int) H[i]));
    else if (family == 3) *upper = gsl_cdf_ugaussian_Pinv(gsl_cdf_negative_binomial_P(Y[i],Xi[1]/(H[i]+Xi[1]),Xi[0]));
    else if (family == 4) *upper = gsl_cdf_ugaussian_Pinv(cdf_beta_binomial_P(H[i],Y[i],Xi[0],Xi[1]));
    else if (family == 5) *upper = gsl_cdf_ugaussian_Pinv(cdf_generalized_poisson_P3(Y[i],H[i]*exp(lp),Xi[2],&normConst));
    else if (family == 6) *upper = gsl_cdf_ugaussian_Pinv(cdf_com_poisson_P(Y[i],H[i]*Xi[0],Xi[1]));
    else if (family == 7) *upper = gsl_cdf_ugaussian_Pinv(cdf_hyper_poisson_P2(Y[i],H[i]*Xi[0],Xi[1]));
    else if (family == 8) *upper = gsl_cdf_ugaussian_Pinv(cdf_tri_parametric_P(Y[i],Xi[0],Xi[1],Xi[2]));
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) limits for ith predictive scenario,
//given offsets H, value, scenario, rates, and two doubles where lower and upper limits are stored.
void calcGLMLimitsPred(double *H, int k, int i, double *Xi, double *lower, double *upper, int family){
    double lmt = 9999.99;
    if (k==0)
        *lower = -lmt;
    else{
        if (family == 1) *lower = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(k-1,H[i]*Xi[0]));
        else if (family == 2) *lower = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(k-1,Xi[0],(int) H[i]));
        else if (family == 3) *lower = gsl_cdf_ugaussian_Pinv(gsl_cdf_negative_binomial_P(k-1,Xi[1]/(H[i]+Xi[1]),Xi[0]));
        else if (family == 4) *lower = gsl_cdf_ugaussian_Pinv(cdf_beta_binomial_P(H[i],k-1,Xi[0],Xi[1]));
    }
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
    if (family == 1) *upper = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(k,H[i]*Xi[0]));
    else if (family == 2) *upper = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(k,Xi[0],(int) H[i]));
    else if (family == 3) *upper = gsl_cdf_ugaussian_Pinv(gsl_cdf_negative_binomial_P(k,Xi[1]/(H[i]+Xi[1]),Xi[0]));
    else if (family == 4) *upper = gsl_cdf_ugaussian_Pinv(cdf_beta_binomial_P(H[i],k,Xi[0],Xi[1]));
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) for Generalized Poisson
void calcGLMLimitsPredGP(int k, double H, double *Xi, double *lower, double *upper,
                         double *CDFL, double *CDFU, double *normConst){
    double lmt = 9999.99;
    double mu = H*Xi[0];
    double f = Xi[1]; 
    double B = -mu/(f-1);
    *normConst = 1.0;
    *CDFL = 0.0;
    *CDFU = 0.0;
    if (k > 0){
        *CDFL = cdf_generalized_poisson_P3(k-1,mu,f,normConst);                
        *CDFU = *CDFL;
        if ((f>=1)||(f < 1 && k < B))
		    *CDFU += exp(log(mu)+(k-1)*log(mu+(f-1)*k)-k*log(f)-(mu+(f-1)*k)/f-gsl_sf_lnfact(k)) / *normConst;
    }                   
    else if (k==0)
        *CDFU = cdf_generalized_poisson_P3(k,mu,f,normConst);    
    if (*CDFL > 1.0) *CDFL = 1.0;
    if (*CDFU > 1.0) *CDFU = 1.0;
    *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
    if (k==0) *lower = -lmt;
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) for Generalized Poisson with fixed effect
void calcGLMLimitsPredGP2(double *H, int q, int i, double *Xi, double *lower, double *upper,
                          double *CDFL, double *CDFU, double *normConst,int *BFE){
    double lmt = 9999.99;
    *normConst = 1.0;
    *CDFL = 0.0;
    *CDFU = 0.0;
    int j = 0;
    double mu = H[i]*exp(Xi[0] + Xi[1] * BFE[i]);
    double f = Xi[2];
    double B = -mu/(f-1);

    if (f == 1.0){
        if (q > 0) *CDFL = gsl_cdf_poisson_P(q-1,mu);
        *CDFU = gsl_cdf_poisson_P(q,mu);
    }
    else if (f > 1.0){
        for (j = 0; j < q; j++)
            *CDFL += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
        *CDFU = *CDFL + exp(log(mu)+(q-1)*log(mu+(f-1)*q)-q*log(f)-(mu+(f-1)*q)/f-gsl_sf_lnfact(q));
    }
    else if ((mu < 1.0 && f < 0.97) || (mu < 2.0 && f < 0.80) || (mu < 3.0 && f < 0.65) || (mu < 4.0 && f < 0.60) ||
            (mu < 5.0 && f < 0.55))
    {
        while (j < q && j < B){
            if ((mu+(f-1)*j)>0) *CDFL += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
            j++;
        }
        *CDFU = *CDFL;
        if (q < B) if ((mu+(f-1)*q)>0) *CDFU += exp(log(mu)+(q-1)*log(mu+(f-1)*q)-q*log(f)-(mu+(f-1)*q)/f-gsl_sf_lnfact(q));
        *normConst = *CDFU;
        while (j < B){
            j++;
            if ((mu+(f-1)*j)>0) *normConst += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
        }
        *CDFL = *CDFL/(*normConst);
        *CDFU = *CDFU/(*normConst);
    }
    else
    {
        while (j < q && j < B){
            if ((mu+(f-1)*j)>0) *CDFL += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
            j++;
        }
        *CDFU = *CDFL;
        if (q < B) if ((mu+(f-1)*q)>0) *CDFU += exp(log(mu)+(q-1)*log(mu+(f-1)*q)-q*log(f)-(mu+(f-1)*q)/f-gsl_sf_lnfact(q));
    }
    if (*CDFL > 1.0) *CDFL = 1.0;
    if (*CDFU > 1.0) *CDFU = 1.0;
    *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
    if (q==0) *lower = -lmt;
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) for C Poisson
void calcGLMLimitsPredCP(double *H, int x, int i, double *Xi, double *lower, double *upper,
                         double *CDFL, double *CDFU, double *normConst){
    double lmt = 9999.99;
    *CDFU = 0.0;

    double lambda = pow(H[i]*Xi[0] + (Xi[1]-1)/(2*Xi[1]),Xi[1]); //Xi[0];
    double nu = Xi[1];

    int j, k, K; //of (32) of COM-Poisson revival paper
    int MIN; //min(x,k)
    double epsilon = 0.99;//of (32) of COM-Poisson revival paper
    double R;
    K = 0;
    while(lambda/pow(K+1,nu) > epsilon) K++;
    k = K + 2;
    MIN = x;
    if (k < MIN) MIN = k;
    for (j = 0; j < MIN+1; j++)
        *CDFU += exp(j*log(lambda)-nu*gsl_sf_lnfact(j));
    *normConst = *CDFU;
    if (x > MIN)
        for (j = MIN+1; j < x+1; j++)
            *CDFU += exp(j*log(lambda)-nu*gsl_sf_lnfact(j));
    *CDFL = *CDFU - exp(x*log(lambda)-nu*gsl_sf_lnfact(x));
    if (k > MIN)
        for (j = MIN+1; j < k+1; j++)
            *normConst += exp(j*log(lambda)-nu*gsl_sf_lnfact(j));
    R = exp((k+1)*log(lambda)-nu*gsl_sf_lnfact(k+1)-log(1-epsilon)-log(*normConst));
    while(R > 0.000001){
        k++;
        *normConst += exp(k*log(lambda)-nu*gsl_sf_lnfact(k));
        R = exp((k+1)*log(lambda)-nu*gsl_sf_lnfact(k+1)-log(1-epsilon)-log(*normConst));
    }
    *CDFL = *CDFL/(*normConst);
    *CDFU = *CDFU/(*normConst);
    if (*CDFL > 1.0) *CDFL = 1.0;
    if (*CDFU > 1.0) *CDFU = 1.0;
    *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
    if (x==0) *lower = -lmt;
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) for hyper Poisson
void calcGLMLimitsPredHP(double *H, int x, int i, double *Xi, double *lower, double *upper, double *CDFL, double *CDFU, double *normConst){
    double lmt = 9999.99;
    *CDFL = 0.0;
    double lambda = solve_hyper_poisson(H[i]*Xi[0],Xi[1]); //H[i]*Xi[0];
    double gamma = Xi[1];
    int j;
    double C = 1/gsl_sf_hyperg_1F1(1.0,gamma,lambda);
    *normConst = C;
    for (j = 0; j < x; j++) *CDFL += C * exp(j*log(lambda) - gsl_sf_lnpoch(gamma,j));
    *CDFU = *CDFL + C * exp(x*log(lambda) - gsl_sf_lnpoch(gamma,x));
    if (*CDFL > 1.0) *CDFL = 1.0;
    if (*CDFU > 1.0) *CDFU = 1.0;
    *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
    if (x==0) *lower = -lmt;
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) for CTP
void calcGLMLimitsPredCTP(double *H, int x, int i, double *Xi, double *lower, double *upper, double *CDFL, double *CDFU, double *normConst){
    double lmt = 9999.99;
    *CDFL = 0.0;
    double beta = sqrt(H[i]*Xi[0]*(Xi[1]-2*Xi[2]-1) - Xi[2]*Xi[2]);
    double gamma = Xi[1];
    double alpha = Xi[2];
    int j;
    gsl_sf_result lnr1, lnr2, lnr3, lnr4, lnr5, lnr6, arg1, arg2, arg3, arg4, arg5, arg6;
    gsl_sf_lngamma_complex_e(gamma-alpha,-beta,&lnr1,&arg1);
    gsl_sf_lngamma_complex_e(gamma-alpha,beta,&lnr2,&arg2);
    gsl_sf_lngamma_complex_e(alpha,-beta,&lnr5,&arg5);
    gsl_sf_lngamma_complex_e(alpha,beta,&lnr6,&arg6);
    double f0 = lnr1.val + lnr2.val - gsl_sf_lngamma(gamma) - gsl_sf_lngamma(gamma-2*alpha);
    f0 = exp(f0) * gsl_sf_gamma(gamma);
    *normConst = f0;
    for (j = 0; j < x; j++){
        gsl_sf_lngamma_complex_e(alpha+j,-beta,&lnr3,&arg3);
        gsl_sf_lngamma_complex_e(alpha+j,beta,&lnr4,&arg4);
        *CDFL += exp(lnr3.val + lnr4.val - lnr5.val - lnr6.val - gsl_sf_lngamma(gamma+j) - gsl_sf_lnfact(j));
    }
    *CDFL *= f0;
    gsl_sf_lngamma_complex_e(alpha+x,-beta,&lnr3,&arg3);
    gsl_sf_lngamma_complex_e(alpha+x,beta,&lnr4,&arg4);
    *CDFU = *CDFL + f0 * exp(lnr3.val + lnr4.val - lnr5.val - lnr6.val - gsl_sf_lngamma(gamma+x) - gsl_sf_lnfact(x));
    if (*CDFL > 1.0) *CDFL = 1.0;
    if (*CDFU > 1.0) *CDFU = 1.0;
    *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
    if (x==0) *lower = -lmt;
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) limits for ith predictive scenario,
//given offsets H, value, scenario, rates, and two doubles where lower and upper limits are stored.
void calcGLMLimitsPredL(int k, double H, double *Xi, double *lower, int family){
    double normConst;
    double lmt = 9999.99;
    if (k==0)
        *lower = -lmt;
    else{
        if (family == 1) *lower = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(k-1,H*Xi[0]));
        else if (family == 2) *lower = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(k-1,Xi[0],(int) H));
        else if (family == 3) *lower = gsl_cdf_ugaussian_Pinv(gsl_cdf_negative_binomial_P(k-1,Xi[1]/(H+Xi[1]),Xi[0]));
        else if (family == 4) *lower = gsl_cdf_ugaussian_Pinv(cdf_beta_binomial_P2(H,k-1,Xi[0],Xi[1]));
        else if (family == 5) *lower = gsl_cdf_ugaussian_Pinv(cdf_generalized_poisson_P3(k-1,H*Xi[0],Xi[1],&normConst));
    }
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
}

//For GP only
void calcGLMLimitsPredLGP(int k, double H, double *Xi, double *lower, double *CDFL, double normConst){
    double lmt = 9999.99;
    double mu = H*Xi[0];
    double f = Xi[1];
    double B = -mu/(f-1);
    double prob = 0.0;
    if (k==0)
        *lower = -lmt;
    else{
        if (f < 1 && k < B && (mu+(f-1)*k) > 0) prob = exp(log(mu)+(k-1)*log(mu+(f-1)*k)-k*log(f)-(mu+(f-1)*k)/f-gsl_sf_lnfact(k))/normConst;
        if (f >= 1) prob = exp(log(mu)+(k-1)*log(mu+(f-1)*k)-k*log(f)-(mu+(f-1)*k)/f-gsl_sf_lnfact(k))/normConst;
        *CDFL = *CDFL - prob;
        if (*CDFL < 0.0) *CDFL = 0.0;
        *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    }
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
}

//For GP only
void calcGLMLimitsPredLGP2(double *H, int k, int i, double *Xi, double *lower, double *CDFL, double normConst, int *BFE){
    double lmt = 9999.99;
    double mu = H[i]*exp(Xi[0] + Xi[1] * BFE[i]);
    double f = Xi[2];
    double B = -mu/(f-1);
    double prob = 0.0;
    if (k==0)
        *lower = -lmt;
    else{
        if (f < 1 && k < B && (mu+(f-1)*k) > 0) prob = exp(log(mu)+(k-1)*log(mu+(f-1)*k)-k*log(f)-(mu+(f-1)*k)/f-gsl_sf_lnfact(k))/normConst;
        if (f >= 1) prob = exp(log(mu)+(k-1)*log(mu+(f-1)*k)-k*log(f)-(mu+(f-1)*k)/f-gsl_sf_lnfact(k))/normConst;
        *CDFL = *CDFL - prob;
        if (*CDFL < 0.0) *CDFL = 0.0;
        *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    }
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
}

//For CP only
void calcGLMLimitsPredLCP(double *H, int k, int i, double *Xi, double *lower, double *CDFL, double normConst){
    double lmt = 9999.99;
    double lambda = pow(H[i]*Xi[0] + (Xi[1]-1)/(2*Xi[1]),Xi[1]); //Xi[0];
    double nu = Xi[1];
    double prob = 0.0;
    if (k==0)
        *lower = -lmt;
    else{
        prob = exp(k*log(lambda)-nu*gsl_sf_lnfact(k))/normConst;
        *CDFL = *CDFL - prob;
        if (*CDFL < 0.0) *CDFL = 0.0;
        *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    }
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
}

//For HP only
void calcGLMLimitsPredLHP(double *H, int k, int i, double *Xi, double *lower, double *CDFL, double normConst){
    double lmt = 9999.99;
    double lambda = solve_hyper_poisson(H[i]*Xi[0],Xi[1]); //H[i]*Xi[0];
    double gamma = Xi[1];
    //double lambda = H[i]*Xi[0];
    //double gamma = Xi[1];
    double prob;
    if (k==0)
        *lower = -lmt;
    else{
        prob = exp(k*log(lambda) - gsl_sf_lngamma(gamma+k) + gsl_sf_lngamma(gamma))/gsl_sf_hyperg_1F1(1.0,gamma,lambda);
        *CDFL = *CDFL - prob;
        if (*CDFL < 0.0) *CDFL = 0.0;
        *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    }
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
}

//For CTP only
void calcGLMLimitsPredLCTP(double *H, int k, int i, double *Xi, double *lower, double *CDFL, double normConst){
    double lmt = 9999.99;
    //double beta = Xi[0];
    double beta = sqrt(H[i]*Xi[0]*(Xi[1]-2*Xi[2]-1) - Xi[2]*Xi[2]);
    double gamma = Xi[1];
    double alpha = Xi[2];
    gsl_sf_result lnr1, lnr2, lnr3, lnr4, lnr5, lnr6, arg1, arg2, arg3, arg4, arg5, arg6;
    double f0;
    double prob;
    if (k==0)
        *lower = -lmt;
    else{
        gsl_sf_lngamma_complex_e(gamma-alpha,-beta,&lnr1,&arg1);
        gsl_sf_lngamma_complex_e(gamma-alpha,beta,&lnr2,&arg2);
        gsl_sf_lngamma_complex_e(alpha+k,-beta,&lnr3,&arg3);
        gsl_sf_lngamma_complex_e(alpha+k,beta,&lnr4,&arg4);
        gsl_sf_lngamma_complex_e(alpha,-beta,&lnr5,&arg5);
        gsl_sf_lngamma_complex_e(alpha,beta,&lnr6,&arg6);
        f0 = lnr1.val + lnr2.val - gsl_sf_lngamma(gamma) - gsl_sf_lngamma(gamma-2*alpha);
        f0 = exp(f0) * gsl_sf_gamma(gamma);
        prob = exp(lnr3.val + lnr4.val - lnr5.val - lnr6.val - gsl_sf_lngamma(gamma+k) - gsl_sf_lnfact(k));
        *CDFL = *CDFL - f0 * prob;
        if (*CDFL < 0.0) *CDFL = 0.0;
        *lower = gsl_cdf_ugaussian_Pinv(*CDFL);
    }
    if (*lower < -lmt) *lower = -lmt;
    if (*lower > lmt) *lower = lmt;
}

//Returns c_{y-1}(gamma,H) (lower) and c_{y}(gamma,H) (upper) limits for ith predictive scenario,
//given offsets H, value, scenario, rates, and two doubles where lower and upper limits are stored.
void calcGLMLimitsPredU(int k, double H, double *Xi, double *upper, int family){
        double normConst;
        double lmt = 9999.99;
        if (family == 1) *upper = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(k,H*Xi[0]));
        else if (family == 2) *upper = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(k,Xi[0],(int) H));
        else if (family == 3) *upper = gsl_cdf_ugaussian_Pinv(gsl_cdf_negative_binomial_P(k,Xi[1]/(H+Xi[1]),Xi[0]));
        else if (family == 4) *upper = gsl_cdf_ugaussian_Pinv(cdf_beta_binomial_P2(H,k,Xi[0],Xi[1]));
        else if (family == 5) *upper = gsl_cdf_ugaussian_Pinv(cdf_generalized_poisson_P3(k,H*Xi[0],Xi[1],&normConst));
        if (*upper < -lmt) *upper = -lmt;
        if (*upper > lmt) *upper = lmt;
}

//For GP only
void calcGLMLimitsPredUGP(int k, double H, double *Xi, double *upper, double *CDFU, double normConst){
        double lmt = 9999.99;
        double mu = H*Xi[0];
        double f = Xi[1];
        double B = -mu/(f-1);
        double prob = 0.0;
        if (f < 1 && k < B && (mu+(f-1)*k) > 0) prob = exp(log(mu)+(k-1)*log(mu+(f-1)*k)-k*log(f)-(mu+(f-1)*k)/f-gsl_sf_lnfact(k))/normConst;
        if (f >= 1) prob = exp(log(mu)+(k-1)*log(mu+(f-1)*k)-k*log(f)-(mu+(f-1)*k)/f-gsl_sf_lnfact(k))/normConst;
        *CDFU = *CDFU + prob;
        if (*CDFU > 1.0) *CDFU = 1.0;
        *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
        if (*upper < -lmt) *upper = -lmt;
        if (*upper > lmt) *upper = lmt;
}

//For GP only
void calcGLMLimitsPredUGP2(double *H, int k, int i, double *Xi, double *upper, double *CDFU, double normConst, int *BFE){
        double lmt = 9999.99;
        double mu = H[i]*exp(Xi[0] + Xi[1] * BFE[i]);
        double f = Xi[2];
        double B = -mu/(f-1);
        double prob = 0.0;
        if (f < 1 && k < B && (mu+(f-1)*k) > 0) prob = exp(log(mu)+(k-1)*log(mu+(f-1)*k)-k*log(f)-(mu+(f-1)*k)/f-gsl_sf_lnfact(k))/normConst;
        if (f >= 1) prob = exp(log(mu)+(k-1)*log(mu+(f-1)*k)-k*log(f)-(mu+(f-1)*k)/f-gsl_sf_lnfact(k))/normConst;
        *CDFU = *CDFU + prob;
        if (*CDFU > 1.0) *CDFU = 1.0;
        *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
        if (*upper < -lmt) *upper = -lmt;
        if (*upper > lmt) *upper = lmt;
}

//For CP only
void calcGLMLimitsPredUCP(double *H, int k, int i, double *Xi, double *upper, double *CDFU, double normConst){
    double lmt = 9999.99;
    double lambda = pow(H[i]*Xi[0] + (Xi[1]-1)/(2*Xi[1]),Xi[1]); //Xi[0];
    double nu = Xi[1];
    double prob;
    prob = exp(k*log(lambda)-nu*gsl_sf_lnfact(k))/normConst;
    *CDFU = *CDFU + prob;
    if (*CDFU > 1.0) *CDFU = 1.0;
    *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//For HP only
void calcGLMLimitsPredUHP(double *H, int k, int i, double *Xi, double *upper, double *CDFU, double normConst){
    double lmt = 9999.99;
    double lambda = solve_hyper_poisson(H[i]*Xi[0],Xi[1]); //H[i]*Xi[0];
    double gamma = Xi[1];
    //double lambda = H[i]*Xi[0];
    //double gamma = Xi[1];
    double prob;
    prob = exp(k*log(lambda) - gsl_sf_lngamma(gamma+k) + gsl_sf_lngamma(gamma))/gsl_sf_hyperg_1F1(1.0,gamma,lambda);
    *CDFU = *CDFU + prob;
    if (*CDFU > 1.0) *CDFU = 1.0;
    *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//For CTP only
void calcGLMLimitsPredUCTP(double *H, int k, int i, double *Xi, double *upper, double *CDFU, double normConst){
    double lmt = 9999.99;
    //double beta = Xi[0];
    double beta = sqrt(H[i]*Xi[0]*(Xi[1]-2*Xi[2]-1) - Xi[2]*Xi[2]);
    double gamma = Xi[1];
    double alpha = Xi[2];
    gsl_sf_result lnr1, lnr2, lnr3, lnr4, lnr5, lnr6, arg1, arg2, arg3, arg4, arg5, arg6;
    double f0;
    double prob;
    gsl_sf_lngamma_complex_e(gamma-alpha,-beta,&lnr1,&arg1);
    gsl_sf_lngamma_complex_e(gamma-alpha,beta,&lnr2,&arg2);
    gsl_sf_lngamma_complex_e(alpha+k,-beta,&lnr3,&arg3);
    gsl_sf_lngamma_complex_e(alpha+k,beta,&lnr4,&arg4);
    gsl_sf_lngamma_complex_e(alpha,-beta,&lnr5,&arg5);
    gsl_sf_lngamma_complex_e(alpha,beta,&lnr6,&arg6);
    f0 = lnr1.val + lnr2.val - gsl_sf_lngamma(gamma) - gsl_sf_lngamma(gamma-2*alpha);
    f0 = exp(f0) * gsl_sf_gamma(gamma);
    prob = exp(lnr3.val + lnr4.val - lnr5.val - lnr6.val - gsl_sf_lngamma(gamma+k) - gsl_sf_lnfact(k));
    *CDFU = *CDFU + f0 * prob;
    if (*CDFU > 1.0) *CDFU = 1.0;
    *upper = gsl_cdf_ugaussian_Pinv(*CDFU);
    if (*upper < -lmt) *upper = -lmt;
    if (*upper > lmt) *upper = lmt;
}

//Returns c_{y-1} (lower) and c_{y} (upper) limits for 1 sampling unit, for one count and one binomial response,
//given design matrix X, responses Y, offsets E, sample size, vector with number of regressors, cumulative number of regressors,
//number of discrete responses, sampling unit, regression coefficients, and two vectors where lower and upper limits are stored.
void calcMixedLimits(double *X, double *Y, double *E, int n, int *nreg, int *cumnreg, int nDres, int i, double *beta, double *lower, double *upper){
    int k, j;
    double lp;
    for (k = 0; k < nDres; k++){
        lp = 0;
        for (j = cumnreg[k]; j < cumnreg[k+1]; j++)
            lp += beta[j] * X[j*n+i];
        if (Y[k*n+i]==0)
            lower[k]=-999.99;
        else{
            if (k == 0) lower[k] = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[k*n+i]-1,E[k*n+i]*exp(lp)));
            if (k == 1) lower[k] = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(Y[k*n+i]-1,invlogit(lp),E[k*n+i]));
        }
        if (lower[k] < -999.99) lower[k] = -999.99;
        if (lower[k] > 999.99) lower[k] = 999.99;
        if (k == 0) upper[k] = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[k*n+i],E[k*n+i]*exp(lp)));
        if (k == 1) upper[k] = gsl_cdf_ugaussian_Pinv(gsl_cdf_binomial_P(Y[k*n+i],invlogit(lp),E[k*n+i]));
        if (upper[k] < -999.99) upper[k] = -999.99;
        if (upper[k] > 999.99) upper[k] = 999.99;
    }
}

//Returns c_{y-1}(E*gamma) (lower) and c_{y}(E*gamma) (upper) limits for n_h sampling units, for nres responses,
//given design matrix X that does Not contain a column of ones, responses Y, offsets E, sample size, number of regressors,
//number of responses, number of components, component of interest, component allocation, regression coefficients,
//and two matrices where lower and upper limits are stored.
void calcLimitsXM1(double *X, int *Y, double *E, int n, int nreg, int nres, int i, double *beta, double *lower, double *upper){
    int k, j;
    double loggamma;
    for (k = 0; k < nres; k++){
        loggamma = beta[(nreg+1)*k];
        for (j = 0; j < nreg; j++)
            loggamma += beta[(nreg+1)*k+1+j] * X[j*n+i];
        if (Y[k*n+i]==0)
            lower[k]=-999.99;
        else
            lower[k] = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[k*n+i]-1,E[k*n+i]*exp(loggamma)));
        if (lower[k] < -999.99) lower[k] = -999.99;
        if (lower[k] > 999.99) lower[k] = 999.99;
        upper[k] = gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(Y[k*n+i],E[k*n+i]*exp(loggamma)));
        if (upper[k] < -999.99) upper[k] = -999.99;
        if (upper[k] > 999.99) upper[k] = 999.99;
    }
}


//Calculates the covariance matrix of cluster h given vectors of latent variables of length nres of mean 0,
//and vectors of confounding variables W of length nconf and mean muh
void setSh(double *W, int n, int nres, int nconf, int ncomp, int h, int *nmembers, int *compAlloc, double Ystar[n][nres],
           double muh[ncomp][nconf], gsl_matrix *Sh){
    gsl_matrix_set_zero(Sh);
    if (nmembers[h] > 0){
        int i, k;
        double baseSh[nres+nconf];
        gsl_matrix_view Storeni;
        for (i = 0; i < n; i++){
            if (compAlloc[i]==h){
                for (k = 0; k < nres; k++)
                    baseSh[k] = Ystar[i][k];
                for (k = 0; k < nconf; k++)
                    baseSh[nres+k] = W[k*n+i] - muh[h][k];
                Storeni = gsl_matrix_view_array(baseSh,nres+nconf,1);
                gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,&Storeni.matrix,&Storeni.matrix,1.0,Sh);
            }
        }
    }
}

//Calculates the covariance matrix of cluster h given vectors of latent variables of length nDres of mean 0,
//vectors of continuous responses of length nCres of mean x^T beta and
//vectors of confounding variables W of length nconf and mean muh
void setShMixed(int n, int nDres, int nCres, int nres, int nconf, int ncomp, int totNreg, int *cumnreg, int h, int *nmembers,
                int *compAlloc, double Ystar[n][nDres], double *Y, double *X, double beta[ncomp][totNreg], double *W,
                double muh[ncomp][nconf], gsl_matrix *Sh){
    gsl_matrix_set_zero(Sh);
    if (nmembers[h] > 0){
        int i, j, k;
        double lp;
        double baseSh[nres+nconf];
        gsl_matrix_view Storeni;
        for (i = 0; i < n; i++){
            if (compAlloc[i]==h){
                for (k = 0; k < nDres; k++)
                    baseSh[k] = Ystar[i][k];
                for (k = 0; k < nCres; k++){
                    lp = 0.0;
                    for (j = cumnreg[nDres+k]; j < cumnreg[nDres+k+1]; j++)
                        lp += beta[h][j] * X[j*n+i];
                    baseSh[nDres+k] = Y[(nDres+k)*n+i] - lp;
                }
                for (k = 0; k < nconf; k++)
                    baseSh[nres+k] = W[k*n+i] - muh[h][k];
                Storeni = gsl_matrix_view_array(baseSh,nres+nconf,1);
                gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,&Storeni.matrix,&Storeni.matrix,1.0,Sh);
            }
        }
    }
}

//Calculates the covariance matrix of cluster h given vectors of latent variables of length nDres of mean 0,
//vectors of continuous responses of length nCres of mean x^T beta and
void setShMixedFG(int n, int nDres, int nCres, int nres, int ncomp, int totNreg, int *cumnreg, int h, int *nmembers,
                int *compAlloc, double Ystar[n][nDres], double *Y, double *X, double beta[ncomp][totNreg], gsl_matrix *Sh){
    gsl_matrix_set_zero(Sh);
    if (nmembers[h] > 0){
        int i, j, k;
        double lp;
        double baseSh[nres];
        gsl_matrix_view Storeni;
        for (i = 0; i < n; i++){
            if (compAlloc[i]==h){
                for (k = 0; k < nDres; k++)
                    baseSh[k] = Ystar[i][k];
                for (k = 0; k < nCres; k++){
                    lp = 0.0;
                    for (j = cumnreg[nDres+k]; j < cumnreg[nDres+k+1]; j++)
                        lp += beta[h][j] * X[j*n+i];
                    baseSh[nDres+k] = Y[(nDres+k)*n+i] - lp;
                }
                Storeni = gsl_matrix_view_array(baseSh,nres,1);
                gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,&Storeni.matrix,&Storeni.matrix,1.0,Sh);
            }
        }
    }
}

//Initializes matrices E_h, Sigma^*_h, and D_h to identity
void initEDS(int ncomp, int totran, double SigmaSh[ncomp][totran*totran], double Eh[ncomp][totran*totran],
             double Dh[ncomp][totran*totran]){
    int h, k;
    for (h = 0; h < ncomp; h++){
        for (k = 0; k < (totran*totran); k++){
            SigmaSh[h][k] = 0.0;
            Eh[h][k] = 0.0;
            Dh[h][k] = 0.0;
        }
    }
    for (h = 0; h < ncomp; h++){
        for (k = 0; k < totran; k++){
            SigmaSh[h][(totran+1)*k] = 1.0;
            Eh[h][(totran+1)*k] = 1.0;
            Dh[h][(totran+1)*k] = 1.0;
        }
    }
}

//Calculates the cluster totals of the latent and observed confounder variables: inputs are
//confounder values, sample size, number of responses, number of confounders, cluster label,
//component allocaions, imputed latent variables and the variable that contains the totals
void calcTotals(double *W, int n, int nres, int nconf, int h, int *nmembers, int *compAlloc, double Ystar[n][nres], double *totalh){
    int i, k;
    int totran = nres + nconf;
    for (k = 0; k < totran; k++)
        totalh[k] = 0.0;
    if (nmembers[h] > 0){
        for (i = 0; i < n; i++){
            if (compAlloc[i]==h){
                for (k = 0; k < nres; k++)
                    totalh[k] += Ystar[i][k];
                for (k = 0; k < nconf; k++)
                    totalh[nres+k] += W[k*n+i];
            }
        }
    }
}

//Calculates the following 2 cluster totals: sum X^T S^{-1} X, sum X^T S^{-1} v. Inputs are
//number of discrete responses, number of continuous responses, number of confounders, number of regressors, sample size,
//component, nmembers, component allocation, latent variables, observed Y, obs confounders, design matrix, store sum 1,
//store sum 2, Sigma inverse, and X_i^{star}
void calcMixedTotals(int nDres, int nCres, int nconf, int totNregC, int *nreg, int *cumnreg, int n, int h, int *nmembers,
                     int * compAlloc, double Ystar[n][nDres], double *Y, double *W, double *X, gsl_vector *TOT1,
                     gsl_matrix *TOT2, gsl_matrix *SigmaInv, gsl_matrix *Xistar){
    int i, j, k, dim;
    dim = nDres+nCres+nconf;
    gsl_matrix *XS = gsl_matrix_alloc(totNregC+nconf,dim);
    gsl_vector *v = gsl_vector_alloc(dim);
    gsl_vector_set_zero(TOT1);
    gsl_matrix_set_zero(TOT2);
    if (nmembers[h] > 0){
        for (i = 0; i < n; i++){
            if (compAlloc[i]==h){
                for (k = 0; k < nCres; k++)
                    for (j = 0; j < nreg[nDres+k]; j++)
                        gsl_matrix_set(Xistar,nDres+k,cumnreg[nDres+k]-cumnreg[nDres]+j,X[n*(cumnreg[nDres+k]+j)+i]);
                for (k = 0; k < nDres; k++)
                    gsl_vector_set(v,k,Ystar[i][k]);
                for (k = 0; k < nCres; k++)
                    gsl_vector_set(v,nDres+k,Y[(nDres+k)*n+i]);
                for (k = 0; k < nconf; k++)
                    gsl_vector_set(v,nDres+nCres+k,W[k*n+i]);
                gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Xistar,SigmaInv,0.0,XS);
                gsl_blas_dgemv(CblasNoTrans,1.0,XS,v,1.0,TOT1);
                gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,XS,Xistar,1.0,TOT2);
            }
        }
    }
    gsl_matrix_free(XS); gsl_vector_free(v);
}

void labelSwitchingA(unsigned long int s, int n, int nconf, int totNreg, int nres, int totran, int ncomp,
                     int *nmembers, double muh[ncomp][nconf], double SigmaSh[ncomp][totran*totran],
                     double Dh[ncomp][totran*totran], double Eh[ncomp][totran*totran],
                     double beta[ncomp][totNreg], int *compAlloc, double *pi){
    int j, h, komp, labelA, labelB;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double nonempty[ncomp];
    unsigned int vecAlloc[ncomp];
    double switchA, decisionA;

    double muhtemp[nconf];
    double vartemp[totran*totran];
    double betatemp[totNreg];
    int compAlloctemp[n];
    int nmemberstemp;

    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0)
            nonempty[h] = 1.0;
        else
            nonempty[h] = 0.0;
    }

    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);

    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelA = komp;

    nonempty[labelA] = 0.0;
    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);
    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelB = komp;

    if ((labelA >= ncomp) ||  (labelB >= ncomp))
        switchA = 0.0;
    else
        switchA = pow((pi[labelB]/pi[labelA]), (double)  (nmembers[labelA]-nmembers[labelB]));

    if (switchA > 1.0) switchA = 1.0;

    decisionA = gsl_ran_flat(r,0.0,1.0);

    if (switchA > decisionA){
        for (j = 0; j < nconf; j++) muhtemp[j] = muh[labelA][j];
        for (j = 0; j < nconf; j++) muh[labelA][j] = muh[labelB][j];
        for (j = 0; j < nconf; j++) muh[labelB][j] = muhtemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = SigmaSh[labelA][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelA][j] = SigmaSh[labelB][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelB][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Dh[labelA][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelA][j] = Dh[labelB][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelB][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Eh[labelA][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelA][j] = Eh[labelB][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelB][j] = vartemp[j];

        for (j = 0; j < totNreg; j++) betatemp[j] = beta[labelA][j];
        for (j = 0; j < totNreg; j++) beta[labelA][j] = beta[labelB][j];
        for (j = 0; j < totNreg; j++) beta[labelB][j] = betatemp[j];

        nmemberstemp = nmembers[labelA];
        nmembers[labelA] = nmembers[labelB];
        nmembers[labelB] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelA) compAlloc[j]=labelB;
            if (compAlloctemp[j] == labelB) compAlloc[j]=labelA;
        }
    }
    gsl_rng_free(r);
}

void SpatialLabelSwitchingA(unsigned long int s, int n, int nconf, int totNreg, int nres, int totran, int ncomp, int *nmembers,
                     double muh[ncomp][nconf], double SigmaSh[ncomp][totran*totran], double Dh[ncomp][totran*totran],
                     double Eh[ncomp][totran*totran], double beta[ncomp][totNreg], int *compAlloc, double pi[ncomp][n]){
    int j, h, komp, labelA, labelB;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double nonempty[ncomp];
    unsigned int vecAlloc[ncomp];
    double switchA, decisionA;

    double muhtemp[nconf];
    double vartemp[totran*totran];
    double betatemp[totNreg];
    int compAlloctemp[n];
    int nmemberstemp;

    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0)
            nonempty[h] = 1.0;
        else
            nonempty[h] = 0.0;
    }

    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);

    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelA = komp;

    nonempty[labelA] = 0.0;
    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);
    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelB = komp;

    if ((labelA >= ncomp) ||  (labelB >= ncomp))
        switchA = 0.0;
    else{
        switchA = 1.0;
        for (j = 0; j < n; j++){
            if (compAlloc[j]==labelA) switchA *= pi[labelB][j]/pi[labelA][j];
            if (compAlloc[j]==labelB) switchA *= pi[labelA][j]/pi[labelB][j];
        }
    }

    if (switchA > 1.0) switchA = 1.0;

    decisionA = gsl_ran_flat(r,0.0,1.0);

    if (switchA > decisionA){
        for (j = 0; j < nconf; j++) muhtemp[j] = muh[labelA][j];
        for (j = 0; j < nconf; j++) muh[labelA][j] = muh[labelB][j];
        for (j = 0; j < nconf; j++) muh[labelB][j] = muhtemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = SigmaSh[labelA][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelA][j] = SigmaSh[labelB][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelB][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Dh[labelA][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelA][j] = Dh[labelB][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelB][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Eh[labelA][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelA][j] = Eh[labelB][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelB][j] = vartemp[j];

        for (j = 0; j < totNreg; j++) betatemp[j] = beta[labelA][j];
        for (j = 0; j < totNreg; j++) beta[labelA][j] = beta[labelB][j];
        for (j = 0; j < totNreg; j++) beta[labelB][j] = betatemp[j];

        nmemberstemp = nmembers[labelA];
        nmembers[labelA] = nmembers[labelB];
        nmembers[labelB] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelA) compAlloc[j]=labelB;
            if (compAlloctemp[j] == labelB) compAlloc[j]=labelA;
        }
    }
    gsl_rng_free(r);
}

void SpatialLabelSwitchingAFG(unsigned long int s, int n, int totNreg, int nres, int totran, int ncomp, int *nmembers,
                     double SigmaSh[ncomp][totran*totran], double Dh[ncomp][totran*totran], double Eh[ncomp][totran*totran],
                     double beta[ncomp][totNreg], int *compAlloc, double pi[ncomp][n]){
    int j, h, komp, labelA, labelB;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double nonempty[ncomp];
    unsigned int vecAlloc[ncomp];
    double switchA, decisionA;

    double vartemp[totran*totran];
    double betatemp[totNreg];
    int compAlloctemp[n];
    int nmemberstemp;

    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0)
            nonempty[h] = 1.0;
        else
            nonempty[h] = 0.0;
    }

    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);

    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelA = komp;

    nonempty[labelA] = 0.0;
    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);
    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelB = komp;

    if ((labelA >= ncomp) ||  (labelB >= ncomp))
        switchA = 0.0;
    else{
        switchA = 1.0;
        for (j = 0; j < n; j++){
            if (compAlloc[j]==labelA) switchA *= pi[labelB][j]/pi[labelA][j];
            if (compAlloc[j]==labelB) switchA *= pi[labelA][j]/pi[labelB][j];
        }
    }

    if (switchA > 1.0) switchA = 1.0;

    decisionA = gsl_ran_flat(r,0.0,1.0);

    if (switchA > decisionA){
        for (j = 0; j < (totran*totran); j++) vartemp[j] = SigmaSh[labelA][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelA][j] = SigmaSh[labelB][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelB][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Dh[labelA][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelA][j] = Dh[labelB][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelB][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Eh[labelA][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelA][j] = Eh[labelB][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelB][j] = vartemp[j];

        for (j = 0; j < totNreg; j++) betatemp[j] = beta[labelA][j];
        for (j = 0; j < totNreg; j++) beta[labelA][j] = beta[labelB][j];
        for (j = 0; j < totNreg; j++) beta[labelB][j] = betatemp[j];

        nmemberstemp = nmembers[labelA];
        nmembers[labelA] = nmembers[labelB];
        nmembers[labelB] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelA) compAlloc[j]=labelB;
            if (compAlloctemp[j] == labelB) compAlloc[j]=labelA;
        }
    }
    gsl_rng_free(r);
}

void labelSwitchingB(unsigned long int s, int n, int nconf, int totNreg, int nres, int totran, int ncomp, int *nmembers,
                     double muh[ncomp][nconf], double SigmaSh[ncomp][totran*totran], double Dh[ncomp][totran*totran],
                     double Eh[ncomp][totran*totran], double beta[ncomp][totNreg], int *compAlloc, double *V){
    int j, h, komp, maxZ, labelC;
    maxZ = 0;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double equalProb, temp, switchB, decisionB;

    double muhtemp[nconf];
    double vartemp[totran*totran];
    double betatemp[totNreg];
    int compAlloctemp[n];
    int nmemberstemp;
    double Vtemp;

    for (h = 0; h < ncomp; h++)
        if (nmembers[h] > 0) maxZ = h;

    equalProb = 1.0/((double) maxZ);
    temp = gsl_ran_flat(r,0.0,1.0);

    komp=0;
    while(equalProb < temp){
        komp++;
        equalProb += 1.0/((double) maxZ);
    }
    labelC = komp;

    if (labelC >= (ncomp-1))
        switchB = 0.0;
    else
        switchB =  pow((1-V[labelC+1]),nmembers[labelC])/pow((1-V[labelC]),nmembers[labelC+1]);

    if (switchB > 1.0) switchB = 1.0;
    decisionB = gsl_ran_flat(r,0.0,1.0);

    if (switchB > decisionB){
        for (j = 0; j < nconf; j++) muhtemp[j] = muh[labelC][j];
        for (j = 0; j < nconf; j++) muh[labelC][j] = muh[labelC+1][j];
        for (j = 0; j < nconf; j++) muh[labelC+1][j] = muhtemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = SigmaSh[labelC][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelC][j] = SigmaSh[labelC+1][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Dh[labelC][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelC][j] = Dh[labelC+1][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Eh[labelC][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelC][j] = Eh[labelC+1][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelC+1][j] = vartemp[j];

        for (j = 0; j < totNreg; j++) betatemp[j] = beta[labelC][j];
        for (j = 0; j < totNreg; j++) beta[labelC][j] = beta[labelC+1][j];
        for (j = 0; j < totNreg; j++) beta[labelC+1][j] = betatemp[j];

        nmemberstemp = nmembers[labelC];
        nmembers[labelC] = nmembers[(labelC+1)];
        nmembers[(labelC+1)] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelC) compAlloc[j]=labelC+1;
            if (compAlloctemp[j] == (labelC+1)) compAlloc[j]=labelC;
        }

        Vtemp = V[labelC];
        V[labelC] = V[labelC+1];
        V[labelC+1] = Vtemp;
    }
    gsl_rng_free(r);
}

void SpatialLabelSwitchingB(unsigned long int s, int n, int nconf, int totNreg, int nres, int totran, int ncomp, int *nmembers,
                     double muh[ncomp][nconf], double SigmaSh[ncomp][totran*totran], double Dh[ncomp][totran*totran],
                     double Eh[ncomp][totran*totran], double beta[ncomp][totNreg], int *compAlloc, double P[ncomp][n], double eta[ncomp][n]){
    int j, h, komp, maxZ, labelC;
    maxZ = 0;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double equalProb, temp, switchB, decisionB;

    double muhtemp[nconf];
    double vartemp[totran*totran];
    double betatemp[totNreg];
    int compAlloctemp[n];
    int nmemberstemp;
    double etatemp[n];

    for (h = 0; h < ncomp; h++)
        if (nmembers[h] > 0) maxZ = h;

    equalProb = 1/((double) maxZ);
    temp = gsl_ran_flat(r,0.0,1.0);

    komp=0;
    while(equalProb < temp){
        komp++;
        equalProb += 1/((double) maxZ);
    }
    labelC = komp;

    if (labelC >= (ncomp-1))
        switchB = 0.0;
    else{
        switchB = 1.0;
        for (j = 0; j < n; j++){
            if (compAlloc[j]==labelC) switchB *= (1-P[labelC+1][j]);
            if (compAlloc[j]==(labelC+1)) switchB *= 1/(1-P[labelC][j]);
        }
    }

    if (switchB > 1.0) switchB = 1.0;
    decisionB = gsl_ran_flat(r,0.0,1.0);

    if (switchB > decisionB){
        for (j = 0; j < nconf; j++) muhtemp[j] = muh[labelC][j];
        for (j = 0; j < nconf; j++) muh[labelC][j] = muh[labelC+1][j];
        for (j = 0; j < nconf; j++) muh[labelC+1][j] = muhtemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = SigmaSh[labelC][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelC][j] = SigmaSh[labelC+1][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Dh[labelC][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelC][j] = Dh[labelC+1][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Eh[labelC][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelC][j] = Eh[labelC+1][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelC+1][j] = vartemp[j];

        for (j = 0; j < totNreg; j++) betatemp[j] = beta[labelC][j];
        for (j = 0; j < totNreg; j++) beta[labelC][j] = beta[labelC+1][j];
        for (j = 0; j < totNreg; j++) beta[labelC+1][j] = betatemp[j];

        nmemberstemp = nmembers[labelC];
        nmembers[labelC] = nmembers[(labelC+1)];
        nmembers[(labelC+1)] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelC) compAlloc[j]=labelC+1;
            if (compAlloctemp[j] == (labelC+1)) compAlloc[j]=labelC;
        }

        for (j = 0; j < n; j++) etatemp[j] = eta[labelC][j];
        for (j = 0; j < n; j++) eta[labelC][j] = eta[labelC+1][j];
        for (j = 0; j < n; j++) eta[labelC+1][j] = etatemp[j];
    }
    gsl_rng_free(r);
}

void SpatialLabelSwitchingBFG(unsigned long int s, int n, int totNreg, int nres, int totran, int ncomp, int *nmembers,
                     double SigmaSh[ncomp][totran*totran], double Dh[ncomp][totran*totran], double Eh[ncomp][totran*totran],
                     double beta[ncomp][totNreg], int *compAlloc, double P[ncomp][n], double eta[ncomp][n]){
    int j, h, komp, maxZ, labelC;
    maxZ = 0;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double equalProb, temp, switchB, decisionB;

    double vartemp[totran*totran];
    double betatemp[totNreg];
    int compAlloctemp[n];
    int nmemberstemp;
    double etatemp[n];

    for (h = 0; h < ncomp; h++)
        if (nmembers[h] > 0) maxZ = h;

    equalProb = 1/((double) maxZ);
    temp = gsl_ran_flat(r,0.0,1.0);

    komp=0;
    while(equalProb < temp){
        komp++;
        equalProb += 1/((double) maxZ);
    }
    labelC = komp;

    if (labelC >= (ncomp-1))
        switchB = 0.0;
    else{
        switchB = 1.0;
        for (j = 0; j < n; j++){
            if (compAlloc[j]==labelC) switchB *= (1-P[labelC+1][j]);
            if (compAlloc[j]==(labelC+1)) switchB *= 1/(1-P[labelC][j]);
        }
    }

    if (switchB > 1.0) switchB = 1.0;
    decisionB = gsl_ran_flat(r,0.0,1.0);

    if (switchB > decisionB){

        for (j = 0; j < (totran*totran); j++) vartemp[j] = SigmaSh[labelC][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelC][j] = SigmaSh[labelC+1][j];
        for (j = 0; j < (totran*totran); j++) SigmaSh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Dh[labelC][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelC][j] = Dh[labelC+1][j];
        for (j = 0; j < (totran*totran); j++) Dh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (totran*totran); j++) vartemp[j] = Eh[labelC][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelC][j] = Eh[labelC+1][j];
        for (j = 0; j < (totran*totran); j++) Eh[labelC+1][j] = vartemp[j];

        for (j = 0; j < totNreg; j++) betatemp[j] = beta[labelC][j];
        for (j = 0; j < totNreg; j++) beta[labelC][j] = beta[labelC+1][j];
        for (j = 0; j < totNreg; j++) beta[labelC+1][j] = betatemp[j];

        nmemberstemp = nmembers[labelC];
        nmembers[labelC] = nmembers[(labelC+1)];
        nmembers[(labelC+1)] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelC) compAlloc[j]=labelC+1;
            if (compAlloctemp[j] == (labelC+1)) compAlloc[j]=labelC;
        }

        for (j = 0; j < n; j++) etatemp[j] = eta[labelC][j];
        for (j = 0; j < n; j++) eta[labelC][j] = eta[labelC+1][j];
        for (j = 0; j < n; j++) eta[labelC+1][j] = etatemp[j];
    }
    gsl_rng_free(r);
}

void labelSwitchingAHannahNoSpace(unsigned long int s, int n, int nreg, int nres, int ncomp, int *nmembers, double *pi,
                     double muh[ncomp][nreg], double Sigmah[ncomp][nreg*nreg], double Dh[ncomp][nres*nres],
                     double Eh[ncomp][nres*nres], double Rh[ncomp][nres*nres], double beta[ncomp][(nreg+1)*nres], int *compAlloc){
    int j, h, komp, labelA, labelB;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double nonempty[ncomp];
    unsigned int vecAlloc[ncomp];
    double switchA, decisionA;

    double muhtemp[nreg];
    double sigmatemp[nreg*nreg];
    double betatemp[(nreg+1)*nres];
    double vartemp[nres*nres];
    int compAlloctemp[n];
    int nmemberstemp;

    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0)
            nonempty[h] = 1.0;
        else
            nonempty[h] = 0.0;
    }

    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);

    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelA = komp;

    nonempty[labelA] = 0.0;
    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);
    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelB = komp;

    if ((labelA >= ncomp) ||  (labelB >= ncomp))
        switchA = 0.0;
    else
        switchA = pow((pi[labelB]/pi[labelA]), (double)  (nmembers[labelA]-nmembers[labelB]));

    if (switchA > 1.0) switchA = 1.0;

    decisionA = gsl_ran_flat(r,0.0,1.0);

    if (switchA > decisionA){
        for (j = 0; j < nreg; j++) muhtemp[j] = muh[labelA][j];
        for (j = 0; j < nreg; j++) muh[labelA][j] = muh[labelB][j];
        for (j = 0; j < nreg; j++) muh[labelB][j] = muhtemp[j];

        for (j = 0; j < (nreg*nreg); j++) sigmatemp[j] = Sigmah[labelA][j];
        for (j = 0; j < (nreg*nreg); j++) Sigmah[labelA][j] = Sigmah[labelB][j];
        for (j = 0; j < (nreg*nreg); j++) Sigmah[labelB][j] = sigmatemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Rh[labelA][j];
        for (j = 0; j < (nres*nres); j++) Rh[labelA][j] = Rh[labelB][j];
        for (j = 0; j < (nres*nres); j++) Rh[labelB][j] = vartemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Dh[labelA][j];
        for (j = 0; j < (nres*nres); j++) Dh[labelA][j] = Dh[labelB][j];
        for (j = 0; j < (nres*nres); j++) Dh[labelB][j] = vartemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Eh[labelA][j];
        for (j = 0; j < (nres*nres); j++) Eh[labelA][j] = Eh[labelB][j];
        for (j = 0; j < (nres*nres); j++) Eh[labelB][j] = vartemp[j];

        for (j = 0; j < ((nreg+1)*nres); j++) betatemp[j] = beta[labelA][j];
        for (j = 0; j < ((nreg+1)*nres); j++) beta[labelA][j] = beta[labelB][j];
        for (j = 0; j < ((nreg+1)*nres); j++) beta[labelB][j] = betatemp[j];

        nmemberstemp = nmembers[labelA];
        nmembers[labelA] = nmembers[labelB];
        nmembers[labelB] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelA) compAlloc[j]=labelB;
            if (compAlloctemp[j] == labelB) compAlloc[j]=labelA;
        }
    }
    gsl_rng_free(r);
}

void labelSwitchingBHannahNoSpace(unsigned long int s, int n, int nreg, int nres, int ncomp, int *nmembers,
                     double muh[ncomp][nreg], double Sigmah[ncomp][nreg*nreg], double Dh[ncomp][nres*nres],
                     double Eh[ncomp][nres*nres], double Rh[ncomp][nres*nres], double beta[ncomp][(nreg+1)*nres], int *compAlloc, double *V){
    int j, h, komp, maxZ, labelC;
    maxZ = 0;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double equalProb, temp, switchB, decisionB;

    double muhtemp[nreg];
    double sigmatemp[nreg*nreg];
    double betatemp[(nreg+1)*nres];
    double vartemp[nres*nres];
    int compAlloctemp[n];
    int nmemberstemp;
    double Vtemp;

    for (h = 0; h < ncomp; h++)
        if (nmembers[h] > 0) maxZ = h;

    equalProb = 1/((double) maxZ);
    temp = gsl_ran_flat(r,0.0,1.0);

    komp=0;
    while(equalProb < temp){
        komp++;
        equalProb += 1/((double) maxZ);
    }
    labelC = komp;

    if (labelC >= (ncomp-1))
        switchB = 0.0;
    else
        switchB =  pow((1-V[labelC+1]),nmembers[labelC])/pow((1-V[labelC]),nmembers[labelC+1]);

    if (switchB > 1.0) switchB = 1.0;
    decisionB = gsl_ran_flat(r,0.0,1.0);

    if (switchB > decisionB){
        for (j = 0; j < nreg; j++) muhtemp[j] = muh[labelC][j];
        for (j = 0; j < nreg; j++) muh[labelC][j] = muh[labelC+1][j];
        for (j = 0; j < nreg; j++) muh[labelC+1][j] = muhtemp[j];

        for (j = 0; j < (nreg*nreg); j++) sigmatemp[j] = Sigmah[labelC][j];
        for (j = 0; j < (nreg*nreg); j++) Sigmah[labelC][j] = Sigmah[labelC+1][j];
        for (j = 0; j < (nreg*nreg); j++) Sigmah[labelC+1][j] = sigmatemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Rh[labelC][j];
        for (j = 0; j < (nres*nres); j++) Rh[labelC][j] = Rh[labelC+1][j];
        for (j = 0; j < (nres*nres); j++) Rh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Dh[labelC][j];
        for (j = 0; j < (nres*nres); j++) Dh[labelC][j] = Dh[labelC+1][j];
        for (j = 0; j < (nres*nres); j++) Dh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Eh[labelC][j];
        for (j = 0; j < (nres*nres); j++) Eh[labelC][j] = Eh[labelC+1][j];
        for (j = 0; j < (nres*nres); j++) Eh[labelC+1][j] = vartemp[j];

        for (j = 0; j < ((nreg+1)*nres); j++) betatemp[j] = beta[labelC][j];
        for (j = 0; j < ((nreg+1)*nres); j++) beta[labelC][j] = beta[labelC+1][j];
        for (j = 0; j < ((nreg+1)*nres); j++) beta[labelC+1][j] = betatemp[j];

        nmemberstemp = nmembers[labelC];
        nmembers[labelC] = nmembers[labelC+1];
        nmembers[labelC+1] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelC) compAlloc[j]=labelC+1;
            if (compAlloctemp[j] == (labelC+1)) compAlloc[j]=labelC;
        }

        Vtemp = V[labelC];
        V[labelC] = V[labelC+1];
        V[labelC+1] = Vtemp;
    }
    gsl_rng_free(r);
}

void labelSwitchingAHannahSpace(unsigned long int s, int n, int nreg, int nres, int ncomp, int *nmembers, double pi[ncomp][n],
                     double muh[ncomp][nreg], double Sigmah[ncomp][nreg*nreg], double Dh[ncomp][nres*nres],
                     double Eh[ncomp][nres*nres], double Rh[ncomp][nres*nres], double beta[ncomp][(nreg+1)*nres], int *compAlloc){
    int j, h, komp, labelA, labelB;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double nonempty[ncomp];
    unsigned int vecAlloc[ncomp];
    double switchA, decisionA;

    double muhtemp[nreg];
    double sigmatemp[nreg*nreg];
    double betatemp[(nreg+1)*nres];
    double vartemp[nres*nres];
    int compAlloctemp[n];
    int nmemberstemp;

    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0)
            nonempty[h] = 1.0;
        else
            nonempty[h] = 0.0;
    }

    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);

    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelA = komp;

    nonempty[labelA] = 0.0;
    gsl_ran_multinomial(r,ncomp,1,nonempty,vecAlloc);
    komp=0;
    while(vecAlloc[komp]==0)
        komp++;
    labelB = komp;

    if ((labelA >= ncomp) ||  (labelB >= ncomp))
        switchA = 0.0;
    else{
        switchA = 1.0;
        for (j = 0; j < n; j++){
            if (compAlloc[j]==labelA) switchA *= pi[labelB][j]/pi[labelA][j];
            if (compAlloc[j]==labelB) switchA *= pi[labelA][j]/pi[labelB][j];
        }
    }

    if (switchA > 1.0) switchA = 1.0;

    decisionA = gsl_ran_flat(r,0.0,1.0);

    if (switchA > decisionA){
        for (j = 0; j < nreg; j++) muhtemp[j] = muh[labelA][j];
        for (j = 0; j < nreg; j++) muh[labelA][j] = muh[labelB][j];
        for (j = 0; j < nreg; j++) muh[labelB][j] = muhtemp[j];

        for (j = 0; j < (nreg*nreg); j++) sigmatemp[j] = Sigmah[labelA][j];
        for (j = 0; j < (nreg*nreg); j++) Sigmah[labelA][j] = Sigmah[labelB][j];
        for (j = 0; j < (nreg*nreg); j++) Sigmah[labelB][j] = sigmatemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Rh[labelA][j];
        for (j = 0; j < (nres*nres); j++) Rh[labelA][j] = Rh[labelB][j];
        for (j = 0; j < (nres*nres); j++) Rh[labelB][j] = vartemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Dh[labelA][j];
        for (j = 0; j < (nres*nres); j++) Dh[labelA][j] = Dh[labelB][j];
        for (j = 0; j < (nres*nres); j++) Dh[labelB][j] = vartemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Eh[labelA][j];
        for (j = 0; j < (nres*nres); j++) Eh[labelA][j] = Eh[labelB][j];
        for (j = 0; j < (nres*nres); j++) Eh[labelB][j] = vartemp[j];

        for (j = 0; j < ((nreg+1)*nres); j++) betatemp[j] = beta[labelA][j];
        for (j = 0; j < ((nreg+1)*nres); j++) beta[labelA][j] = beta[labelB][j];
        for (j = 0; j < ((nreg+1)*nres); j++) beta[labelB][j] = betatemp[j];

        nmemberstemp = nmembers[labelA];
        nmembers[labelA] = nmembers[labelB];
        nmembers[labelB] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelA) compAlloc[j]=labelB;
            if (compAlloctemp[j] == labelB) compAlloc[j]=labelA;
        }
    }
    gsl_rng_free(r);
}

void labelSwitchingBHannahSpace(unsigned long int s, int n, int nreg, int nres, int ncomp, int *nmembers,
                     double muh[ncomp][nreg], double Sigmah[ncomp][nreg*nreg], double Dh[ncomp][nres*nres],
                     double Eh[ncomp][nres*nres], double Rh[ncomp][nres*nres], double beta[ncomp][(nreg+1)*nres],
                     int *compAlloc, double P[ncomp][n], double eta[ncomp][n]){
    int j, h, komp, maxZ, labelC;
    maxZ = 0;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double equalProb, temp, switchB, decisionB;

    double muhtemp[nreg];
    double sigmatemp[nreg*nreg];
    double betatemp[(nreg+1)*nres];
    double vartemp[nres*nres];
    int compAlloctemp[n];
    int nmemberstemp;
    double etatemp[n];

    for (h = 0; h < ncomp; h++)
        if (nmembers[h] > 0) maxZ = h;

    equalProb = 1/((double) maxZ);
    temp = gsl_ran_flat(r,0.0,1.0);

    komp=0;
    while(equalProb < temp){
        komp++;
        equalProb += 1/((double) maxZ);
    }
    labelC = komp;

    if (labelC >= (ncomp-1))
        switchB = 0.0;
    else{
        switchB = 1.0;
        for (j = 0; j < n; j++){
            if (compAlloc[j]==labelC) switchB *= (1-P[labelC+1][j]);
            if (compAlloc[j]==(labelC+1)) switchB *= 1/(1-P[labelC][j]);
        }
    }

    if (switchB > 1.0) switchB = 1.0;
    decisionB = gsl_ran_flat(r,0.0,1.0);

    if (switchB > decisionB){
        for (j = 0; j < nreg; j++) muhtemp[j] = muh[labelC][j];
        for (j = 0; j < nreg; j++) muh[labelC][j] = muh[labelC+1][j];
        for (j = 0; j < nreg; j++) muh[labelC+1][j] = muhtemp[j];

        for (j = 0; j < (nreg*nreg); j++) sigmatemp[j] = Sigmah[labelC][j];
        for (j = 0; j < (nreg*nreg); j++) Sigmah[labelC][j] = Sigmah[labelC+1][j];
        for (j = 0; j < (nreg*nreg); j++) Sigmah[labelC+1][j] = sigmatemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Rh[labelC][j];
        for (j = 0; j < (nres*nres); j++) Rh[labelC][j] = Rh[labelC+1][j];
        for (j = 0; j < (nres*nres); j++) Rh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Dh[labelC][j];
        for (j = 0; j < (nres*nres); j++) Dh[labelC][j] = Dh[labelC+1][j];
        for (j = 0; j < (nres*nres); j++) Dh[labelC+1][j] = vartemp[j];

        for (j = 0; j < (nres*nres); j++) vartemp[j] = Eh[labelC][j];
        for (j = 0; j < (nres*nres); j++) Eh[labelC][j] = Eh[labelC+1][j];
        for (j = 0; j < (nres*nres); j++) Eh[labelC+1][j] = vartemp[j];

        for (j = 0; j < ((nreg+1)*nres); j++) betatemp[j] = beta[labelC][j];
        for (j = 0; j < ((nreg+1)*nres); j++) beta[labelC][j] = beta[labelC+1][j];
        for (j = 0; j < ((nreg+1)*nres); j++) beta[labelC+1][j] = betatemp[j];

        nmemberstemp = nmembers[labelC];
        nmembers[labelC] = nmembers[labelC+1];
        nmembers[labelC+1] = nmemberstemp;

        for (j = 0; j < n; j++)
            compAlloctemp[j] = compAlloc[j];

        for (j = 0; j < n; j++){
            if (compAlloctemp[j] == labelC) compAlloc[j]=labelC+1;
            if (compAlloctemp[j] == (labelC+1)) compAlloc[j]=labelC;
        }

        for (j = 0; j < n; j++) etatemp[j] = eta[labelC][j];
        for (j = 0; j < n; j++) eta[labelC][j] = eta[labelC+1][j];
        for (j = 0; j < n; j++) eta[labelC+1][j] = etatemp[j];
    }
    gsl_rng_free(r);
}

double numberOfClusters(int ncomp, int *nmembers){
    int h;
    double ns; ns = 0.0;

    for (h = 0; h < ncomp; h++)
        if (nmembers[h] > 0)
            ns += 1.0;

    return(ns);
}

void calcSpatialSums(int n, int ncomp, int *nmembers, double alphau, double eta[ncomp][n], double *BS, gsl_matrix *copyqij){
    int h, j, r, c;
    double BS1, BS2;
    BS1 = 0.0; BS2 = 0.0;

    for (h = 0; h < ncomp; h++){
        if (nmembers[h] > 0){
            for (j = 0; j < n; j++)
                BS1 += (eta[h][j]-alphau)*(eta[h][j]-alphau);
            for (r = 0; r < (n-1); r++)
                for (c = (r+1); c < n; c++)
                     if (gsl_matrix_get(copyqij,r,c)==-1) BS2 += (eta[h][r]-eta[h][c])*(eta[h][r]-eta[h][c]);
        }
    }
    BS[0] = BS1; BS[1] = BS2;
}


//Sum of squares of elements of vector x that belong to cluster h
double SSQh(int n, int h, int *compAlloc, double *x){
    int i;
    double sumssq = 0.0;
    for (i = 0; i < n; i++){
        if (compAlloc[i]==h){
            sumssq += x[i]*x[i];
        }
    }
return(sumssq);
}

//Print gsl matrix
void print_matrix(gsl_matrix *A)
{
    int i, j;
    for (i = 0; i < A->size1; i++) {
        for (j = 0; j < A->size2; j++)
            Rprintf("%g\t", gsl_matrix_get(A, i, j));
        Rprintf("\n");
    }
    Rprintf("\n");
}

//Print gsl matrix
void print_vector(gsl_vector *V)
{
    int i;
    for (i = 0; i < V->size; i++) 
        Rprintf("%g\t", V->data[i * V->stride]);
    Rprintf("\n");
}
