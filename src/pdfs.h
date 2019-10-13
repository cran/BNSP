/* pdfs.h Implements probability density functions
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

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

//Bivariate Normal pdf with arguments (dim=2,x,parameters,1,return.value)
int bivNormalpdf(unsigned dim, const double *x, void *parameters, unsigned fdim, double *fval){
    double *params = parameters;
    double mu0 = params[0];
    double mu1 = params[1];
    double sigma1 = params[2];
    double sigma2 = params[3];
    double rho = params[4];
    double exponent;

    double sigma1s, sigma2s;           //sigma1^2 and sigma2^2
    sigma1s = sigma1*sigma1;
    sigma2s = sigma2*sigma2;

    //double eigenval1;
    double eigenval2;
    //eigenval1 = (sigma1s + sigma2s)/2 + sqrt(4*sigma1s*sigma2s*rho*rho + (sigma1s-sigma2s)*(sigma1s-sigma2s))/2;
    eigenval2 = (sigma1s + sigma2s)/2 - sqrt(4*sigma1s*sigma2s*rho*rho + (sigma1s-sigma2s)*(sigma1s-sigma2s))/2;

    exponent = (1/(1-rho*rho))*(
               ((x[0]-mu0)*(x[0]-mu0)/(sigma1*sigma1)) +
               ((x[1]-mu1)*(x[1]-mu1)/(sigma2*sigma2)) -
               (2*rho*(x[0]-mu0)*(x[1]-mu1)/(sigma1*sigma2)));

    fval[0] = exp(-exponent/2)/(2*M_PI*sigma1*sigma2*sqrt(1-rho*rho));
    if (eigenval2 < 0.001) fval[0] = 0;

    return 0;
}

//params are means, diag elems of Sigma^{-1}, off diag elems of Sigma^{-1}, log|Sigma^-1| 
int MultiNormalPDF(unsigned dim, const double *x, void *parameters, unsigned fdim, double *fval){
    double *params = parameters;
    int ncovs = dim*(dim-1)/2; 
    int i, j, k;
    double mu[dim];
    double vars[dim];
    double covs[ncovs+1];//add 1 to avoid problems with dim 0
    double logdet, Q, logval;
    for (i = 0; i < dim; i++) 
        mu[i] = params[i];    
    for (i = 0; i < dim; i++) 
        vars[i] = params[dim+i];
    for (i = 0; i < ncovs; i++) 
        covs[i] = params[2*dim+i];
    logdet = params[2*dim+ncovs]; 
    Q = 0.0;
    for (i = 0; i < dim; i++)
        Q += vars[i]*(x[i]-mu[i])*(x[i]-mu[i]);
    k = 0;
    for (i = 0; i < (dim-1); i++)
        for (j = i+1; j < dim; j++)
            Q += 2*covs[k++]*(x[i]-mu[i])*(x[j]-mu[j]);
    logval = -(dim*log(2*M_PI)+logdet+Q)/2; 
    fval[0] = exp(logval);
    return 0;
}

//Generalized bivariate Normal pdf with arguments (dim=2,x,parameters,1,return.value)
int gBivNormalpdf(unsigned dim, const double *x, void *parameters, unsigned fdim, double *fval){
    double tol = 0.001; //tolerance level for the smallest eigenvalue.
    double *params = parameters;
    double mu0 = params[0];
    double mu1 = params[1];
    double sigma1 = params[2];
    double sigma2 = params[3];
    double rho = params[4];

    double sigma1s, sigma2s;           //sigma1^2 and sigma2^2
    sigma1s = sigma1*sigma1;
    sigma2s = sigma2*sigma2;

    double eigenval1, eigenval2;       //store the 2 eigen values
    eigenval1 = (sigma1s + sigma2s)/2 + sqrt(4*sigma1s*sigma2s*rho*rho + (sigma1s-sigma2s)*(sigma1s-sigma2s))/2;
    eigenval2 = (sigma1s + sigma2s)/2 - sqrt(4*sigma1s*sigma2s*rho*rho + (sigma1s-sigma2s)*(sigma1s-sigma2s))/2;

    double eigenvec11, eigenvec12, l1; // store 2 elements of eigenvector 1

    if (rho==0){
        if (sigma1 > sigma2) {
            eigenvec11 = 1;
            eigenvec12 = 0;
        }
        else{
            eigenvec11 = 0;
            eigenvec12 = 1;
        }
    }
    else{
        eigenvec11 = 1;
        eigenvec12 = (eigenval1-sigma1s)/(sigma1*sigma2*rho);
        l1 = sqrt(eigenvec11 + eigenvec12*eigenvec12);
        eigenvec11 *= 1/l1;
        eigenvec12 *= 1/l1;
    }

    double eigenvec21, eigenvec22, l2; // store 2 elements of eigenvector 2

    if (rho==0){
        if (sigma1 > sigma2) {
            eigenvec21 = 0;
            eigenvec22 = 1;
        }
        else{
            eigenvec21 = 1;
            eigenvec22 = 0;
        }
    }
    else{
        eigenvec21 = 1;
        eigenvec22 = (eigenval2-sigma1s)/(sigma1*sigma2*rho);
        l2 = sqrt(eigenvec21 + eigenvec22*eigenvec22);
        eigenvec21 *= 1/l2;
        eigenvec22 *= 1/l2;
    }

    double inv11, inv12, inv22;
    inv11 = eigenvec11*eigenvec11/eigenval1;
    inv12 = eigenvec11*eigenvec12/eigenval1;
    inv22 = eigenvec12*eigenvec12/eigenval1;
    if (eigenval2 > tol){
        inv11 += eigenvec21*eigenvec21/eigenval2;
        inv12 += eigenvec21*eigenvec22/eigenval2;
        inv22 += eigenvec22*eigenvec22/eigenval2;
    }

    double det;
    det = eigenval1;
    if (eigenval2 > tol) det *= eigenval2;

    double exponent;
    exponent = (x[0]-mu0)*(x[0]-mu0)*inv11 +
               (x[1]-mu1)*(x[1]-mu1)*inv22 +
               2*(x[0]-mu0)*(x[1]-mu1)*inv12;

    fval[0] = exp(-exponent/2)/(2*M_PI*sqrt(det));

    return 0;
}

//Multivariate Normal pdf with arguments (dim,x,mean,Sigma,tol) f_p(x|mu,Sigma)
double logMVNormalpdf(int dim, gsl_vector *x, gsl_vector *mu, gsl_matrix *S, double tol){
    int i;
    double det, temp, QF;

    gsl_vector *copyX = gsl_vector_alloc(dim);
    gsl_matrix *copyS = gsl_matrix_alloc(dim,dim);
    gsl_vector *eval = gsl_vector_alloc(dim);
    gsl_matrix *evec = gsl_matrix_alloc(dim,dim);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(dim);
    gsl_matrix *D = gsl_matrix_calloc(dim,dim);
    gsl_matrix *VD = gsl_matrix_alloc(dim,dim);
    gsl_matrix *Sinv = gsl_matrix_alloc(dim,dim);
    gsl_vector *Prod1 = gsl_vector_alloc(dim);

    gsl_vector_memcpy(copyX,x);
    gsl_vector_sub(copyX,mu);
    gsl_matrix_memcpy(copyS,S);
    gsl_eigen_symmv(copyS,eval,evec,w);

    det = 1.0;
    for (i=0; i < dim; i++){
        temp = gsl_vector_get(eval,i);
	    if (temp > tol){
	        gsl_matrix_set(D,i,i,1/temp);
	        det *= temp;
	    }
	    else
	        gsl_matrix_set(D,i,i,0.0);
    }

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,D,0.0,VD);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,VD,evec,0.0,Sinv);
    gsl_blas_dgemv(CblasNoTrans,1.0,Sinv,copyX,0.0,Prod1);
    gsl_blas_ddot(copyX,Prod1,&QF);

    gsl_vector_free(copyX); gsl_matrix_free(copyS);
    gsl_vector_free(eval); gsl_matrix_free(evec); gsl_eigen_symmv_free(w);
    gsl_matrix_free(D); gsl_matrix_free(VD); gsl_matrix_free(Sinv);
    gsl_vector_free(Prod1);

    return(-0.5*QF-0.5*log(det)-((double)dim/2.0)*log(2*M_PI));
}

//Multivariate Normal pdf with arguments (dim,x,mean,Sigma) f_p(x|mu,Sigma)
double logMVNormalpdf2(int dim, gsl_vector *x, gsl_vector *mu, gsl_matrix *S){
    double det, QF;

    gsl_vector *copyX = gsl_vector_alloc(dim);
    gsl_matrix *copyS = gsl_matrix_alloc(dim,dim);
    gsl_matrix *Sinv = gsl_matrix_alloc(dim,dim);
    gsl_vector *Prod1 = gsl_vector_alloc(dim);
    gsl_permutation *p = gsl_permutation_alloc(dim); // Permutation
    int sign;                                        // Sign of the permutation

    gsl_vector_memcpy(copyX,x);
    gsl_vector_sub(copyX,mu);
    gsl_matrix_memcpy(copyS,S);
    gsl_linalg_LU_decomp(copyS,p,&sign);
    det = gsl_linalg_LU_det(copyS,sign);
    gsl_linalg_LU_invert(copyS,p,Sinv);
    gsl_blas_dgemv(CblasNoTrans,1.0,Sinv,copyX,0.0,Prod1);
    gsl_blas_ddot(copyX,Prod1,&QF);

    gsl_vector_free(copyX); gsl_matrix_free(copyS); gsl_matrix_free(Sinv);
    gsl_vector_free(Prod1); gsl_permutation_free(p);
    return(-0.5*QF-0.5*log(det)-((double)dim/2.0)*log(2*M_PI));
}

//Multivariate Normal pdf with arguments (dim,x,mean,Sigma^{-1}) f_p(x|mu,Sigma^{-1})
double logMVNormalpdf3(int dim, gsl_vector *x, gsl_vector *mu, gsl_matrix *S){
    double det, QF;
    gsl_vector *copyX = gsl_vector_alloc(dim);
    gsl_matrix *copyS = gsl_matrix_alloc(dim,dim);
    gsl_vector *Prod1 = gsl_vector_alloc(dim);
    gsl_permutation *p = gsl_permutation_alloc(dim); // Permutation
    int sign;                                        // Sign of the permutation

    gsl_vector_memcpy(copyX,x);
    gsl_vector_sub(copyX,mu);
    gsl_matrix_memcpy(copyS,S);
    gsl_linalg_LU_decomp(copyS,p,&sign);
    det = gsl_linalg_LU_det(copyS,sign);
    gsl_blas_dgemv(CblasNoTrans,1.0,S,copyX,0.0,Prod1);
    gsl_blas_ddot(copyX,Prod1,&QF);

    gsl_vector_free(copyX); gsl_matrix_free(copyS); //gsl_matrix_free(Sinv);
    gsl_vector_free(Prod1); gsl_permutation_free(p);
    return(-0.5*QF+0.5*log(det)-((double)dim/2.0)*log(2*M_PI));
}

//Joint posterior pdf of D_h and Sigma^*_h: log of
//|D|^(eta/2-1) |Sigma|^((eta-s-1-nmembers)/2) etr(-0.5*(eta Hinv Eh + SigmaInv Sh))
//nres denotes the number of non-identifiable varianes and nconf the number of identiable ones
double logPostPdfDSigma(gsl_matrix *D, gsl_matrix *Sigma, gsl_matrix *Eh, gsl_matrix *Hinv, gsl_matrix *Sh,
                    int nres, int nconf, int nmembers, double eta){
    int i, dim;
    double detD, detSigma, trace, result, temp;
    dim = nres+nconf;

    gsl_matrix *CopySigma = gsl_matrix_alloc(dim,dim);
    gsl_matrix *HinvEh = gsl_matrix_alloc(dim,dim);
    gsl_matrix *EigenD = gsl_matrix_calloc(dim,dim);
    gsl_vector *eval = gsl_vector_alloc(dim);
    gsl_matrix *evec = gsl_matrix_alloc(dim,dim);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(dim);
    gsl_matrix *VD = gsl_matrix_alloc(dim,dim);
    gsl_matrix *VDVT = gsl_matrix_alloc(dim,dim);

    gsl_matrix_memcpy(CopySigma,Sigma);

    detD = 1.0;
    for (i = 0; i < nres; i++) detD *= gsl_matrix_get(D,i,i);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Hinv,Eh,0.0,HinvEh);

    gsl_eigen_symmv(CopySigma,eval,evec,w);
    detSigma = 1.0;
    for (i=0; i < dim; i++){
        temp = gsl_vector_get(eval,i);
        detSigma *= temp;
	    gsl_matrix_set(EigenD,i,i,1/temp);
    }
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,EigenD,0.0,VD);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,VD,evec,0.0,VDVT);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,VDVT,Sh,0.0,VD);

    gsl_matrix_add(VD,HinvEh);

    trace = 0.0;
    for (i=0; i < dim; i++) trace += gsl_matrix_get(VD,i,i);

    result = ((eta/2.0)-1.0)*log(detD) + ((double)(eta-dim-1-nmembers)/2.0)*log(detSigma) - 0.5*trace;

    gsl_matrix_free(CopySigma); gsl_matrix_free(HinvEh);
    gsl_matrix_free(EigenD); gsl_vector_free(eval); gsl_matrix_free(evec); gsl_eigen_symmv_free(w);
    gsl_matrix_free(VD); gsl_matrix_free(VDVT);

    return(result);
}

//Computes log of |D|^(num1/2) |E|^(-num2/2) |K|^(num3/2) etr(-num4*E^{-1}*M/2) where D is diagonal 
double logPropPdfDR(gsl_matrix *D, gsl_matrix *E, gsl_matrix *M, gsl_matrix *K, int p, double num1, double num2, double num3, double num4){
    int i;
    double temp, detD, detE, detK, trace, result;
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(p);
    gsl_matrix *CopyE = gsl_matrix_alloc(p,p);
    gsl_vector *eval = gsl_vector_alloc(p);
    gsl_matrix *evec = gsl_matrix_alloc(p,p);
    gsl_matrix *EigenE = gsl_matrix_calloc(p,p);
    gsl_matrix *VD = gsl_matrix_alloc(p,p);
    gsl_matrix *VDVT = gsl_matrix_alloc(p,p);

    detD = 1.0;
    if (fabs(num1) > 0) 
        for (i = 0; i < p; i++) 
	        detD *= gsl_matrix_get(D,i,i);
    
    gsl_matrix_memcpy(CopyE,E);
    gsl_eigen_symmv(CopyE,eval,evec,w);
    detE = 1.0;
    for (i=0; i < p; i++){
        temp = gsl_vector_get(eval,i);
        detE *= temp;
        if (fabs(num4) > 0) gsl_matrix_set(EigenE,i,i,1/temp);
    }
    trace = 0.0;
    if (fabs(num4) > 0){
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,EigenE,0.0,VD);
        gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,VD,evec,0.0,VDVT);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,VDVT,M,0.0,VD);    
        for (i = 0; i < p; i++) 
            trace += gsl_matrix_get(VD,i,i);
    }
    detK = 1.0;
    if (fabs(num3) > 0){
        gsl_eigen_symmv_workspace * w2 = gsl_eigen_symmv_alloc(p);
        gsl_matrix *CopyK = gsl_matrix_alloc(p,p);
        gsl_vector *eval2 = gsl_vector_alloc(p);
        gsl_matrix *evec2 = gsl_matrix_alloc(p,p);
        gsl_matrix_memcpy(CopyK,K);
        gsl_eigen_symmv(CopyK,eval2,evec2,w2);    

        for (i=0; i < p; i++)          
            detK *= gsl_vector_get(eval2,i);
        
        gsl_eigen_symmv_free(w2); gsl_matrix_free(CopyK); gsl_vector_free(eval2); gsl_matrix_free(evec2);   
    }
    
    result = (num1/2.0)*log(detD) + (num3/2.0)*log(detK) - (num2/2.0)*log(detE) - num4*0.5*trace;

    gsl_eigen_symmv_free(w); gsl_matrix_free(CopyE); gsl_vector_free(eval); gsl_matrix_free(evec);
    gsl_matrix_free(EigenE); gsl_matrix_free(VD); gsl_matrix_free(VDVT);        
    
    return(result);
}

//Log Transition Ratio: current matrix, proposed matrix, number of non-identifiable variances, number of identifiable variances,
//degrees of freedom
double logtrnsR(gsl_matrix *Ehc, gsl_matrix *Ehp, int nres, int nconf, double nu){
    int i, dim;
    double detEhc, detEhp, tempc, tempp, trace, prod, result;
    dim = nres+nconf;

    gsl_matrix *EigenDc = gsl_matrix_calloc(dim,dim);
    gsl_matrix *EigenDp = gsl_matrix_calloc(dim,dim);

    gsl_vector *evalc = gsl_vector_alloc(dim);
    gsl_matrix *evecc = gsl_matrix_alloc(dim,dim);
    gsl_eigen_symmv_workspace * wc = gsl_eigen_symmv_alloc(dim);

    gsl_vector *evalp = gsl_vector_alloc(dim);
    gsl_matrix *evecp = gsl_matrix_alloc(dim,dim);
    gsl_eigen_symmv_workspace * wp = gsl_eigen_symmv_alloc(dim);

    gsl_matrix *CopyEhc = gsl_matrix_alloc(dim,dim);
    gsl_matrix *CopyEhp = gsl_matrix_alloc(dim,dim);

    gsl_matrix *VD = gsl_matrix_alloc(dim,dim);
    gsl_matrix *VD2 = gsl_matrix_alloc(dim,dim);
    gsl_matrix *VDVT = gsl_matrix_alloc(dim,dim);

    gsl_matrix_memcpy(CopyEhc,Ehc);
    gsl_matrix_memcpy(CopyEhp,Ehp);

    gsl_eigen_symmv(CopyEhc,evalc,evecc,wc);
    gsl_eigen_symmv(CopyEhp,evalp,evecp,wp);

    detEhc = 1.0;
    detEhp = 1.0;
    for (i=0; i < dim; i++){
        tempc = gsl_vector_get(evalc,i);
        tempp = gsl_vector_get(evalp,i);
        detEhc *= tempc;
        detEhp *= tempp;
	    gsl_matrix_set(EigenDc,i,i,1/tempc);
	    gsl_matrix_set(EigenDp,i,i,1/tempp);
    }

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evecc,EigenDc,0.0,VD);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,VD,evecc,0.0,VDVT);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,VDVT,Ehp,0.0,VD); //Ehc^{-1} Ehp

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evecp,EigenDp,0.0,VD2);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,VD2,evecp,0.0,VDVT);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,VDVT,Ehc,0.0,VD2); //Ehp^{-1} Ehc

    gsl_matrix_sub(VD,VD2);

    trace = 0.0;
    for (i=0; i < dim; i++) trace += gsl_matrix_get(VD,i,i);

    prod = 0.0;
    for (i=0; i < nres; i++) prod += log(gsl_matrix_get(Ehc,i,i)/gsl_matrix_get(Ehp,i,i));

    result = (-nu+((double)dim+1.0)/2.0)*(log(detEhp)-log(detEhc)) + (nu/2.0)*trace + ((dim-1.0)/2.0)*prod;

    gsl_matrix_free(EigenDc); gsl_matrix_free(EigenDp);
    gsl_vector_free(evalc); gsl_matrix_free(evecc); gsl_eigen_symmv_free(wc);
    gsl_vector_free(evalp); gsl_matrix_free(evecp); gsl_eigen_symmv_free(wp);
    gsl_matrix_free(CopyEhc); gsl_matrix_free(CopyEhp);
    gsl_matrix_free(VD); gsl_matrix_free(VD2); gsl_matrix_free(VDVT);

    return(result);
}

//Beta-Binomial cdf
double cdf_beta_binomial_P(int n, int q, double a, double b){
    double result = 0.0;
    if (q >= n) result = 1.0;
    else {
        int j;
        double Bab = gsl_sf_lnbeta(a,b);
        for (j = 0; j < q+1; j++) result += exp(gsl_sf_lnchoose(n,j) + gsl_sf_lnbeta(j+a,n-j+b) - Bab);
        if (result > 1.0) result = 1.0;
    }
    return(result);
}

double cdf_beta_binomial_P2(int n, int q, double a, double b){
    double result, prob;
    int j = 0;
    double Tot = 0.0;    
    prob = exp(gsl_sf_lnbeta(a,n+b) - gsl_sf_lnbeta(a,b));
    Tot = prob;
    for (j = 0; j < q; j++){
        prob *= ((j+a)*(n-j)/((j+1)*(n-j-1+b)));
        Tot += prob;            
	}
	result = Tot;        
    if (result > 1.0) result=1.0;
    return(result);
}

//A Generalized-Poisson cdf
double cdf_generalized_poisson_P1(int q, double lambda1, double lambda2){
    int j = 0;
    double B = -lambda1/lambda2; //for making sure that pmf is always >=0, (3.1) of Consul and Jain 1973
    double result = 0.0;
    if (lambda2 >= 0) for (j = 0; j < q+1; j++) result += exp((j-1)*log(lambda1+j*lambda2)-(lambda1+j*lambda2)-gsl_sf_lnfact(j));
    else if (lambda2 < 0 ){
        while (j < q+1 && j < B){
            result += exp((j-1)*log(lambda1+j*lambda2)-(lambda1+j*lambda2)-gsl_sf_lnfact(j));
            j++;
        }
    }
    return(lambda1*result);
}

//A Generalized-Poisson cdf that can include offset
double cdf_generalized_poisson_P2(int q, double mu, double f){
    double normConst = 0.0; //for when f < 1 normalization is needed
    double result;
    int j = 0;
    double B = -mu/(f-1); //for making sure that pmf is always >=0, (2.3) of Consul and Famoye 1992
    double Tot = 0.0;
    if (f == 1) Tot = gsl_cdf_poisson_P(q,mu);
    else if (f > 1) for (j = 0; j < q+1; j++) Tot += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
    else if (f < 1 ){
        while (j < q+1 && j < B){
            if ((mu+(f-1)*j)>0) Tot += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
            j++;
        }
        normConst = Tot;
        while (j < B){
            if ((mu+(f-1)*j)>0) normConst += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
            j++;
        }
    }
    if (f>1 && Tot>1) Tot=1.0;
    if (f>=1) result = Tot;
    else result = Tot/normConst;
    return(result);
}

//A Generalized-Poisson cdf that can include offset and doesn't always calculate norm const
double cdf_generalized_poisson_P3(int q, double mu, double f, double *normConst){
    *normConst = 1.0;
    double result, prob;
    int j = 0;
    double B = -mu/(f-1); //for making sure that pmf is always >=0, (2.3) of Consul and Famoye 1992
    double Tot = 0.0;    
    if (f == 1.0){
        Tot = gsl_cdf_poisson_P(q,mu);
        result = Tot;
    }
    else if (f > 1.0){        
        prob = exp(-mu/f);
        Tot = prob;
        for (j = 1; j < q+1; j++){
            prob *= exp((j-1)*log(mu+(f-1)*j)-(j-2)*log(mu+(f-1)*(j-1))-log(f)+(1-f)/f-log(j));
            Tot +=prob;            
		}
        result = Tot;        
    }
    else if ((mu < 1.0 && f < 0.97) || (mu < 2.0 && f < 0.80) || (mu < 3.0 && f < 0.65) || (mu < 4.0 && f < 0.60) ||
        (mu < 5.0 && f < 0.55))
    {
		prob = exp(-mu/f);
        Tot = prob;
        j = 1;
        while (j < q+1 && j < B){
            //if ((mu+(f-1)*j)>0){ 
            prob *= exp((j-1)*log(mu+(f-1)*j)-(j-2)*log(mu+(f-1)*(j-1))-log(f)+(1-f)/f-log(j));
            Tot +=prob;
            //Tot += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
			//}
            j++;
        }        
        *normConst = Tot;
        while (j < B){
            //if ((mu+(f-1)*j)>0){ 
            //normConst += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
            prob *= exp((j-1)*log(mu+(f-1)*j)-(j-2)*log(mu+(f-1)*(j-1))-log(f)+(1-f)/f-log(j));
            *normConst +=prob;
			//}
            j++;
        }
        result = Tot / *normConst;
        //Rprintf("%f %f %i %i %f \n",Tot/normConst2,Tot2/normConst,j,q+1,B);
    }
    else
    {
        prob = exp(-mu/f);
        Tot = prob;
        j = 1;
        while (j < q+1 && j < B){
            //if ((mu+(f-1)*j)>0) 
            //Tot += exp(log(mu)+(j-1)*log(mu+(f-1)*j)-j*log(f)-(mu+(f-1)*j)/f-gsl_sf_lnfact(j));
            prob *= exp((j-1)*log(mu+(f-1)*j)-(j-2)*log(mu+(f-1)*(j-1))-log(f)+(1-f)/f-log(j));
            Tot +=prob;
            j++;
        }
        result = Tot;
    }
    if (result > 1.0) result=1.0;
    return(result);
}

//COM_Poisson cdf without offset
double cdf_com_poisson_P(int x, double lambda, double nu){
    int j, k, K; //of (32) of COM-Poisson revival paper
    int MIN; //min(x,k)
    double unnormalized = 0.0;
    double normalized;
    double epsilon = 0.99;//of (32) of COM-Poisson revival paper
    double Z;
    double R;
    K = 0;
    while(lambda/pow(K+1,nu) > epsilon) K++;
    k = K + 2;
    MIN = x;
    if (k < MIN) MIN = k;
    for (j = 0; j < MIN+1; j++)
        unnormalized += exp(j*log(lambda)-nu*gsl_sf_lnfact(j));
    Z = unnormalized;
    if (x > MIN)
        for (j = MIN+1; j < x+1; j++)
            unnormalized += exp(j*log(lambda)-nu*gsl_sf_lnfact(j));
    if (k > MIN)
        for (j = MIN+1; j < k+1; j++)
            Z += exp(j*log(lambda)-nu*gsl_sf_lnfact(j));
    R = exp((k+1)*log(lambda)-nu*gsl_sf_lnfact(k+1)-log(1-epsilon)-log(Z));
    while(R > 0.000001){
        k++;
        Z += exp(k*log(lambda)-nu*gsl_sf_lnfact(k));
        R = exp((k+1)*log(lambda)-nu*gsl_sf_lnfact(k+1)-log(1-epsilon)-log(Z));
    }
    normalized = unnormalized/Z;
    if (normalized > 1.0) normalized = 1.0;
    //Rprintf("%s %i %f %f %f \n","com:",x,lambda,nu,normalized);
    return(normalized);
}

//COM_Poisson cdf with offset
double cdf_com_poisson_P2(int x, double mu, double nu){
    double lambda = pow(mu + (nu-1)/(2*nu),nu);
    int j, k, K; //of (32) of COM-Poisson revival paper
    int MIN; //min(x,k)
    double unnormalized = 0.0;
    double normalized;
    double epsilon = 0.99;//of (32) of COM-Poisson revival paper
    double Z;
    double R;
    K = 0;
    while(lambda/pow(K+1,nu) > epsilon) K++;
    k = K + 2;
    MIN = x;
    if (k < MIN) MIN = k;
    for (j = 0; j < MIN+1; j++)
        unnormalized += exp(j*log(lambda)-nu*gsl_sf_lnfact(j));
    Z = unnormalized;
    if (x > MIN)
        for (j = MIN+1; j < x+1; j++)
            unnormalized += exp(j*log(lambda)-nu*gsl_sf_lnfact(j));
    if (k > MIN)
        for (j = MIN+1; j < k+1; j++)
            Z += exp(j*log(lambda)-nu*gsl_sf_lnfact(j));
    R = exp((k+1)*log(lambda)-nu*gsl_sf_lnfact(k+1)-log(1-epsilon)-log(Z));
    while(R > 0.000001){
        k++;
        Z += exp(k*log(lambda)-nu*gsl_sf_lnfact(k));
        R = exp((k+1)*log(lambda)-nu*gsl_sf_lnfact(k+1)-log(1-epsilon)-log(Z));
    }
    normalized = unnormalized/Z;
    if (normalized > 1.0) normalized = 1.0;
    if (gsl_isinf(unnormalized)) normalized = 1.0;
    if (gsl_isnan(normalized)) normalized = 0.0;
    //Rprintf("%s %i %f %f %f %f %f \n","com:",x,lambda,nu,unnormalized,Z,normalized);
    return(normalized);
}

//Hyper Poisson cdf without mean model(Castillo & Sanchez 2013)
double cdf_hyper_poisson_P(int y, double lambda, double gamma){
    int j;
    double result = 0.0;
    for (j = 0; j < y+1; j++) result += exp(j*log(lambda) - gsl_sf_lnpoch(gamma,j)-log(gsl_sf_hyperg_1F1(1.0,gamma,lambda)));
    if (result > 1.0) result = 1.0;
    return(result);
}

//Hyper Poisson cdf with mean model(Castillo & Sanchez 2013)
struct hP_params{
    double m, g;
};

double hPmean(double lambda, void *params){
    struct hP_params *p = (struct hP_params *) params;
    double mu = p->m;
    double gamma = p->g;
    return mu - lambda + (gamma-1) * (gsl_sf_hyperg_1F1(1.0,gamma,lambda)-1)/gsl_sf_hyperg_1F1(1.0,gamma,lambda);
}

//Given mu and gamma solve for lambda
double solve_hyper_poisson(double mu, double gamma){
    double lambda;
    double r = 0;
    r = 0;
    int status;
    int iter = 0, max_iter = 100;
    struct hP_params params = {mu,gamma};
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double x_lo = MIN(mu,MAX(mu+gamma-1,gamma*mu));
    double x_hi = MAX(mu,MIN(mu+gamma-1,gamma*mu));
    gsl_function F;
    F.function = &hPmean;
    F.params = &params;
    T = gsl_root_fsolver_bisection;// brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    do{
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);        
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
        lambda = x_hi/2 + x_lo/2;
      }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    return(lambda);
}

double cdf_hyper_poisson_P2(int y, double mu, double gamma){
    double lambda;
    lambda =  solve_hyper_poisson(mu,gamma);
    int j;
    double result = 0.0;
    for (j = 0; j < y+1; j++) result += exp(j*log(lambda) - gsl_sf_lnpoch(gamma,j) - log(gsl_sf_hyperg_1F1(1.0,gamma,lambda)));
    //for (j = 0; j < y+1; j++) result += exp(j*log(lambda) - gsl_sf_lngamma(gamma+j) + gsl_sf_lngamma(gamma) -log(gsl_sf_hyperg_1F1(1.0,gamma,lambda)));
    if (result > 1.0) result = 1.0;
    return(result);
}

//Triparametric Poisson cdf without mean model (Avi et al 2004)
double cdf_tri_parametric_P(int y, double beta, double gamma, double alpha){
    gsl_sf_result lnr1, lnr2, lnr3, lnr4, lnr5, lnr6, arg1, arg2, arg3, arg4, arg5, arg6;
    int j;
    double result = 0.0;
    gsl_sf_lngamma_complex_e(gamma-alpha,-beta,&lnr1,&arg1);
    gsl_sf_lngamma_complex_e(gamma-alpha,beta,&lnr2,&arg2);
    gsl_sf_lngamma_complex_e(alpha,-beta,&lnr5,&arg5);
    gsl_sf_lngamma_complex_e(alpha,beta,&lnr6,&arg6);
    double f0 = lnr1.val + lnr2.val - gsl_sf_lngamma(gamma) - gsl_sf_lngamma(gamma-2*alpha);
    f0 = exp(f0) * gsl_sf_gamma(gamma);
    for (j = 0; j < y+1; j++){
        gsl_sf_lngamma_complex_e(alpha+j,-beta,&lnr3,&arg3);
        gsl_sf_lngamma_complex_e(alpha+j,beta,&lnr4,&arg4);
        result += exp(lnr3.val + lnr4.val - lnr5.val - lnr6.val - gsl_sf_lngamma(gamma+j) - gsl_sf_lnfact(j));
    }
    result *= f0;
    if (result > 1.0) result = 1.0;
    return(result);
}

//Triparametric Poisson cdf with mean model
double cdf_tri_parametric_P2(int y, double mu, double gamma, double alpha){
    double beta = sqrt(-alpha*alpha+mu*(gamma-2*alpha-1));
    gsl_sf_result lnr1, lnr2, lnr3, lnr4, lnr5, lnr6, arg1, arg2, arg3, arg4, arg5, arg6;
    int j;
    double result = 0.0;
    gsl_sf_lngamma_complex_e(gamma-alpha,-beta,&lnr1,&arg1);
    gsl_sf_lngamma_complex_e(gamma-alpha,beta,&lnr2,&arg2);
    gsl_sf_lngamma_complex_e(alpha,-beta,&lnr5,&arg5);
    gsl_sf_lngamma_complex_e(alpha,beta,&lnr6,&arg6);
    double f0 = lnr1.val + lnr2.val - gsl_sf_lngamma(gamma) - gsl_sf_lngamma(gamma-2*alpha);
    f0 = exp(f0) * gsl_sf_gamma(gamma);
    for (j = 0; j < y+1; j++){
        gsl_sf_lngamma_complex_e(alpha+j,-beta,&lnr3,&arg3);
        gsl_sf_lngamma_complex_e(alpha+j,beta,&lnr4,&arg4);
        result += exp(lnr3.val + lnr4.val - lnr5.val - lnr6.val - gsl_sf_lngamma(gamma+j) - gsl_sf_lnfact(j));
        //Rprintf("%s %f\n","result: ",result);
    }
    //Rprintf("%s %f\n","result: ",result);
    if (f0 > 0) result *= f0;
    else result = 0.0;
    //Rprintf("%s %f %f %f %f %f %f %f %f \n","test cpd",
    //lnr3.val,lnr4.val,lnr5.val,lnr6.val,gsl_sf_lngamma(gamma+y),gsl_sf_lnfact(y),result,f0);
    if (result > 1.0) result = 1.0;
    return(result);
}


//Univariate Normal and Binary pdf
int NormalBpdf(unsigned dim, const double *x, void *parameters, unsigned fdim, double *fval){
    double *params = parameters;
    double mu1 = params[0];
    double beta1 = params[1];
    double beta2 = params[2];
    double var1 = params[3];
    double var2 = params[4];
    double cut1 = params[5];
    double cut2 = params[6];
    double x2 = params[7];
    double m1 = params[8];
    double m2 = params[9];
    fval[0] = exp(-(x[0]-mu1)*(x[0]-mu1)/(2*var1))/sqrt(2*M_PI*var1)*
              (gsl_cdf_ugaussian_P((cut1-beta1*(x[0]-m1)-beta2*(x2-m2))/sqrt(var2))-
                                 gsl_cdf_ugaussian_P((cut2-beta1*(x[0]-m1)-beta2*(x2-m2))/sqrt(var2)));
    return 0;
}
