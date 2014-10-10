/* mathm.h Includes math functions
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

//Function that computes the parameters of the conditional distribution of y^*|w,
//where y^* is of length nres, and w of length nconf. Other inputs are the joint covariance,
//mean, realized w of length nconf, tolerance, and a vector where conditional parameters are store in this order:
//mu_1, ..., mu_nres,
//sigma_1, ..., sigma_nres,
//cor_1.2, cor_1.3, ..., cor_1.nres,
//cor_2.3, ..., cor_2.nres,
// all the way to cor_nres-1.nres
void MNCondParams(int nres, int nconf, gsl_matrix *SigmaS, gsl_vector *muS, gsl_vector *W, double tol, double *params){
    int i, j, move;
    int tot = nres + nconf;
    gsl_matrix *CopySigma = gsl_matrix_alloc(tot,tot);
    gsl_matrix *Prod1 = gsl_matrix_alloc(nres,nconf);
    gsl_vector *Diff = gsl_vector_alloc(nconf);
    gsl_vector *StoreCmean = gsl_vector_alloc(nres);
    gsl_matrix_memcpy(CopySigma,SigmaS);
    gsl_vector_memcpy(Diff,W);
    gsl_vector_view muw = gsl_vector_subvector(muS,nres,nconf);
    gsl_vector_sub(Diff,&muw.vector);
    gsl_matrix_view R = gsl_matrix_submatrix(CopySigma,0,0,nres,nres);
    gsl_matrix_view Sigma = gsl_matrix_submatrix(CopySigma,nres,nres,nconf,nconf);
    gsl_matrix_view C = gsl_matrix_submatrix(CopySigma,0,nres,nres,nconf);
    ginv(nconf,tol,&Sigma.matrix);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&C.matrix,&Sigma.matrix,0.0,Prod1);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,-1.0,Prod1,&C.matrix,1.0,&R.matrix);
    gsl_blas_dgemv(CblasNoTrans,1.0,Prod1,Diff,0.0,StoreCmean);
    gsl_vector_view muy = gsl_vector_subvector(muS,0,nres);
    gsl_vector_add(StoreCmean,&muy.vector);
    for (i = 0; i < nres; i++) params[i] = gsl_vector_get(StoreCmean,i);
    for (i = 0; i < nres; i++) params[i+nres] = sqrt(gsl_matrix_get(&R.matrix,i,i));
    move = 0;
    for (i = 0; i < nres-1; i++){
        for (j = i+1; j < nres; j++){
            params[2*nres+move] = gsl_matrix_get(&R.matrix,i,j)/(params[nres+i]*params[nres+j]);
            move++;
        }
    }
    gsl_matrix_free(CopySigma);
    gsl_matrix_free(Prod1);
    gsl_vector_free(Diff);
    gsl_vector_free(StoreCmean);
}

//Function that computes 1. the covariance matrix, and 2. the part of the mean that does not depend on the data
//nor on the mean, of the conditional distribution of y^*|w, where y^* is of length L1, and w of length L2.
//Other inputs are the joint covariance, tolerance, output mean part, output covariance part, and vector with outputs as above for cuba.
void MNCondParams1of2(int L1, int L2, gsl_matrix *JSigma, double tol, gsl_matrix *PartMean, gsl_matrix *CondCov, double *params){
    int i, j, move;
    int tot = L1 + L2;
    gsl_matrix *CopySigma = gsl_matrix_alloc(tot,tot);
    gsl_matrix_memcpy(CopySigma,JSigma);
    gsl_matrix_view R = gsl_matrix_submatrix(CopySigma,0,0,L1,L1);
    gsl_matrix_view Sigma = gsl_matrix_submatrix(CopySigma,L1,L1,L2,L2);
    gsl_matrix_view C = gsl_matrix_submatrix(CopySigma,0,L1,L1,L2);
    //ginv(L2,tol,&Sigma.matrix);
    Inverse(L2,&Sigma.matrix);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&C.matrix,&Sigma.matrix,0.0,PartMean);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,-1.0,PartMean,&C.matrix,1.0,&R.matrix);
    gsl_matrix_memcpy(CondCov,&R.matrix);
    for (i = 0; i < L1; i++)
        params[i+L1] = sqrt(gsl_matrix_get(&R.matrix,i,i));
    move = 0;
    for (i = 0; i < L1-1; i++){
        for (j = i+1; j < L1; j++){
            params[2*L1+move] = gsl_matrix_get(&R.matrix,i,j)/(params[L1+i]*params[L1+j]);
            move++;
        }
    }
    gsl_matrix_free(CopySigma);
}

//Computes E(y|w) = mu_{y} + Matrix (w-mu_{w}). y is of length L1 and w of length L2.
//Inputes are the 2 lengths, mu = E(y,w), W, Matrix, vector to store the conditional mean
void MNCondParams2of2(int L1, int L2, gsl_vector *mu, gsl_vector *W, gsl_matrix *Matrix, double *params){
    int i;
    gsl_vector *Diff = gsl_vector_alloc(L2);
    gsl_vector *StoreCmean = gsl_vector_alloc(L1);
    gsl_vector_memcpy(Diff,W);
    gsl_vector_view muw = gsl_vector_subvector(mu,L1,L2);
    gsl_vector_sub(Diff,&muw.vector);
    gsl_blas_dgemv(CblasNoTrans,1.0,Matrix,Diff,0.0,StoreCmean);
    gsl_vector_view muy = gsl_vector_subvector(mu,0,L1);
    gsl_vector_add(StoreCmean,&muy.vector);
    for (i = 0; i < L1; i++)
        params[i] = gsl_vector_get(StoreCmean,i);
    gsl_vector_free(Diff);
    gsl_vector_free(StoreCmean);
}

//Computes E(y|w) = mu_{y} + Matrix (w-mu_{w}). y is of length L1 and w of length L2.
//Inputes are the 2 lengths, mu = E(y,w), W, Matrix, vector to store the conditional mean
//mu_{y} is know to be zero
void MNCondParams2of2B(int L1, int L2, gsl_vector *mu, gsl_vector *W, gsl_matrix *Matrix, double *params){
    int i;
    gsl_vector *Diff = gsl_vector_alloc(L2);
    gsl_vector *StoreCmean = gsl_vector_alloc(L1);
    gsl_vector_memcpy(Diff,W);
    gsl_vector_view muw = gsl_vector_subvector(mu,L1,L2);
    gsl_vector_sub(Diff,&muw.vector);
    gsl_blas_dgemv(CblasNoTrans,1.0,Matrix,Diff,0.0,StoreCmean);
    for (i = 0; i < L1; i++)
        params[i] = gsl_vector_get(StoreCmean,i);
    gsl_vector_free(Diff);
    gsl_vector_free(StoreCmean);
}

//Computes E(y|w) = mu_{y} + Matrix (w-mu_{w}). y is of length L1 and w of length L2.
//Inputes are the 2 lengths, mu = E(y,w), W, Matrix, vector to store the conditional mean
//mu_{y} is know to be zero. From above there is just a difference in storing
void MNCondParams2of2C(int L1, int L2, gsl_vector *mu, gsl_vector *W, gsl_matrix *Matrix, gsl_vector *Cmean){
    gsl_vector *Diff = gsl_vector_alloc(L2);
    gsl_vector_memcpy(Diff,W);
    gsl_vector_view muw = gsl_vector_subvector(mu,L1,L2);
    gsl_vector_sub(Diff,&muw.vector);
    gsl_blas_dgemv(CblasNoTrans,1.0,Matrix,Diff,0.0,Cmean);
    gsl_vector_free(Diff);
}

//Computes E(y|w) = mu_{y} + Matrix (w-mu_{w}) when mu_{y} is know to be zero. y is of length L1 and w of length L2.
//Inputes are the 2 lengths, mu = E(w), W, Matrix, vector to store the conditional mean
void MNCondParams2of2CB(int L1, int L2, gsl_vector *mu, gsl_vector *W, gsl_matrix *Matrix, gsl_vector *Cmean){
    gsl_vector *Diff = gsl_vector_alloc(L2);
    gsl_vector_memcpy(Diff,W);
    gsl_vector_sub(Diff,mu);
    gsl_blas_dgemv(CblasNoTrans,1.0,Matrix,Diff,0.0,Cmean);
    gsl_vector_free(Diff);
}


//Computes E(y|w) = mu_{y} + Matrix (w-mu_{w}) when mu_{y} is known to be a vec of zeros.
//y is of length L1 and w of length L2. Inputes are the 2 lengths, mu = E(w), W, Matrix,
//vector to store the conditional mean
void MNCondParams2of2D(int L1, int L2, gsl_vector *mu, gsl_vector *W, gsl_matrix *Matrix, double *params){
    int i;
    gsl_vector *Diff = gsl_vector_alloc(L2);
    gsl_vector *StoreCmean = gsl_vector_alloc(L1);
    gsl_vector_memcpy(Diff,W);
    gsl_vector_sub(Diff,mu);
    gsl_blas_dgemv(CblasNoTrans,1.0,Matrix,Diff,0.0,StoreCmean);
    for (i = 0; i < L1; i++)
        params[i] = gsl_vector_get(StoreCmean,i);
    gsl_vector_free(Diff);
    gsl_vector_free(StoreCmean);
}
