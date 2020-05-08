/* matalg.h Includes functions for matrix algebra
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

//Matrix generalized inverse function with arguments (dimension, tolerance, Matrix). After this function has been called
//matrix A becomes its inverse
void ginv(int p, double tol, gsl_matrix *A){
    int i; 
    double temp, max;
    gsl_matrix *D = gsl_matrix_calloc(p,p);
    gsl_matrix *M = gsl_matrix_alloc(p,p);
    gsl_matrix *N = gsl_matrix_alloc(p,p);
    gsl_vector *eval = gsl_vector_alloc(p);
    gsl_matrix *evec = gsl_matrix_alloc(p,p);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(p);
    gsl_eigen_symmv(A,eval,evec,w);
    gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
    max = gsl_vector_get(eval,0);
    for (i=0; i < p; i++){ 
        temp = gsl_vector_get(eval,i);        
	    if (temp > tol * max) gsl_matrix_set(D,i,i,1/temp);else{gsl_matrix_set(D,i,i,0.0);}
    }
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,D,0.0,M);  //D1 = V D
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,M,evec,0.0,N); //D2 = D1 A = ginv(A)
    gsl_matrix_memcpy(A,N);
    gsl_matrix_free(D); gsl_matrix_free(M); gsl_matrix_free(N);
    gsl_vector_free(eval); gsl_matrix_free(evec); gsl_eigen_symmv_free(w);
}

//Matrix generalized inverse function with arguments (dimension, tolerance, Matrix). After this function has 
//been called matrix A becomes its inverse and det it's determinant
void ginv2(int p, double tol, gsl_matrix *A, double *det){
    int i; 
    double temp, max;
    gsl_matrix *D = gsl_matrix_calloc(p,p);
    gsl_matrix *M = gsl_matrix_alloc(p,p);
    gsl_matrix *N = gsl_matrix_alloc(p,p);
    gsl_vector *eval = gsl_vector_alloc(p);
    gsl_matrix *evec = gsl_matrix_alloc(p,p);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(p);
    gsl_eigen_symmv(A,eval,evec,w);
    gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
    max = gsl_vector_get(eval,0);
    *det = 1.0;
    for (i=0; i < p; i++){ //D = diag{1/eigen for eigen > tol and 0 otherwise}
        temp = gsl_vector_get(eval,i);
        *det *= temp;
	    if (temp > tol * max) gsl_matrix_set(D,i,i,1/temp);else{gsl_matrix_set(D,i,i,0.0);}
    }
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,D,0.0,M);  //D1 = V D
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,M,evec,0.0,N); //D2 = D1 A = ginv(A)
    gsl_matrix_memcpy(A,N);
    gsl_matrix_free(D); gsl_matrix_free(M); gsl_matrix_free(N);
    gsl_vector_free(eval); gsl_matrix_free(evec); gsl_eigen_symmv_free(w);
}

//Matrix inverse function with arguments (dimension, Matrix). After this function has been called
//matrix A becomes its inverse
void Inverse(int dim, gsl_matrix *A){
    /* Definitions needed in the LU decomposition */
    gsl_permutation *p = gsl_permutation_alloc(dim); // Permutation
    int sign;                                        // Sign of the permutation
    gsl_matrix* Ainv = gsl_matrix_alloc(dim,dim);    // Matrix to store the inverse

    gsl_linalg_LU_decomp(A,p,&sign);
    gsl_linalg_LU_invert(A,p,Ainv);
    gsl_matrix_memcpy(A,Ainv);

    gsl_matrix_free(Ainv);
    gsl_permutation_free(p);
}

//Matrix inverse root function with arguments (dimension, tolerance, Matrix). After this function has been called
//matrix A becomes A^{-1/2}
void gHalfInv(int p, double tol, gsl_matrix *A){
    int i; 
    double temp, max;
    gsl_matrix *D = gsl_matrix_calloc(p,p);
    gsl_matrix *M = gsl_matrix_alloc(p,p);
    gsl_matrix *N = gsl_matrix_alloc(p,p);
    gsl_vector *eval = gsl_vector_alloc(p);
    gsl_matrix *evec = gsl_matrix_alloc(p,p);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(p);
    gsl_eigen_symmv(A,eval,evec,w);
    gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
    max = gsl_vector_get(eval,0);
    for (i=0; i < p; i++){ //D = diag{1/eigen for eigen > tol and 0 otherwise}
        temp = gsl_vector_get(eval,i);
	    if (temp > tol * max) gsl_matrix_set(D,i,i,1/sqrt(temp));else{gsl_matrix_set(D,i,i,0.0);}
    }
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,D,0.0,M);  //D1 = V D
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,M,evec,0.0,N); //D2 = D1 A = A^{-0.5}
    gsl_matrix_memcpy(A,N);
    gsl_matrix_free(D); gsl_matrix_free(M); gsl_matrix_free(N);
    gsl_vector_free(eval); gsl_matrix_free(evec); gsl_eigen_symmv_free(w);
}

//Matrix root function with arguments (dimension, tolerance, Matrix). After this function has been called
//matrix A becomes A^{1/2}
void matHalf(int p, double tol, gsl_matrix *A){
    int i; 
    double temp, max;
    gsl_matrix *D = gsl_matrix_calloc(p,p);
    gsl_matrix *M = gsl_matrix_alloc(p,p);
    gsl_matrix *N = gsl_matrix_alloc(p,p);
    gsl_vector *eval = gsl_vector_alloc(p);
    gsl_matrix *evec = gsl_matrix_alloc(p,p);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(p);
    gsl_eigen_symmv(A,eval,evec,w);
    gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
    max = gsl_vector_get(eval,0);
    for (i=0; i < p; i++){ //D = diag{sqrt(eigen) for eigen > tol and 0 otherwise}
        temp = gsl_vector_get(eval,i);
	    if (temp > tol * max) gsl_matrix_set(D,i,i,sqrt(temp));else{gsl_matrix_set(D,i,i,0.0);}
	    //if (temp <= tol * max) puts("hi!");
    }
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,evec,D,0.0,M);  //D1 = V D
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,M,evec,0.0,N); //D2 = D1 A = ginv(A)
    gsl_matrix_memcpy(A,N);
    gsl_matrix_free(D); gsl_matrix_free(M); gsl_matrix_free(N);
    gsl_vector_free(eval); gsl_matrix_free(evec); gsl_eigen_symmv_free(w);
}

//Given an overparametrized matrix E, the function returns in D the nres nonidentified variances,
//and S becomes a matrix that is partly correlation and partly covariance matrix: E = D^{1/2} S D^{1/2}
void decomposeEtoDS(int nres, int nconf, gsl_matrix *E, gsl_matrix *D, gsl_matrix *S){
    int k,l;
    double temp1, temp2;
    gsl_matrix_memcpy(S,E);

    for (k = 0; k < nres; k++)
        gsl_matrix_set(S,k,k,1.0);

    for (k = 0; k < nres; k++){
        temp1 = gsl_matrix_get(E,k,k);
        gsl_matrix_set(D,k,k,temp1);
        for (l = 0; l < nres+nconf; l++){
            if (k != l){
                temp2 = gsl_matrix_get(S,k,l)/sqrt(temp1);
                gsl_matrix_set(S,k,l,temp2);
                gsl_matrix_set(S,l,k,temp2);
            }
        }
    }
}

//Computes |E|
double det(int p, gsl_matrix *E){
    int i;
    double detE = 1.0;
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(p);
    gsl_matrix *CopyE = gsl_matrix_alloc(p,p);
    gsl_vector *eval = gsl_vector_alloc(p);
    gsl_matrix *evec = gsl_matrix_alloc(p,p);
    gsl_matrix_memcpy(CopyE,E);
    gsl_eigen_symmv(CopyE,eval,evec,w);
    for (i=0; i < p; i++)
        detE *= gsl_vector_get(eval,i);
    gsl_eigen_symmv_free(w); gsl_matrix_free(CopyE); gsl_vector_free(eval); gsl_matrix_free(evec);
    return(detE);
}
