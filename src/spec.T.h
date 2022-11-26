extern void ginv(int p, double tol, gsl_matrix *A);

//Compute psi^hat_ksi and Delta
void DeltaPsiHat(int m, int p, int l, int LC, double tol,  
                 double omegaij[m][p], double U[m][p], double psi[p][LC], double cpsi, 
                 int start, int end, double *AllBases, int *delta, int Ndelta, 
                 gsl_vector *alphaHat, gsl_matrix *D){
    int i, j, move;
	double vecu[m];
	double vecp[m];
    double base[m*Ndelta];
    double base2[m*Ndelta];
    double base3[Ndelta];
    gsl_matrix *I = gsl_matrix_alloc(Ndelta,Ndelta);
    gsl_vector *V = gsl_vector_alloc(Ndelta);
    gsl_matrix_set_identity(I);
    gsl_matrix_view W, PW;
    gsl_vector_view vecU, vecPsi;          
    for (i = 0; i < m; i++){
        //Rprintf("%i %i %f \n",i,l,omegaij[i][l]);
        //Rprintf("%i %i %f \n",i,l,gsl_sf_psi(omegaij[i][l]/2));
        //Rprintf("%i %i %f \n",i,l,gsl_sf_psi_1(omegaij[i][l]/2));
        vecu[i] = omegaij[i][l] * (-gsl_sf_psi(omegaij[i][l]/2) + log(omegaij[i][l]/2) + 1 + log(U[i][l]) - U[i][l]) / 2;
	}
    for (i = 0; i < m; i++)
        vecp[i] = omegaij[i][l]*omegaij[i][l]*(gsl_sf_psi_1(omegaij[i][l]/2)/2 - 1/omegaij[i][l])/2;
    move = 0;    
    for (i = 0; i < m; i++){
        for (j = start; j < end; j++){
            if (delta[j-start]==1){
				base[move] = AllBases[m+j*m+i]; 
				base2[move] = vecp[i]*base[move];
				move++;
			}
		}
	}
    move = 0;
    for (j = start; j < end; j++)
        if (delta[j-start]==1) 
           base3[move++] = psi[l][j];       
    W = gsl_matrix_view_array(base,m,Ndelta);
    PW = gsl_matrix_view_array(base2,m,Ndelta);
    vecPsi = gsl_vector_view_array(base3,Ndelta);
    vecU = gsl_vector_view_array(vecu,m);
    //print_matrix(&W.matrix);            
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&W.matrix,&PW.matrix,0.0,D);    
    gsl_blas_dgemv(CblasTrans,1.0,&W.matrix,&vecU.vector,0.0,V);
    gsl_blas_dgemv(CblasNoTrans,1.0,D,&vecPsi.vector,1.0,V);
    gsl_matrix_scale(I,1/cpsi);
    gsl_matrix_add(D,I);
    ginv(Ndelta,tol,D);    
    gsl_blas_dgemv(CblasNoTrans,1.0,D,V,0.0,alphaHat);
    gsl_matrix_free(I); gsl_vector_free(V);
}
