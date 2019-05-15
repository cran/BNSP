extern void ginv(int p, double tol, gsl_matrix *A);
extern void print_matrix(gsl_matrix *A);
extern void print_vector(gsl_vector *V);
extern void sampleMN(unsigned long int s, int p, gsl_vector *y, gsl_vector *mu, gsl_matrix *Sigma, double tol);
extern double logMVNormalpdf(int dim, gsl_vector *x, gsl_vector *mu, gsl_matrix *S, double tol);

//Creates the base for matrix Z Bar Start _Gamma. Inputes: 1. number of time points, 2. dimension,
//3. vec gamma, 4. All bases, 5. output of function - only vizible bases.
void setBaseZBSg(int T, int d, int *gamma, int LG, double *AllBases, double *BaseZBSg){
    int i, j, move;
    move = 0;
    for (i = 0; i < (T*d); i++){
        for (j = 0; j < (LG+1); j++){
            if (j==0) BaseZBSg[move++] = AllBases[j*T*d+i];
            if (j>0 && gamma[j-1]==1) BaseZBSg[move++] = AllBases[j*T*d+i];
		}
	}
}

//Creates the base for matrix Z Tilta_Gamma. Inputes: 1. Time points, 2. dimension, 3. vec gamma, 4. LG
//5. All bases, 6.linear predictors for variance, 7. output of function - only vizible bases.
void setBaseZtg(int T, int d, int *gamma, int Ngamma, int LG, double *AllBases, double *LPV, double *BaseZtg){
    //double BaseZBSg[(T*d)*(LG+1)];
    int i, j, t, move;
    move = 0;
    for (i = 0; i < (T*d); i++){
        for (j = 0; j < (LG+1); j++){
            //if (j==0) BaseZBSg[move++] = AllBases[j*T*d+i];
            //if (j>0 && gamma[j-1]==1) BaseZBSg[move++] = AllBases[j*T*d+i];
            if (j==0) BaseZtg[move++] = AllBases[j*T*d+i];
            if (j>0 && gamma[j-1]==1) BaseZtg[move++] = AllBases[j*T*d+i];
		}
	}
    for (t = 0; t < T; t++)
        for (i = 0; i < d; i++)
            for (j = 0; j < (Ngamma+1); j++) 
                BaseZtg[j+i*(Ngamma+1)+t*d*(Ngamma+1)] *= exp(-LPV[t]/2);                
}

//The likelihood of theta with eta integrated out, S.
double SPcalc(int T, int d, double tol, double *thetaTilde, int *gamma, int Ngamma, int LG, double ceta, double *AllBases, 
              double *LPV, double * qf2){
    double vizZt[(T*d)*(Ngamma+1)];
    double S, qf1;
    gsl_matrix *ZTZ = gsl_matrix_alloc(LG+1,LG+1);
    gsl_matrix *ZTZinv = gsl_matrix_alloc(LG+1,LG+1);
    gsl_vector *ZtTtT = gsl_vector_alloc(LG+1);
    gsl_vector *ZZIZt = gsl_vector_alloc(LG+1);
    gsl_matrix_view Ztg, ZtgTZtg, ZtgTZtgInv;
    gsl_vector_view thetaTild, ZtgTtT, ZZIZtV;
    thetaTild = gsl_vector_view_array(thetaTilde,T*d);
    gsl_blas_ddot(&thetaTild.vector,&thetaTild.vector,&qf1);
    S = qf1;
    //Rprintf("%f\n",qf1);
    setBaseZtg(T,d,gamma,Ngamma,LG,AllBases,LPV,vizZt);
	Ztg = gsl_matrix_view_array(vizZt,T*d,Ngamma+1);
	ZtgTZtg = gsl_matrix_submatrix(ZTZ,0,0,Ngamma+1,Ngamma+1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Ztg.matrix,&Ztg.matrix,0.0,&ZtgTZtg.matrix);
    ZtgTZtgInv = gsl_matrix_submatrix(ZTZinv,0,0,Ngamma+1,Ngamma+1);
	gsl_matrix_memcpy(&ZtgTZtgInv.matrix,&ZtgTZtg.matrix);
    //print_matrix(&ZtgTZtgInv.matrix);
    ginv(Ngamma+1,tol,&ZtgTZtgInv.matrix);
    //print_matrix(&ZtgTZtgInv.matrix);
    //Rprintf("%f %f %f %f %f \n",gsl_matrix_get(&ZtgTZtgInv.matrix,1,0),gsl_matrix_get(&ZtgTZtgInv.matrix,1,1),
    //                            gsl_matrix_get(&ZtgTZtgInv.matrix,1,2),gsl_matrix_get(&ZtgTZtgInv.matrix,1,3),
    //                            gsl_matrix_get(&ZtgTZtgInv.matrix,1,4));
    
    ZtgTtT = gsl_vector_subvector(ZtTtT,0,Ngamma+1);
    gsl_blas_dgemv(CblasTrans,1.0,&Ztg.matrix,&thetaTild.vector,0.0,&ZtgTtT.vector);
	ZZIZtV = gsl_vector_subvector(ZZIZt,0,Ngamma+1);
    gsl_blas_dgemv(CblasNoTrans,1.0,&ZtgTZtgInv.matrix,&ZtgTtT.vector,0.0,&ZZIZtV.vector);
    gsl_blas_ddot(&ZtgTtT.vector,&ZZIZtV.vector,qf2);
    S -= (ceta/(1+ceta))*(*qf2);
    gsl_matrix_free(ZTZ); gsl_matrix_free(ZTZinv); gsl_vector_free(ZtTtT); gsl_vector_free(ZZIZt);
    return(S);
}

//returns proposed values for the vector of binary indicators.
//Inputs: 1. random number generator seed, 2. vec of current indicators, 3. length of vec of ind,
//4. block number, 5. block size, 6. shuffled indeces, 7. & 8. prior parameters, 9. store prop vec ind
void proposeBlockInd(unsigned long int s, int *vecInd, int L, int B, int BS, int *shufInd,
                     double c, double d, int *vecIndP){
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);
    int i, NS, NBC, TBS; //NS(block), NS(complement of block), true block size
    double denom;
    TBS = BS;
    if ((B+1)*TBS > L) TBS = L-B*BS;
    double priorProbs[TBS+1];
    unsigned int vecNS[TBS+1];
    int proposal[TBS];
    //puts("in prop block 1");
    NBC = 0;
    for (i = 0; i < (B*BS); i++)
        NBC += vecInd[shufInd[i]];
    for (i = ((B+1)*BS); i < L; i++)
        NBC += vecInd[shufInd[i]];    
    denom = gsl_sf_beta((double) (NBC+c),(double) (L-TBS-NBC+d));
    for (i = 0; i < (TBS+1); i++)
        priorProbs[i] = gsl_sf_choose(TBS,i) *
            gsl_sf_beta((double) (NBC+i+c),(double) (L-NBC-i+d))/denom;
    gsl_ran_multinomial(r,TBS+1,1,priorProbs,vecNS);
    NS = 0;
    while(vecNS[NS]==0) NS++;
    for (i = 0; i < TBS; i++) proposal[i] = 0;
    for (i = 0; i < NS; i++) proposal[i] = 1;
    gsl_ran_shuffle(r,proposal,TBS,sizeof (int));

    for (i = 0; i < L; i++)
        vecIndP[i] = vecInd[i];

    for (i = 0; i < TBS; i++)
        vecIndP[shufInd[B*BS+i]] = proposal[i];
        
    //for (i = 0; i < (TBS+1); i++)        
    //    Rprintf("%i %i %i %f %i %i %i \n",B,TBS,NBC,priorProbs[i],vecNS[i],NS,shufInd[B*BS+i]);
    
    //for (i = 0; i < TBS; i++)        
    //    Rprintf("%i %i %i \n",B,TBS,proposal[i]);
       
    gsl_rng_free(r);
}

//Posterior inverse covariance matrix of theta.
void iCovPostTheta(int T, int d, double tol, int *gamma, int Ngamma, int LG, double ceta, double sigma2,
		  double tau2, double *AllBases, double *LPV, gsl_matrix *A){
    int r, c, DIM;
    DIM = Ngamma+1;
    //if (DIM == 0) DIM = 1;
    double vizZt[(T*d)*(Ngamma+1)];
    double expandLPV[T*d];
    gsl_matrix_view Z;
    gsl_matrix *ZTZ = gsl_matrix_alloc(DIM,DIM);
    gsl_matrix *ZTZinv = gsl_matrix_alloc(DIM,DIM);
    gsl_matrix *ZZTZinv = gsl_matrix_alloc(T*d,DIM);
    gsl_matrix *ZZTZinvZT = gsl_matrix_alloc(T*d,T*d);
    gsl_matrix *I = gsl_matrix_alloc(T*d,T*d);
    gsl_matrix_set_identity(I);
    gsl_matrix_set_identity(A);
    gsl_matrix_scale(A,1/tau2);
    //if (Ngamma > 0){
    setBaseZtg(T,d,gamma,Ngamma,LG,AllBases,LPV,vizZt);
    Z = gsl_matrix_view_array(vizZt,T*d,DIM);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Z.matrix,&Z.matrix,0.0,ZTZ);
    gsl_matrix_memcpy(ZTZinv,ZTZ);
    ginv(DIM,tol,ZTZinv);        
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&Z.matrix,ZTZinv,0.0,ZZTZinv);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,ZZTZinv,&Z.matrix,0.0,ZZTZinvZT);	
    //} //if Ngamma = 0 then ZZTZinvZT = 0;
    gsl_matrix_scale(ZZTZinvZT,ceta/(1+ceta));
    gsl_matrix_sub(I,ZZTZinvZT);
    for (r = 0; r < T; r++)
        for (c = 0; c < d; c++)
            expandLPV[c+r*d] = LPV[r];
    for (r = 0; r < (T*d); r++)
        for (c = 0; c < (T*d); c++)
            I->data[r * I->tda + c] *= exp(-expandLPV[r]/2)*exp(-expandLPV[c]/2);
    gsl_matrix_scale(I,1/sigma2);
    //ginv(T*d,tol,I);
    gsl_matrix_add(A,I);
    ginv(T*d,tol,A);    
    gsl_matrix_free(ZTZ); gsl_matrix_free(ZTZinv); gsl_matrix_free(ZZTZinv);
    gsl_matrix_free(ZZTZinvZT); gsl_matrix_free(I);
}

//Compute posterior mean of eta. Arguments: 1. dimension, tolerance, Matrix
void postMeanVarEta(int T, int d, double tol, int *gamma, int Ngamma, int LG, double sigma2, double ceta, 
		    double *LPV, double *AllBases, double *thetaTilde, gsl_vector *MeanEta, gsl_matrix *varEta, int sw){
    double vizZt[(T*d)*(Ngamma+1)];
    gsl_matrix_view Z;
    gsl_vector_view thetaTild;
    thetaTild = gsl_vector_view_array(thetaTilde,T*d);
    gsl_matrix *ZTZ = gsl_matrix_alloc(Ngamma+1,Ngamma+1);
    gsl_matrix *ZTZinv = gsl_matrix_alloc(Ngamma+1,Ngamma+1);
    gsl_matrix *ZTZinvZT = gsl_matrix_alloc(Ngamma+1,T*d);
    //if (Ngamma > 0){
    setBaseZtg(T,d,gamma,Ngamma,LG,AllBases,LPV,vizZt);
    Z = gsl_matrix_view_array(vizZt,T*d,Ngamma+1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Z.matrix,&Z.matrix,0.0,ZTZ);
    gsl_matrix_memcpy(ZTZinv,ZTZ);
    ginv(Ngamma+1,tol,ZTZinv);
    
    gsl_matrix_memcpy(varEta,ZTZinv);
    gsl_matrix_scale(varEta,sigma2*ceta/(1+ceta));
    
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,ZTZinv,&Z.matrix,0.0,ZTZinvZT);
    gsl_blas_dgemv(CblasNoTrans,1.0,ZTZinvZT,&thetaTild.vector,0.0,MeanEta);
    gsl_vector_scale(MeanEta,ceta/(1+ceta));
    
    //} else gsl_vector_set_all(MeanEta,0.0);
    //if (sw==0) {Rprintf("%f %f %f %f %f \n",LPV[0],LPV[1],LPV[2],LPV[3],LPV[4]);}
    //if (sw==0) print_matrix(&Z.matrix);
    //if (sw==0) print_matrix(ZTZ);
    //if (sw==0) print_matrix(ZTZinv);
    //if (sw==0) print_matrix(ZTZinvZT);
    //if (sw==0) {Rprintf("%f ",gsl_vector_get(MeanEta,0));Rprintf(" %f ",gsl_vector_get(MeanEta,1));Rprintf(" %f ",gsl_vector_get(MeanEta,2));
      //          Rprintf(" %f ",gsl_vector_get(MeanEta,3));Rprintf(" %f ",gsl_vector_get(MeanEta,4));Rprintf(" %f \n",gsl_vector_get(MeanEta,5));}
    gsl_matrix_free(ZTZ); gsl_matrix_free(ZTZinv); gsl_matrix_free(ZTZinvZT);
}

//Compute vector of squared residuals
void cSqRes(int T, int d, int *gamma, int Ngamma, int LG, double *AllBases, gsl_vector *MeanEta, double *theta, double *sqRes){
    int i;
    double BaseZBSg[T*d*(Ngamma+1)];
    setBaseZBSg(T,d,gamma,LG,AllBases,BaseZBSg);
    gsl_vector *thetaHat = gsl_vector_alloc(T*d);
    gsl_matrix_view Z;
	Z = gsl_matrix_view_array(BaseZBSg,T*d,Ngamma+1);
    gsl_blas_dgemv(CblasNoTrans,1.0,&Z.matrix,MeanEta,0.0,thetaHat);
    for (i = 0; i < (T*d); i++){
        sqRes[i] = pow(theta[i] - thetaHat->data[i*thetaHat->stride],2);
        //if (i<5) Rprintf("%i %f %f %f \n",i,theta[i],thetaHat->data[i*thetaHat->stride],sqRes[i]);
    }    
    gsl_vector_free(thetaHat);
}

//Compute alpha^hat_delta and Delta_delta
void DeltaAlphaHat(int T, int d, double tol, double *LPV, double *sqRes, int *delta, int Ndelta, int start, int end,
        double *AllBases, double sigma2, double *sigma2t, double calpha, gsl_matrix *D, gsl_vector *alphaHat){
    double baseZBd[T*d*Ndelta];
    int t, i, j, move;
    double vecd[T*d];
    gsl_matrix *I = gsl_matrix_alloc(Ndelta,Ndelta);
    gsl_matrix *V = gsl_matrix_alloc(Ndelta,T*d);
    gsl_matrix_set_identity(I);
    gsl_matrix_view Z;
    gsl_vector_view smallD;
    for (t = 0; t < T; t++)
        for (i = 0; i < d; i++) 
            vecd[i+t*d] = log(sigma2) + LPV[t] + (sqRes[i+t*d] - sigma2t[t])/sigma2t[t];
            //vecd[i+t*d] = LPV[t] + (sqRes[i+t*d] - sigma2t[t])/sigma2t[t];            
            //vecd[i+t*d] = log(sigma2t[t]) + (sqRes[i+t*d] - sigma2t[t])/sigma2t[t];
    move = 0;
    for (i = 0; i < (T*d); i++)
        for (j = start; j < end; j++)
            if (delta[j]==1) baseZBd[move++] = AllBases[T*d+j*T*d+i];    
    Z = gsl_matrix_view_array(baseZBd,T*d,Ndelta);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Z.matrix,&Z.matrix,0.0,D);
    //gsl_matrix_scale(D,0.5);
    gsl_matrix_scale(I,1/calpha);
    gsl_matrix_add(D,I);
    //Inverse(Ndelta,D);
    ginv(Ndelta,tol,D);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,D,&Z.matrix,0.0,V);
    smallD = gsl_vector_view_array(vecd,T*d);
    gsl_blas_dgemv(CblasNoTrans,1.0,V,&smallD.vector,0.0,alphaHat);
    gsl_matrix_free(I); gsl_matrix_free(V);
}


//Clustering on the correlations

//Creates the base for matrix Z Bar Star_{Gamma_h}. 
//Inputes: 1. Time points 2. dimension 3. number of surfaces 4. surface 5. LG 6. matrix of gammas 7. compAlloc 
//8. All bases 9. output of function - only vizible bases.
void setBaseZBSgh(int T, int d, int H, int h, int LG, int gamma[H][LG], int *compAlloc,   
                  double *AllBases, double *BaseZBSgh){
    int i, j, t, move;
    move = 0;
    for (t = 0; t < T; t++){
        for (i = 0; i < d; i++){
		    if (compAlloc[i]==h){				
                for (j = 0; j < (LG+1); j++){
                    if (j==0) BaseZBSgh[move++] = AllBases[j*T*d+i+t*d];
                    if (j>0 && gamma[h][j-1]==1) BaseZBSgh[move++] = AllBases[j*T*d+i+t*d];                                          
		        }
		    }		    
	    }
	}
}
     
//Creates the base for matrix Z Tilta_Gamma_h. 
//Inputes: 1. Time points 2. dimension 3. number of surfaces 4. surface 5. LG 6. matrix of gammas 7. compAlloc 
//8. nmembers 9. All bases 10.linear predictors for variance 11. output of function - only vizible bases.
void setBaseZtgh(int T, int d, int H, int h, int LG, int gamma[H][LG], int Ngamma, int *compAlloc, int nmembers,  
                 double *AllBases, double *LPV, double *BaseZtgh){
    //double BaseZBSg[(T*nmembers)*(Ngamma+1)];
    int i, j, t, move;
    move = 0;
    for (t = 0; t < T; t++){
        for (i = 0; i < d; i++){
		    if (compAlloc[i]==h){				
                for (j = 0; j < (LG+1); j++){
                    //if (j==0) BaseZBSg[move++] = AllBases[j*T*d+i+t*d];
                    //if (j>0 && gamma[h][j-1]==1) BaseZBSg[move++] = AllBases[j*T*d+i+t*d];                                          
                    if (j==0) BaseZtgh[move++] = AllBases[j*T*d+i+t*d];
                    if (j>0 && gamma[h][j-1]==1) BaseZtgh[move++] = AllBases[j*T*d+i+t*d];                                          
		        }
		    }		    
	    }
	}
    for (t = 0; t < T; t++)
        for (i = 0; i < nmembers; i++)
            for (j = 0; j < (Ngamma+1); j++)
                BaseZtgh[j+i*(Ngamma+1)+t*nmembers*(Ngamma+1)] *= exp(-LPV[t]/2);//*BaseZBSg[j+i*(Ngamma+1)+t*nmembers*(Ngamma+1)];
}

//The likelihood of theta_h with eta_h integrated out, Sh.                
double SPh(int T, int d, int H, int h, double tol, double *thetaTilde, int LG, int gamma[H][LG], int Ngamma,
           int *compAlloc, int nmembers, double ceta, double *AllBases, double *LPV, double *qf2){
    double vizZt[(T*nmembers)*(Ngamma+1)];
    double vizThetaTild[T*nmembers];
    double S, qf1;
    int i, t, move;
    gsl_matrix *ZTZ = gsl_matrix_alloc(LG+1,LG+1);
    gsl_matrix *ZTZinv = gsl_matrix_alloc(LG+1,LG+1);
    gsl_vector *ZtTtT = gsl_vector_alloc(LG+1);
    gsl_vector *ZZIZt = gsl_vector_alloc(LG+1);
    gsl_matrix_view Ztg, ZtgTZtg, ZtgTZtgInv;
    gsl_vector_view thetaTild, ZtgTtT, ZZIZtV;    
    move = 0; 
    for (t = 0; t < T; t++)
        for (i = 0; i < d; i++)
            if (compAlloc[i] == h) vizThetaTild[move++] = thetaTilde[i+t*d];     
    thetaTild = gsl_vector_view_array(vizThetaTild,T*nmembers);
    gsl_blas_ddot(&thetaTild.vector,&thetaTild.vector,&qf1);
    S = qf1;
    setBaseZtgh(T,d,H,h,LG,gamma,Ngamma,compAlloc,nmembers,AllBases,LPV,vizZt);	
	Ztg = gsl_matrix_view_array(vizZt,T*nmembers,Ngamma+1);
	ZtgTZtg = gsl_matrix_submatrix(ZTZ,0,0,Ngamma+1,Ngamma+1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Ztg.matrix,&Ztg.matrix,0.0,&ZtgTZtg.matrix);
    ZtgTZtgInv = gsl_matrix_submatrix(ZTZinv,0,0,Ngamma+1,Ngamma+1);
	gsl_matrix_memcpy(&ZtgTZtgInv.matrix,&ZtgTZtg.matrix);
    ginv(Ngamma+1,tol,&ZtgTZtgInv.matrix);
    ZtgTtT = gsl_vector_subvector(ZtTtT,0,Ngamma+1);
    gsl_blas_dgemv(CblasTrans,1.0,&Ztg.matrix,&thetaTild.vector,0.0,&ZtgTtT.vector);
	ZZIZtV = gsl_vector_subvector(ZZIZt,0,Ngamma+1);
    gsl_blas_dgemv(CblasNoTrans,1.0,&ZtgTZtgInv.matrix,&ZtgTtT.vector,0.0,&ZZIZtV.vector);
    gsl_blas_ddot(&ZtgTtT.vector,&ZZIZtV.vector,qf2);
    S -= (ceta/(1+ceta))*(*qf2);
    gsl_matrix_free(ZTZ); gsl_matrix_free(ZTZinv); gsl_vector_free(ZtTtT); gsl_vector_free(ZZIZt);
    return(S);
}

//Posterior inverse covariance matrix of theta.
void iCovPostThetah(int T, int d, int H, int h, double tol, int LG, int gamma[H][LG], 
           int *compAlloc, int nmembers, int *Ngamma, double ceta, double sigma2,
		   double tau2, double *AllBases, double *LPV, gsl_matrix *A){
    int r, c, DIM;
    DIM = Ngamma[h]+1;
    double vizZt[(T*nmembers)*(Ngamma[h]+1)];
    double expandLPV[T*nmembers];
    gsl_matrix_view Z;
    gsl_matrix *ZTZ = gsl_matrix_alloc(DIM,DIM);
    gsl_matrix *ZTZinv = gsl_matrix_alloc(DIM,DIM);
    gsl_matrix *ZZTZinv = gsl_matrix_alloc(T*nmembers,DIM);
    gsl_matrix *ZZTZinvZT = gsl_matrix_alloc(T*nmembers,T*nmembers);
    gsl_matrix *I = gsl_matrix_alloc(T*nmembers,T*nmembers);
    gsl_matrix_set_identity(I);
    gsl_matrix_set_identity(A);
    gsl_matrix_scale(A,1/tau2);
    setBaseZtgh(T,d,H,h,LG,gamma,Ngamma[h],compAlloc,nmembers,AllBases,LPV,vizZt);
    Z = gsl_matrix_view_array(vizZt,T*nmembers,DIM);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Z.matrix,&Z.matrix,0.0,ZTZ);
    gsl_matrix_memcpy(ZTZinv,ZTZ);
    ginv(DIM,tol,ZTZinv);        
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&Z.matrix,ZTZinv,0.0,ZZTZinv);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,ZZTZinv,&Z.matrix,0.0,ZZTZinvZT);	
    gsl_matrix_scale(ZZTZinvZT,ceta/(1+ceta));
    gsl_matrix_sub(I,ZZTZinvZT);
    for (r = 0; r < T; r++)
        for (c = 0; c < nmembers; c++)
            expandLPV[c+r*nmembers] = LPV[r];
    for (r = 0; r < (T*nmembers); r++)
        for (c = 0; c < (T*nmembers); c++)
            I->data[r * I->tda + c] *= exp(-expandLPV[r]/2)*exp(-expandLPV[c]/2);
    gsl_matrix_scale(I,1/sigma2);
    gsl_matrix_add(A,I);
    ginv(T*nmembers,tol,A);    
    gsl_matrix_free(ZTZ); gsl_matrix_free(ZTZinv); gsl_matrix_free(ZZTZinv);
    gsl_matrix_free(ZZTZinvZT); gsl_matrix_free(I);
}

//Posterior inverse covariance matrix of theta_kl.
void iCovPostThetahkl(int T, int d, int H, int h, double tol, int LG, int gamma[H][LG], 
           int *compAlloc, int nmembers, int *Ngamma, double ceta, double sigma2,
		   double tau2, double *AllBases, double *LPV, gsl_matrix *A){
    int r, c, DIM;
    DIM = Ngamma[h]+1;      
    double vizZt2[T*(Ngamma[h]+1)];
    int compAlloc2[d];    
    gsl_matrix *ZTZ = gsl_matrix_alloc(DIM,DIM);
    gsl_matrix *ZTZ2 = gsl_matrix_alloc(DIM,DIM);        
    gsl_matrix *ZTZinv = gsl_matrix_alloc(DIM,DIM);  
    gsl_matrix *ZZTZinv = gsl_matrix_alloc(T,DIM);
    gsl_matrix *ZZTZinvZT = gsl_matrix_alloc(T,T);
    gsl_matrix *I = gsl_matrix_alloc(T,T);
    gsl_matrix_set_identity(I);
    gsl_matrix_view Z, Z2;
    if (nmembers>0){
        double vizZt[(T*nmembers)*(Ngamma[h]+1)];
        setBaseZtgh(T,d,H,h,LG,gamma,Ngamma[h],compAlloc,nmembers,AllBases,LPV,vizZt);
        Z = gsl_matrix_view_array(vizZt,T*nmembers,DIM);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Z.matrix,&Z.matrix,0.0,ZTZ);
    }
    else
        gsl_matrix_set_identity(ZTZ);     
    gsl_matrix_scale(ZTZ,1/ceta);
    for (r = 0; r < d; r++)
        compAlloc2[r] = H+2;
    compAlloc2[0] = h;    
    setBaseZtgh(T,d,H,h,LG,gamma,Ngamma[h],compAlloc2,1,AllBases,LPV,vizZt2);
    Z2 = gsl_matrix_view_array(vizZt2,T,DIM);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Z2.matrix,&Z2.matrix,0.0,ZTZ2);
    gsl_matrix_add(ZTZ,ZTZ2);
    gsl_matrix_memcpy(ZTZinv,ZTZ);
    ginv(DIM,tol,ZTZinv);      
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&Z2.matrix,ZTZinv,0.0,ZZTZinv);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,ZZTZinv,&Z2.matrix,0.0,ZZTZinvZT);	
    gsl_matrix_sub(I,ZZTZinvZT);
    for (r = 0; r < T; r++)
        for (c = 0; c < T; c++)
            I->data[r * I->tda + c] *= exp(-LPV[r]/2)*exp(-LPV[c]/2);
    gsl_matrix_scale(I,1/sigma2);
    gsl_matrix_memcpy(A,I); 
    gsl_matrix_free(ZTZ); gsl_matrix_free(ZTZ2); gsl_matrix_free(ZTZinv); gsl_matrix_free(ZZTZinv);
    gsl_matrix_free(ZZTZinvZT); gsl_matrix_free(I);
}

//Compute posterior mean of eta_h.
void postMeanVarEtaH(int T, int d, int H, int h, double tol, int LG, int gamma[H][LG], int *compAlloc, int *nmembers, 
            int *Ngamma, double sigma2, double ceta, double *LPV, double *AllBases, double *thetaTilde, 
            gsl_vector *MeanEta, gsl_matrix *varEta, int sw, int kk){
    
    double vizZt[(T*nmembers[h])*(Ngamma[h]+1)];
    gsl_matrix_view Z;
    
    gsl_vector_view thetaTild;
    double vizThetaTilda[T*nmembers[h]];
    int i, t, move;    
    move = 0; 
    for (t = 0; t < T; t++)
        for (i = 0; i < d; i++)
            if (compAlloc[i] == h) vizThetaTilda[move++] = thetaTilde[t*d+i];     
    thetaTild = gsl_vector_view_array(vizThetaTilda,T*nmembers[h]);
    
    gsl_matrix *ZTZ = gsl_matrix_alloc(Ngamma[h]+1,Ngamma[h]+1);
    gsl_matrix *ZTZinv = gsl_matrix_alloc(Ngamma[h]+1,Ngamma[h]+1);
    gsl_matrix *ZTZinvZT = gsl_matrix_alloc(Ngamma[h]+1,T*nmembers[h]);    
    setBaseZtgh(T,d,H,h,LG,gamma,Ngamma[h],compAlloc,nmembers[h],AllBases,LPV,vizZt);    
    Z = gsl_matrix_view_array(vizZt,T*nmembers[h],Ngamma[h]+1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Z.matrix,&Z.matrix,0.0,ZTZ);
    gsl_matrix_memcpy(ZTZinv,ZTZ);
    ginv(Ngamma[h]+1,tol,ZTZinv);   
    gsl_matrix_memcpy(varEta,ZTZinv);
    gsl_matrix_scale(varEta,sigma2*ceta/(1+ceta));    
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,ZTZinv,&Z.matrix,0.0,ZTZinvZT);
    gsl_blas_dgemv(CblasNoTrans,1.0,ZTZinvZT,&thetaTild.vector,0.0,MeanEta);
    gsl_vector_scale(MeanEta,ceta/(1+ceta));
     
    if (sw>9999 && kk==-1) {Rprintf("%f %f %f %f %f %f %f \n",LPV[0],LPV[1],LPV[2],LPV[3],LPV[4],sigma2,ceta);}
    if (sw>9999 && kk==-1) print_matrix(&Z.matrix);
    if (sw>9999 && kk==-1) print_matrix(ZTZ);
    if (sw>9999 && kk==-1) print_matrix(ZTZinv);
    if (sw>9999 && kk==-1) print_matrix(ZTZinvZT);
    if (sw>9999 && kk==-1) print_matrix(varEta);
    if (sw>9999 && kk==-1) print_vector(&thetaTild.vector);
    if (sw>9999 && kk==-1){ 
        for (i = 0; i < (Ngamma[h]+1); i++) Rprintf("%f ", gsl_vector_get(MeanEta,i));
        Rprintf("\n ");
    }
    gsl_matrix_free(ZTZ); gsl_matrix_free(ZTZinv); gsl_matrix_free(ZTZinvZT);
}

//Compute vector of squared residuals for cluster h
void cSqResh(int T, int d, int h, int *compAlloc, int *gamma, int *Ngamma, int LG, double *AllBases, 
             gsl_vector *MeanEta, double *theta, double *sqRes){
    int i, t;
    double BaseZBSg[T*d*(Ngamma[h]+1)];
    setBaseZBSg(T,d,gamma,LG,AllBases,BaseZBSg);
    gsl_vector *thetaHat = gsl_vector_alloc(T*d);
    gsl_matrix_view Z;
	Z = gsl_matrix_view_array(BaseZBSg,T*d,Ngamma[h]+1);
    gsl_blas_dgemv(CblasNoTrans,1.0,&Z.matrix,MeanEta,0.0,thetaHat);
    for (t = 0; t < T; t++)
        for (i = 0; i < d; i++)
            if (compAlloc[i] == h) 
                sqRes[i+t*d] = pow(theta[i+t*d] - thetaHat->data[(i+t*d)*thetaHat->stride],2);   
    gsl_vector_free(thetaHat);
}

//Clustering on the variables

void compAllocVtoCompAlloc(int G, int p, int *compAllocV, int * compAlloc){
    //int d = 10;
    int i, j, k, l, temp, move, move2;
    //gsl_matrix *R = gsl_matrix_alloc(p,p);
    //gsl_matrix_set_all(R,-99);
    //for (k = 0; k < d; k++) compAlloc[k] = -99;
    move = 0;
    for (k = 0; k < G; k++){
        temp = 0;
        for (i = 0; i < (p-1); i++){
            for (j = (i+1); j < p; j++){
                if (compAllocV[i] == compAllocV[j] && compAllocV[i] == k){ 
                    //gsl_matrix_set(R,i,j,move);
                    compAlloc[i*p+j-(i+1)*(i+2)/2] = move;
                    //Rprintf("%s %i %i %i %i %i \n","i*p+j-(i+1)",i,j,p,i*p+j-(i+1),move);
                    temp++;
		        }
            }
        }
        if (temp > 0) move++;
    }      
    move2=0;
    for (k = 0; k < (G-1); k++){
        for (l = (k+1); l < G; l++){
            temp = 0;
            for (i = 0; i < (p-1); i++){
                for (j = (i+1); j < p; j++){                 
			        if (((compAllocV[i] == k && compAllocV[j] == l) || (compAllocV[i] == l && compAllocV[j] == k))){ 
                            //gsl_matrix_set(R,i,j,move+move2);
                            //compAlloc[j-1-i+(p*(p-1)-(p-i)*(p-1-i))/2] = move+move2;
                            compAlloc[i*p+j-(i+1)*(i+2)/2] = move+move2;
                            //Rprintf("%s %i %i %i %i %i \n","j-1-i+(p*(p-1)-(p-i)*(p-1-i))/2",i,j,p,j-1-i+(p*(p-1)-(p-i)*(p-1-i))/2,move+move2);
                            temp++;		            
	                }
                }
            }
            if (temp > 0) move2++;
        }
    }
    //print_matrix(R);
    //for (k = 0; k < d; k++) 
    //    Rprintf("%i ",compAlloc[k]);
    //Rprintf("\n");
    
    //gsl_matrix_free(R);
}

//Fisher transformation
double FisherTr(double r, int I){
	double comp;
	if (I==0) comp = r;
	if (I==1) comp = 0.5*log((1+r)/(1-r));
	return(comp);
}

// St^* = sum_{i=1}^n_t y_{it} y_{it}^T. Inputs are: vec of all responses sorted by id and time, 
// time vec, length of Y, time of interest, dim of Yit, store St^* 
void computeStStar(double *Y, int *time, int N, int t, int p, gsl_matrix *StStar){
    int i;
    int c = 0;
    double Yit[p];
    gsl_matrix_view YitVec;    
    gsl_matrix_set_zero(StStar);  
    for (i = 0; i < N; i++){       
        if (time[i]==t) Yit[c++] = Y[i];
        if (c == p){
		    c=0;
		    YitVec = gsl_matrix_view_array(Yit,p,1);
            gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,&YitVec.matrix,&YitVec.matrix,1.0,StStar);      
		}         
    }    
}

//Sets up X_{i,gamma} tilde transpose columnwise. Inputs are: 1. dim of y-vec, 2. number of sampling units, 3. sampling unit,
// 4. length of gamma, 5. N(gamma) total, 6. LPV, 7. X matrix: one row per sampling unit, 8. gamma matrix, 9.vec to store
void setXigammaStarT(int p, int m, int i, int LG, int Ngamma, double sigma2ij[m][p], double *X, int gamma[p][LG], double *base){
    int j, k, move;    
    move = 0;    
    for (k = 0; k < p; k++)
        for (j = 0; j < (LG+1); j++)
            if ((j==0) || (j>0 && gamma[k][j-1]==1)) base[k*(Ngamma+p)+move++] = X[i*(LG+1)+j]/sqrt(sigma2ij[i][k]);
}

// Function S for multivariate responses: dim of response, # sampling units, length of gammas per regression, tol, ceta,
// total Ngamma, vec of all Y-tilde, LPV, X = [X1,X2,...Xm], gamma mat, R^(-1), Sstar.  
double ScalcMult(int p, int m, int LG, double tol, double ceta, int Ngamma, double *Ytilde, double sigma2ij[m][p], double *X, 
                 int gamma[p][LG], gsl_matrix *Ri, gsl_matrix *St, double *qf2){
    int i, k;
    double S;
    double trace = 0.0;
    double Yi[p];
    double base[p*(Ngamma+p)];                                      
    for (k = 0; k < (p*(Ngamma+p)); k++) 
        base[k] = 0;
    gsl_matrix *XRi = gsl_matrix_alloc(Ngamma+p,p); 
    gsl_matrix *XRiX = gsl_matrix_calloc(Ngamma+p,Ngamma+p);
    gsl_matrix *RiS = gsl_matrix_alloc(p,p);
    gsl_vector *XRiY = gsl_vector_calloc(Ngamma+p);
    gsl_vector *AB = gsl_vector_alloc(Ngamma+p);                    
    gsl_matrix_view Xig;
    gsl_vector_view YiVec;                        
    for (i = 0; i < m; i++){
	    for (k = 0; k < p; k++)        
            Yi[k] = Ytilde[i*p+k];
    	YiVec = gsl_vector_view_array(Yi,p);		
        setXigammaStarT(p,m,i,LG,Ngamma,sigma2ij,X,gamma,base);
	    Xig = gsl_matrix_view_array(base,p,(Ngamma+p));	           
        if (0==1){  
            Rprintf("%s %i \n", "base:", i);
            for (k = 0; k < (p*(Ngamma+p)); k++)
                Rprintf(" %f ", base[k]);                
			//print_matrix(&Xig.matrix);
		}
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Xig.matrix,Ri,0.0,XRi);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,XRi,&Xig.matrix,1.0,XRiX);
        gsl_blas_dgemv(CblasNoTrans,1.0,XRi,&YiVec.vector,1.0,XRiY);
	}
	//print_matrix(XRiX);
	ginv(Ngamma+p,tol,XRiX);
	//print_matrix(XRiX);
	//Rprintf("%f %f %f %f %f \n",gsl_matrix_get(XRiX,1,0)/0.014,gsl_matrix_get(XRiX,1,1)/0.014,
    //                            gsl_matrix_get(XRiX,1,2)/0.014,gsl_matrix_get(XRiX,1,3)/0.014,
    //                            gsl_matrix_get(XRiX,1,4)/0.014);
	
	//Rprintf("%f %f %f %f %f \n",gsl_matrix_get(XRiX,1,0)*0.014,gsl_matrix_get(XRiX,1,1)*0.014,
    //                            gsl_matrix_get(XRiX,1,2)*0.014,gsl_matrix_get(XRiX,1,3)*0.014,
    //                            gsl_matrix_get(XRiX,1,4)*0.014);
	
		        
	gsl_blas_dgemv(CblasTrans,1.0,XRiX,XRiY,0.0,AB);
	
	gsl_blas_ddot(XRiY,AB,qf2);
	S = -(ceta/(1+ceta))*(*qf2);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Ri,St,0.0,RiS);	            
    for (i = 0; i < p; i++) 
        trace += gsl_matrix_get(RiS,i,i);
    //print_matrix(Ri);print_matrix(St);
    //Rprintf("%f \n",trace);		            
    S += trace;		            
    gsl_matrix_free(XRi); gsl_matrix_free(XRiX); gsl_matrix_free(RiS); 
    gsl_vector_free(XRiY); gsl_vector_free(AB); 
    //Rprintf("%s %f %f %f %f \n","from func: ",S,ceta,trace,qf2);
    return(S);
}

// Function that computes tr(R^{-1} sum_{i=1}^m (Ytil_i-Xtil_i^* beta^*)(Ytil_i-Xtil_i^* beta^*)^T) for multivariate responses
double NormalQuadr(int p, int m, int LG, int Ngamma, double *Ytilde, double sigma2ij[m][p], double *X, 
                   int gamma[p][LG], gsl_matrix *Ri, double *beta){        
    int i, k;
    double trace = 0.0;
    double Yi[p];
    double base[p*(Ngamma+p)];                                      
    for (k = 0; k < (p*(Ngamma+p)); k++) base[k] = 0;
    gsl_matrix *St = gsl_matrix_calloc(p,p);     
    gsl_matrix *RiS = gsl_matrix_alloc(p,p);
    gsl_matrix *XBeta = gsl_matrix_alloc(p,1);                    
    gsl_matrix_view YiVec, Xig, Beta;                            
    Beta = gsl_matrix_view_array(beta,(Ngamma+p),1);        
    for (i = 0; i < m; i++){
	    for (k = 0; k < p; k++)        
            Yi[k] = Ytilde[i*p+k];
    	YiVec = gsl_matrix_view_array(Yi,p,1);		
        setXigammaStarT(p,m,i,LG,Ngamma,sigma2ij,X,gamma,base);
	    Xig = gsl_matrix_view_array(base,p,(Ngamma+p));	         	              
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&Xig.matrix,&Beta.matrix,0.0,XBeta);
		gsl_matrix_sub(&YiVec.matrix,XBeta);		
		gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,&YiVec.matrix,&YiVec.matrix,1.0,St);        
	}
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Ri,St,0.0,RiS);	            
    for (i = 0; i < p; i++) 
        trace += gsl_matrix_get(RiS,i,i);        
    gsl_matrix_free(St); gsl_matrix_free(RiS); gsl_matrix_free(XBeta); 
    return(trace);
}

// Function that computes the posterior mean and variance of \ubeta: dim of response, # sampling units, length of gammas per regression, tol, ceta,
// total Ngamma, vec of all Y-tilde, LPV, X = [X1,X2,...Xm], gamma mat, R^(-1), mean, var.  
void postMeanVarEta2(int p, int m, int LG, double tol, double ceta, int Ngamma, double *Ytilde, double sigma2ij[m][p], double *X, 
                       int gamma[p][LG], gsl_matrix *Ri, gsl_vector *MeanEta, gsl_matrix *varEta){
    int i, k;
    double Yi[p];
    double base[p*(Ngamma+p)];
    for (k = 0; k < (p*(Ngamma+p)); k++) 
        base[k] = 0;                                                     
    gsl_matrix *XRi = gsl_matrix_alloc(Ngamma+p,p); 
    gsl_matrix *XRiX = gsl_matrix_calloc(Ngamma+p,Ngamma+p);
    gsl_vector *XRiY = gsl_vector_calloc(Ngamma+p);                        
    gsl_matrix_view Xig;
    gsl_vector_view YiVec;                        
    for (i = 0; i < m; i++){	    
	    for (k = 0; k < p; k++)        
            Yi[k] = Ytilde[i*p+k];
    	YiVec = gsl_vector_view_array(Yi,p);		
        setXigammaStarT(p,m,i,LG,Ngamma,sigma2ij,X,gamma,base);
	    Xig = gsl_matrix_view_array(base,p,Ngamma+p);	           
	    
	    //if (i<2) print_matrix(&Xig.matrix);
	    
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Xig.matrix,Ri,0.0,XRi);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,XRi,&Xig.matrix,1.0,XRiX);
        gsl_blas_dgemv(CblasNoTrans,1.0,XRi,&YiVec.vector,1.0,XRiY);
	}
	ginv(Ngamma+p,tol,XRiX);	        
    gsl_matrix_memcpy(varEta,XRiX);
    gsl_matrix_scale(varEta,ceta/(1+ceta));
    gsl_blas_dgemv(CblasNoTrans,1.0,varEta,XRiY,0.0,MeanEta);
    gsl_matrix_free(XRi); gsl_matrix_free(XRiX); gsl_vector_free(XRiY);
}

//Sets up X_{gamma} transpose. Inputs are: dim of y-vec, number of sampling units, length of gamma, 
// N(gamma) total, X matrix: on row per sampling unit, gamma matrix, vec to store
void setBaseXg(int p, int m, int LG, int Ngamma, double *X, int gamma[p][LG], double *base){
    int i, j, k, move;        
    for (i = 0; i < m; i++){
        move = 0;
        for (k = 0; k < p; k++)
            for (j = 0; j < (LG+1); j++)
                if ((j==0) || (j>0 && gamma[k][j-1]==1)) base[i*p*(Ngamma+p)+k*(Ngamma+p)+move++] = X[i*(LG+1)+j];   
	}
}

//Compute vector of squared residuals
void cSqRes2(int p, int m, int LG, int gamma[p][LG], int Ngamma, double *X, gsl_vector *MeanEta, double *Y, double *sqRes){
    int i;
    double BaseXg[m*p*(Ngamma+p)];    
    for (i = 0; i < (m*p*(Ngamma+p)); i++) 
        BaseXg[i] = 0;                                     
    gsl_vector *yHat = gsl_vector_alloc(p*m);
    gsl_matrix_view Z;
    setBaseXg(p,m,LG,Ngamma,X,gamma,BaseXg);
	Z = gsl_matrix_view_array(BaseXg,m*p,Ngamma+p);
	
	//for (i = 0; i < (Ngamma+p); i++)
      //  Rprintf("%s %f \n",gsl_matrix_get(&Z.matrix,,)); 	
	
    gsl_blas_dgemv(CblasNoTrans,1.0,&Z.matrix,MeanEta,0.0,yHat);
    for (i = 0; i < (p*m); i++){
        sqRes[i] = pow(Y[i] - yHat->data[i*yHat->stride],2);
        //if (i<5) Rprintf("%i %f %f %f \n",i,theta[i],thetaHat->data[i*thetaHat->stride],sqRes[i]);
    }
    gsl_vector_free(yHat);
}

//Compute alpha^hat_delta and Delta_delta
void DeltaAlphaHatExp(int m, int p, int l, double tol, double LPV[m][p], double *sqRes, int *delta, int Ndelta, 
    int start, int end, double *AllBases, double sigma2, double sigma2ij[m][p], double calpha, gsl_matrix *D, gsl_vector *alphaHat){    
    double base[m*Ndelta];
    int i, j, move;
    double vecd[m];
    gsl_matrix *I = gsl_matrix_alloc(Ndelta,Ndelta);
    gsl_matrix *V = gsl_matrix_alloc(Ndelta,m);
    gsl_matrix_set_identity(I);
    gsl_matrix_view Z;
    gsl_vector_view smallD;    
    for (i = 0; i < m; i++){                         
        vecd[i] = log(sigma2) + LPV[i][l] + (sqRes[i*p+l] - sigma2ij[i][l])/sigma2ij[i][l];
	    //Rprintf("%s %f %f %f %f \n","vec d: ",vecd[i],log(sigma2),LPV[i][l],sqRes[i*p+l]);
	}
    move = 0;    
    for (i = 0; i < m; i++)
        for (j = start; j < end; j++)
            if (delta[j-start]==1) base[move++] = AllBases[m+j*m+i];    
    Z = gsl_matrix_view_array(base,m,Ndelta);
    //print_matrix(&Z.matrix);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&Z.matrix,&Z.matrix,0.0,D);
    gsl_matrix_scale(I,1/calpha);
    gsl_matrix_add(D,I);
    ginv(Ndelta,tol,D);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,D,&Z.matrix,0.0,V);
    smallD = gsl_vector_view_array(vecd,m);
    gsl_blas_dgemv(CblasNoTrans,1.0,V,&smallD.vector,0.0,alphaHat);
    gsl_matrix_free(I); gsl_matrix_free(V);
}
