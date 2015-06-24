/* SpecOneRL.h Includes functions specific to one res ltnt
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

//Set Sh for step 1 of mcmc. Arguments are: # of covariates, sample size, component, # of components,
//covariates, cluster mean of covariates, imputed latent variables, cov(y^*,X), allocations, Sh
void SetShOneResLtnt(int p, int n, int h, int ncomp, double *X, double muh[ncomp][p],
           double *latentx, double nuh[ncomp][p], int *compAlloc, gsl_matrix *Sh){
    int i, k;
    gsl_matrix *StoreX = gsl_matrix_alloc(p,1);
    gsl_matrix_set_zero(Sh);
    for (i = 0; i < n; i++){
        if (compAlloc[i]==h){
            for (k = 0; k < p; k++){
                gsl_matrix_set(StoreX,k,0,(X[k*n+i]-muh[h][k]-latentx[i]*nuh[h][k]));
            }
            gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,StoreX,StoreX,1.0,Sh);
        }
    }
    gsl_matrix_free(StoreX);
}

//Set Sh for step 1 of mcmc. Arguments are: # of variables, sample size, component, # of components,
//covariates, cluster mean of variables, imputed latent variables, allocations, Sh
void SetShOneResLtntFx(int p, int n, int h, int ncomp, double *X, double muh[ncomp][p], double *latentx,
                       int *compAlloc, gsl_matrix *Sh){
    gsl_matrix_set_zero(Sh);
    int i, k;
    double baseZs[p];
    gsl_matrix_view Zs;
    for (i = 0; i < n; i++){
        if (compAlloc[i]==h){
            baseZs[0] = latentx[i] - muh[h][0];
            for (k = 0; k < (p-1); k++)
                baseZs[k+1] = X[k*n+i]-muh[h][k+1];
            Zs = gsl_matrix_view_array(baseZs,p,1);
            gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,&Zs.matrix,&Zs.matrix,1.0,Sh);
        }
    }
}

//Set sum_{ki=h} y_i*(xi-muh). Arguments are: # of covariate, sample size, component, # of comps,
//variable to be calculated, allocations, covariate design mat, cluster means, vec of latent vars.
void SetSampleTotNu(int p, int n, int h, int ncomp, double *sampleTot, int *compAlloc, double *X,
                    double muh[ncomp][p], double *latentx){
    int i, k;
    for (k = 0; k < p; k++)
        sampleTot[k] = 0.0;
    for (i = 0; i < n; i++){
        if (compAlloc[i]==h){
            for (k = 0; k < p; k++){
                sampleTot[k] += (X[k*n+i]-muh[h][k])*latentx[i];
            }
        }
    }
}

//Set sum_{ki=h} (xi-nuh*yi^*). Arguments are: # of covariate, sample size, component, # of comps,
//variable to be calculated, allocations, covariate design mat, cluster covs, vec of latent vars.
void SetSampleTotMu(int p, int n, int h, int ncomp, double *sampleTot, int *compAlloc, double *X,
                    double nuh[ncomp][p], double *latentx){
    int i, k;
    for (k = 0; k < p; k++)
        sampleTot[k] = 0.0;
    for (i = 0; i < n; i++){
        if (compAlloc[i]==h){
            for (k = 0; k < p; k++){
                sampleTot[k] += X[k*n+i] - latentx[i]*nuh[h][k];
            }
        }
    }
}

//Set sum_{ki=h} z_i^*. Arguments are: # of variables, sample size, component, # of comps,
//variable to be calculated, allocations, covariate design mat, vec of latent vars.
void SetSampleTotMuFx(int p, int n, int h, int ncomp, double *sampleTot, int *compAlloc, double *X, double *latentx){
    int i, k;
    for (k = 0; k < p; k++)
        sampleTot[k] = 0.0;
    for (i = 0; i < n; i++){
        if (compAlloc[i]==h){
            sampleTot[0] += latentx[i];
            for (k = 0; k < p-1; k++)
                sampleTot[k+1] += X[k*n+i];
        }
    }
}


void labelSwtchA(unsigned long int s, int n, int p, int ncomp, int nRespPars, double Th[ncomp][p*p], double Sigmah[ncomp][p*p],
                 double SigmahI[ncomp][p*p], double nuh[ncomp][p], double muh[ncomp][p], double Xi[ncomp][nRespPars],
                 double storenuSInu[ncomp], double storenuSI[ncomp][p], int *nmembers, int *compAlloc, double *pi){
    int i, j, h, komp, labelA, labelB;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double nonempty[ncomp];
    unsigned int vecAlloc[ncomp];
    double switchA, decisionA;

    double mnhtemp[p];
    double covtemp[p*p];
    double Xitemp[nRespPars];
    int compAlloctemp[n];
    int nmemberstemp;
    double tempnuSInu;
    double tempnuSI[p];

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
        for (j = 0; j < (p*p); j++) covtemp[j] = Th[labelA][j];
        for (j = 0; j < (p*p); j++) Th[labelA][j] = Th[labelB][j];
        for (j = 0; j < (p*p); j++) Th[labelB][j] = covtemp[j];

        for (j = 0; j < (p*p); j++) covtemp[j] = Sigmah[labelA][j];
        for (j = 0; j < (p*p); j++) Sigmah[labelA][j] = Sigmah[labelB][j];
        for (j = 0; j < (p*p); j++) Sigmah[labelB][j] = covtemp[j];

        for (j = 0; j < (p*p); j++) covtemp[j] = SigmahI[labelA][j];
        for (j = 0; j < (p*p); j++) SigmahI[labelA][j] = SigmahI[labelB][j];
        for (j = 0; j < (p*p); j++) SigmahI[labelB][j] = covtemp[j];

        for (j = 0; j < p; j++) mnhtemp[j] = nuh[labelA][j];
        for (j = 0; j < p; j++) nuh[labelA][j] = nuh[labelB][j];
        for (j = 0; j < p; j++) nuh[labelB][j] = mnhtemp[j];

        for (j = 0; j < p; j++) tempnuSI[j] = storenuSI[labelA][j];
        for (j = 0; j < p; j++) storenuSI[labelA][j] = storenuSI[labelB][j];
        for (j = 0; j < p; j++) storenuSI[labelB][j] = tempnuSI[j];

        tempnuSInu = storenuSInu[labelA];
        storenuSInu[labelA] = storenuSInu[labelB];
        storenuSInu[labelB] = tempnuSInu;

        for (j = 0; j < p; j++) mnhtemp[j] = muh[labelA][j];
        for (j = 0; j < p; j++) muh[labelA][j] = muh[labelB][j];
        for (j = 0; j < p; j++) muh[labelB][j] = mnhtemp[j];

        for (j = 0; j < nRespPars; j++) Xitemp[j] = Xi[labelA][j];
        for (j = 0; j < nRespPars; j++) Xi[labelA][j] = Xi[labelB][j];
        for (j = 0; j < nRespPars; j++) Xi[labelB][j] = Xitemp[j];

        nmemberstemp = nmembers[labelA];
        nmembers[labelA] = nmembers[labelB];
        nmembers[labelB] = nmemberstemp;

        for (i = 0; i < n; i++) compAlloctemp[i] = compAlloc[i];

        for (i = 0; i < n; i++){
            if (compAlloctemp[i] == labelA) compAlloc[i]=labelB;
            if (compAlloctemp[i] == labelB) compAlloc[i]=labelA;
        }
    }
    gsl_rng_free(r);
}

void labelSwtchB(unsigned long int s, int n, int p, int ncomp, int nRespPars, double Th[ncomp][p*p], double Sigmah[ncomp][p*p],
                 double SigmahI[ncomp][p*p], double nuh[ncomp][p], double muh[ncomp][p], double Xi[ncomp][nRespPars],
                 double storenuSInu[ncomp], double storenuSI[ncomp][p], int *nmembers, int *compAlloc, double *V){
    int i, j, h, komp, maxZ, labelC;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double equalProb, temp, switchB, decisionB;

    double mnhtemp[p];
    double covtemp[p*p];
    double Xitemp[nRespPars];
    int compAlloctemp[n];
    int nmemberstemp;
    double Vtemp;
    double tempnuSInu;
    double tempnuSI[p];

    for (h = 0; h < ncomp; h++)
        if (nmembers[h] > 0) maxZ = h;

    equalProb = (double) 1/maxZ;
    temp = gsl_ran_flat(r,0.0,1.0);

    komp=0;
    while(equalProb < temp){
        komp++;
        equalProb += (double) 1/maxZ;
    }
    labelC = komp;

    if (labelC >= (ncomp-1))
        switchB = 0.0;
    else
        switchB =  pow((1-V[labelC+1]),nmembers[labelC])/pow((1-V[labelC]),nmembers[labelC+1]);

    if (switchB > 1.0) switchB = 1.0;
    decisionB = gsl_ran_flat(r,0.0,1.0);
    //Rprintf("%s %i %i %f %f %li \n","B:",maxZ,labelC,switchB,decisionB,s);
    if (switchB > decisionB){
        for (j = 0; j < (p*p); j++) covtemp[j] = Th[labelC][j];
        for (j = 0; j < (p*p); j++) Th[labelC][j] = Th[labelC+1][j];
        for (j = 0; j < (p*p); j++) Th[labelC+1][j] = covtemp[j];

        for (j = 0; j < (p*p); j++) covtemp[j] = Sigmah[labelC][j];
        for (j = 0; j < (p*p); j++) Sigmah[labelC][j] = Sigmah[labelC+1][j];
        for (j = 0; j < (p*p); j++) Sigmah[labelC+1][j] = covtemp[j];

        for (j = 0; j < (p*p); j++) covtemp[j] = SigmahI[labelC][j];
        for (j = 0; j < (p*p); j++) SigmahI[labelC][j] = SigmahI[labelC+1][j];
        for (j = 0; j < (p*p); j++) SigmahI[labelC+1][j] = covtemp[j];

        for (j = 0; j < p; j++) mnhtemp[j] = nuh[labelC][j];
        for (j = 0; j < p; j++) nuh[labelC][j] = nuh[labelC+1][j];
        for (j = 0; j < p; j++) nuh[labelC+1][j] = mnhtemp[j];

        for (j = 0; j < p; j++) tempnuSI[j] = storenuSI[labelC][j];
        for (j = 0; j < p; j++) storenuSI[labelC][j] = storenuSI[labelC+1][j];
        for (j = 0; j < p; j++) storenuSI[labelC+1][j] = tempnuSI[j];

        tempnuSInu = storenuSInu[labelC];
        storenuSInu[labelC] = storenuSInu[labelC+1];
        storenuSInu[labelC+1] = tempnuSInu;

        for (j = 0; j < p; j++) mnhtemp[j] = muh[labelC][j];
        for (j = 0; j < p; j++) muh[labelC][j] = muh[labelC+1][j];
        for (j = 0; j < p; j++) muh[labelC+1][j] = mnhtemp[j];

        for (j = 0; j < nRespPars; j++) Xitemp[j] = Xi[labelC][j];
        for (j = 0; j < nRespPars; j++) Xi[labelC][j] = Xi[labelC+1][j];
        for (j = 0; j < nRespPars; j++) Xi[labelC+1][j] = Xitemp[j];

        nmemberstemp = nmembers[labelC];
        nmembers[labelC] = nmembers[labelC+1];
        nmembers[labelC+1] = nmemberstemp;

        for (i = 0; i < n; i++) compAlloctemp[i] = compAlloc[i];

        for (i = 0; i < n; i++){
            if (compAlloctemp[i] == labelC) compAlloc[i]=labelC+1;
            if (compAlloctemp[i] == labelC+1) compAlloc[i]=labelC;
        }

        Vtemp = V[labelC];
        V[labelC] = V[labelC+1];
        V[labelC+1] = Vtemp;
    }
    gsl_rng_free(r);
}

void labelSwtchAFx(unsigned long int s, int n, int p, int ncomp, double Th[ncomp][p*p], double Sigmah[ncomp][p*p],
                   double storeSxxI[ncomp][(p-1)*(p-1)], double storenuSI[ncomp][p-1], double storenuSInu[ncomp],
                   double muh[ncomp][p], int *nmembers, int *compAlloc, double *pi){
    int i, j, h, komp, labelA, labelB;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double nonempty[ncomp];
    unsigned int vecAlloc[ncomp];
    double switchA, decisionA;

    double mnhtemp[p];
    double covtemp[p*p];
    double covtemp2[(p-1)*(p-1)];
    int compAlloctemp[n];
    int nmemberstemp;
    double tempnuSInu;
    double tempnuSI[p-1];

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
        for (j = 0; j < (p*p); j++) covtemp[j] = Th[labelA][j];
        for (j = 0; j < (p*p); j++) Th[labelA][j] = Th[labelB][j];
        for (j = 0; j < (p*p); j++) Th[labelB][j] = covtemp[j];

        for (j = 0; j < (p*p); j++) covtemp[j] = Sigmah[labelA][j];
        for (j = 0; j < (p*p); j++) Sigmah[labelA][j] = Sigmah[labelB][j];
        for (j = 0; j < (p*p); j++) Sigmah[labelB][j] = covtemp[j];

        for (j = 0; j < (p-1)*(p-1); j++) covtemp2[j] = storeSxxI[labelA][j];
        for (j = 0; j < (p-1)*(p-1); j++) storeSxxI[labelA][j] = storeSxxI[labelB][j];
        for (j = 0; j < (p-1)*(p-1); j++) storeSxxI[labelB][j] = covtemp2[j];

        for (j = 0; j < p-1; j++) tempnuSI[j] = storenuSI[labelA][j];
        for (j = 0; j < p-1; j++) storenuSI[labelA][j] = storenuSI[labelB][j];
        for (j = 0; j < p-1; j++) storenuSI[labelB][j] = tempnuSI[j];

        tempnuSInu = storenuSInu[labelA];
        storenuSInu[labelA] = storenuSInu[labelB];
        storenuSInu[labelB] = tempnuSInu;

        for (j = 0; j < p; j++) mnhtemp[j] = muh[labelA][j];
        for (j = 0; j < p; j++) muh[labelA][j] = muh[labelB][j];
        for (j = 0; j < p; j++) muh[labelB][j] = mnhtemp[j];

        nmemberstemp = nmembers[labelA];
        nmembers[labelA] = nmembers[labelB];
        nmembers[labelB] = nmemberstemp;

        for (i = 0; i < n; i++) compAlloctemp[i] = compAlloc[i];

        for (i = 0; i < n; i++){
            if (compAlloctemp[i] == labelA) compAlloc[i]=labelB;
            if (compAlloctemp[i] == labelB) compAlloc[i]=labelA;
        }
    }
    gsl_rng_free(r);
}

void labelSwtchBFx(unsigned long int s, int n, int p, int ncomp, double Th[ncomp][p*p], double Sigmah[ncomp][p*p],
                   double storeSxxI[ncomp][(p-1)*(p-1)], double storenuSI[ncomp][p-1], double storenuSInu[ncomp],
                   double muh[ncomp][p], int *nmembers, int *compAlloc, double *V){
    int i, j, h, komp, maxZ, labelC;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,s);

    double equalProb, temp, switchB, decisionB;

    double mnhtemp[p];
    double covtemp[p*p];
    double covtemp2[(p-1)*(p-1)];
    int compAlloctemp[n];
    int nmemberstemp;
    double Vtemp;
    double tempnuSInu;
    double tempnuSI[p-1];

    for (h = 0; h < ncomp; h++)
        if (nmembers[h] > 0) maxZ = h;

    equalProb = (double) 1/maxZ;
    temp = gsl_ran_flat(r,0.0,1.0);

    komp=0;
    while(equalProb < temp){
        komp++;
        equalProb += (double) 1/maxZ;
    }
    labelC = komp;

    if (labelC >= (ncomp-1))
        switchB = 0.0;
    else
        switchB =  pow((1-V[labelC+1]),nmembers[labelC])/pow((1-V[labelC]),nmembers[labelC+1]);

    if (switchB > 1.0) switchB = 1.0;
    decisionB = gsl_ran_flat(r,0.0,1.0);

    if (switchB > decisionB){
        for (j = 0; j < (p*p); j++) covtemp[j] = Th[labelC][j];
        for (j = 0; j < (p*p); j++) Th[labelC][j] = Th[labelC+1][j];
        for (j = 0; j < (p*p); j++) Th[labelC+1][j] = covtemp[j];

        for (j = 0; j < (p*p); j++) covtemp[j] = Sigmah[labelC][j];
        for (j = 0; j < (p*p); j++) Sigmah[labelC][j] = Sigmah[labelC+1][j];
        for (j = 0; j < (p*p); j++) Sigmah[labelC+1][j] = covtemp[j];

        for (j = 0; j < (p-1)*(p-1); j++) covtemp2[j] = storeSxxI[labelC][j];
        for (j = 0; j < (p-1)*(p-1); j++) storeSxxI[labelC][j] = storeSxxI[labelC+1][j];
        for (j = 0; j < (p-1)*(p-1); j++) storeSxxI[labelC+1][j] = covtemp2[j];

        for (j = 0; j < p-1; j++) tempnuSI[j] = storenuSI[labelC][j];
        for (j = 0; j < p-1; j++) storenuSI[labelC][j] = storenuSI[labelC+1][j];
        for (j = 0; j < p-1; j++) storenuSI[labelC+1][j] = tempnuSI[j];

        tempnuSInu = storenuSInu[labelC];
        storenuSInu[labelC] = storenuSInu[labelC+1];
        storenuSInu[labelC+1] = tempnuSInu;

        for (j = 0; j < p; j++) mnhtemp[j] = muh[labelC][j];
        for (j = 0; j < p; j++) muh[labelC][j] = muh[labelC+1][j];
        for (j = 0; j < p; j++) muh[labelC+1][j] = mnhtemp[j];

        nmemberstemp = nmembers[labelC];
        nmembers[labelC] = nmembers[labelC+1];
        nmembers[labelC+1] = nmemberstemp;

        for (i = 0; i < n; i++) compAlloctemp[i] = compAlloc[i];

        for (i = 0; i < n; i++){
            if (compAlloctemp[i] == labelC) compAlloc[i]=labelC+1;
            if (compAlloctemp[i] == labelC+1) compAlloc[i]=labelC;
        }

        Vtemp = V[labelC];
        V[labelC] = V[labelC+1];
        V[labelC+1] = Vtemp;
    }
    gsl_rng_free(r);
}
