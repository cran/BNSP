#include <stdlib.h>
#include <Rconfig.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "bnsp.h"

/* .C calls */
extern void mvrmC(int *seed1, char **WorkingDir, int *WF1, int *sweeps1, int *burn1, int *thin1, double *Y, double *X, double *Z, int *n1, int *LG1, int *LD1, double *blockSizeProb1, int *maxBSG1, double *blockSizeProb2, int *maxBSD1, double *f1, double *g1, double *h, int *NG1, int *ND1, int *vecLG, int *vecLD, int *cusumVecLG, int *cusumVecLD, int *MVLD, double *cetaParams, double *calphaParams, double *pimu, double *pisigma, int *HN1, double *sigmaParams);
extern void OneResLtnt(int *seed1, double *X, int *Y, double *H, int *sweeps1, int *burn1, int *thin1, int *ncomp1, int *n1, int *p1, int *NBC1, double *Vprior, double *eta1, double *mumu, double *sigmamu, double *alphaXi, double *betaXi, double *Alphaa1, double *Alphab1, double *TruncAlpha1, double *xbar, double *xsd, double *Ymean, int *family1, int *sampler1, int *npred1, double *Xpred, double *Hpred, int *maxy1, int *c, double *meanReg, double *medianReg, double *q1Reg, double *q3Reg, double *modeReg, double *denReg, double *denVar, char **WorkingDir, int *WF1);

static const R_CMethodDef CEntries[] = {
    {"mvrmC",      (DL_FUNC) &mvrmC,      33},
    {"OneResLtnt", (DL_FUNC) &OneResLtnt, 39},
    {NULL, NULL, 0}
};

void R_init_BNSP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE); 
    
	                    
    p_pcubature = (int (*) (unsigned, integrand, void *, unsigned, const double *, const double *,
	                    size_t, double, double, error_norm, double *, double *)) R_GetCCallable("cubature","pcubature");
    p_hcubature = (int (*) (unsigned, integrand, void *, unsigned, const double *, const double *,
	                    size_t, double, double, error_norm, double *, double *)) R_GetCCallable("cubature","hcubature");
}
















