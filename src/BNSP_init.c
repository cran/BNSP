#include <stdlib.h>
#include <Rconfig.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "bnsp.h"

/* .C calls */
extern void mult(int *seed1, char **WorkingDir, int *WF1, int *sweeps1, int *burn1, int *thin1, int *m1, int *p1, double *Y, double *X, double *Z, int *LG1, int *LD1, double *blockSizeProb1, int *maxBSG1, double *blockSizeProb2, int *maxBSD1, int *NG1, int *ND1, int *vecLG, int *vecLD, int *cusumVecLG, int *cusumVecLD, int *MVLD1, double *tuneSigma2R, double *tuneCa, double *tuneSigma2, double *tuneCb, double *h, double *tuneR, double *pimu, double *cetaParams, double *pisigma, int *HNca, double *calphaParams, double *Rparams, int *HNszk, double *szkParams, double *tau1, int *FT, double *dev, int *cont, int *LASTgamma, int *LASTdelta, double *LASTalpha, double *LASTsigma2zk, double *LASTR, double *LASTmuR, double *LASTsigma2R, double *LASTceta, double *LASTcalpha, int *LASTsw, double *LASTDE, int *LASTWB);
extern void multg(int *seed1, char **WorkingDir, int *WF1, int *sweeps1, int *burn1, int *thin1, int *m1, int *p1, double *Y, double *X, double *Z, int *LG1, int *LD1, double *blockSizeProb1, int *maxBSG1, double *blockSizeProb2, int *maxBSD1, int *NG1, int *ND1, int *vecLG, int *vecLD, int *cusumVecLG, int *cusumVecLD, int *MVLD1, double *tuneSigma2R, double *tuneCa, double *tuneSigma2, double *tuneCb, double *h, double *tuneR, double *pimu, double *cetaParams, double *pisigma, int *HNca, double *calphaParams, double *Rparams, int *HNszk, double *szkParams, double *tau1, int *FT, double *dev, int *H1, double *DPparams,int *cont, int *LASTgamma, int *LASTdelta, double *LASTalpha, double *LASTsigma2zk, double *LASTR, double *LASTmuR, double *LASTsigma2R, double *LASTceta, double *LASTcalpha, int *LASTcompAlloc, double *LASTalphaDP, int *LASTsw, double *LASTDE, int *LASTWB);
extern void multgv(int *seed1, char **WorkingDir, int *WF1, int *sweeps1, int *burn1, int *thin1, int *m1, int *p1, double *Y, double *X, double *Z, int *LG1, int *LD1, double *blockSizeProb1, int *maxBSG1, double *blockSizeProb2, int *maxBSD1, int *NG1, int *ND1, int *vecLG, int *vecLD, int *cusumVecLG, int *cusumVecLD, int *MVLD1, double *tuneSigma2R, double *tuneCa, double *tuneSigma2, double *tuneCb, double *h, double *tuneR, double *pimu, double *cetaParams, double *pisigma, int *HNca, double *calphaParams, double *Rparams, int *HNszk, double *szkParams, double *tau1, int *FT, double *dev, int *G1, double *DPparams, int *cont, int *LASTgamma, int *LASTdelta, double *LASTalpha, double *LASTsigma2zk, double *LASTR, double *LASTmuR, double *LASTsigma2R, double *LASTceta, double *LASTcalpha, int *LASTcompAllocV, double *LASTalphaDP, int *LASTsw, double *LASTDE, int *LASTWB);
extern void mvrmC(int *seed1, char **WorkingDir, int *WF1, int *sweeps1, int *burn1, int *thin1, double *Y, double *X, double *Z, int *n1, int *LG1, int *LD1, double *blockSizeProb1, int *maxBSG1, double *blockSizeProb2, int *maxBSD1, double *f1, double *g1, double *h, double *fca1, int *NG1, int *ND1, int *vecLG, int *vecLD, int *cusumVecLG, int *cusumVecLD, int *MVLD, double *cetaParams, int *HNca1, double *calphaParams, double *pimu, double *pisigma, int *HNsg1, double *sigmaParams, double *dev, int *cont, int *LASTgamma, int *LASTdelta, double *LASTalpha, double *LASTsigma2zk, double *LASTceta, double *LASTcalpha, int *LASTWB);
extern void OneResLtnt(int *seed1, double *X, int *Y, double *H, int *sweeps1, int *burn1, int *thin1, int *ncomp1, int *n1, int *p1, int *NBC1, double *Vprior, double *eta1, double *mumu, double *sigmamu, double *alphaXi, double *betaXi, double *Alphaa1, double *Alphab1, double *TruncAlpha1, double *xbar, double *xsd, double *Ymean, int *family1, int *sampler1, int *npred1, double *Xpred, double *Hpred, int *maxy1, int *c, double *meanReg, double *medianReg, double *q1Reg, double *q3Reg, double *modeReg, double *denReg, double *denVar, char **WorkingDir, int *WF1);

static const R_CMethodDef CEntries[] = {    
    {"mult",       (DL_FUNC) &mult,       54},
    {"multg",       (DL_FUNC) &multg,     58},
    {"multgv",       (DL_FUNC) &multgv,   58},
    {"mvrmC",      (DL_FUNC) &mvrmC,      43},
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
