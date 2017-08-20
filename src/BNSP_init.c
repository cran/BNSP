#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void mvrmC(int *seed1, char **WorkingDir, int *WF1, int *sweeps1, int *burn1, int *thin1, double *Y, double *X, double *Z, int *n1, int *LG1, int *LD1, double *blockSizeProb1, int *maxBSG1, double *blockSizeProb2, int *maxBSD1, double *f1, double *g1, double *h, int *NG1, int *ND1, int *vecLG, int *vecLD, int *cusumVecLG, int *cusumVecLD, int *MVLD, double *cetaParams, double *calphaParams, double *pimu, double *pisigma, int *HN1, double *sigmaParams);
extern void OneResLtnt(int *seed1, double *X, int *Y, double *H, int *sweeps1, int *burn1, int *thin1, int *ncomp1, int *n1, int *p1, double *Vprior, double *Vdf1, double *munu, double *sigmanu,                 double *mumu, double *sigmamu, double *alphaXi, double *betaXi, double *Alphaa1, double *Alphab1, double *TruncAlpha1, double *xbar, double *xsd, double *Ymean, double *prec1, int *family1, int *sampler1, int *npred1, double *Xpred, double *Hpred, int *allEqlI, int *maxy1, double *meanReg, double *modeReg, double *Q05Reg, double *Q10Reg, double *Q15Reg, double *Q20Reg, double *Q25Reg, double *Q50Reg, double *Q75Reg, double *Q80Reg, double *Q85Reg, double *Q90Reg, double *Q95Reg, double *denReg, double *denVar, char **WorkingDir, int *WF1, int *compA);

static const R_CMethodDef CEntries[] = {
    {"mvrmC",      (DL_FUNC) &mvrmC,      32},
    {"OneResLtnt", (DL_FUNC) &OneResLtnt, 50},
    {NULL, NULL, 0}
};

void R_init_BNSP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

