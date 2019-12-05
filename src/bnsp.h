#ifndef BNSP_H_
#define BNSP_H_

#include "cubature.h"

int (*p_pcubature)(unsigned, integrand, void *, unsigned, const double *, const double *,
	                    size_t, double, double, error_norm, double *, double *);
int (*p_hcubature)(unsigned, integrand, void *, unsigned, const double *, const double *,
	                    size_t, double, double, error_norm, double *, double *);
#endif
