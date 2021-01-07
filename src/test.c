#include <R.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void t(int *p, int *N, int *V)
{
	int i;    
    
    double *BaseXg = (double*) malloc (p[0]*N[0]*V[0] * sizeof(double) );
    
    Rprintf("%i \n",p[0]*N[0]*V[0]);
    
    
    
    for (i = 0; i < 1000; i++)
        BaseXg[i] = 10;
            
    Rprintf("%s\n","Hello");
    
    free(BaseXg);
}
