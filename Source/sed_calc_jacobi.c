#include "hpc.h"

index sed_calc_jacobi (sed **A, double *x, sed *DF)
/* Calculate numerical Jacobian of A, at point x, and save to DF */
{
    double delta = 1.0e-7;

    double tmp* = calloc(A->n);
    /* H-Method: DF = ( A*(x + delta) - A*(x - delta) )/(2delta) */

}
