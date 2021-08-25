#include "hpc.h"

index hpc_broyden(sed **A, double *b, double *x, double tol, index maxit,
             mesh **H)          
{
    /* Idea behind broydens method:
       Don't calculate (numerical) jacobian matrix DF(x_n) at every iter. 
       step, instead only calculate DF(x_0) at the beginning, and then 
       make small corrections to DF(x_0) at each iter. step */ 
    
    /* Calculate initial Jacobi */
    sed **DF;
    sed_calc_jacobi(A, x, DF);
    
    /* iterate until convergence */
    for ( j = 0 ; j <= maxit; j++) 
    {



    }

	return (maxit);
}