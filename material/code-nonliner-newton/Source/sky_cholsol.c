#include "hpc.h"
/* solve Ax=b for x, using forward/backward substitution.
   A has L*D*L' representation and stored in skyline format */
index sky_cholsol(sky *A, double *x) 
{
  index   i, j, n, p, *Ap;
  double  *Ax, *Ad;
  
  if (!A) { printf ("(null)\n") ; return (0) ; }
  n = A->n; Ap = A->p; Ad = A->d;  Ax = A->x; 
    
  for( i = 1; i < n; i++ ){                       /* forward loop */
    for (p = Ap[i-1], j= Ap[i-1]-Ap[i]+i; p < Ap[i]; p++, j++) {  
      x[i] -= Ax[p] * x[j];
    }
  }
  for( i = 0; i < n ; i++ ) x[i] /= Ad[i];        /* diagonal scaling */
  for( i = n-1; i > 0; i-- ){                     /* backward loop */
    for (p = Ap[i-1], j= Ap[i-1]-Ap[i]+i; p < Ap[i]; p++, j++) 
    {  
      x[j] -= Ax[p] * x[i];
    }
  }
  return(1);
}
