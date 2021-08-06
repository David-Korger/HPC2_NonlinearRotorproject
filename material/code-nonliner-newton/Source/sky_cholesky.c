#include "hpc.h"
/* perform Cholesky decomposition on positive definite matrix stored 
 * in skyline format, the original matrix destroyed */

index sky_cholesky(sky *A) 
{
  index   i, j, k, n, p, *Ap, pij, pki,pkj, jmin, kmin, tmax;
  double  x, *Ax, *Ad, *w, tmp;
  
  if (!A) { printf ("(null)\n") ; return (0) ; }
  n = A->n; Ap = A->p; Ad = A->d;  Ax = A->x; 

  for( i = 1; i < n; i++ ){               /* loop over all rows */
    jmin = i + (Ap[i-1] - Ap[i]);         /* first col index in i'th row */
    tmax = HPC_MAX(1,jmin);
    for (j = tmax, pij=Ap[i]-i+tmax; j < i; j++, pij++) {        
      tmax = HPC_MAX(jmin,j + (Ap[j-1] - Ap[j]));
      for (k = tmax, pki=Ap[i]-i+tmax, pkj=Ap[j]-j+tmax;k<j;k++,pki++,pkj++)
      { 
        Ax[pij] -= Ax[pki] * Ax[pkj];
      } 
    }
    for (j = jmin, pij=Ap[i]-i+jmin; j < i; j++,pij++) {   
      tmp = Ax[pij];  
      Ax[pij] /= Ad[j];                    /* scale i'th row */
      Ad[i] -= Ax[pij] * tmp ;             /* update diagonal entry */
    }    
  }
  return(1);
}
