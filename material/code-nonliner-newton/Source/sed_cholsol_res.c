#include "hpc.h"

index sed_cholsol_res(sed *A, double *b, index *fixed, index nfixed)
{  
	index k, p, n, *Ai;
  double *Ax, *uD;
  bool *mask;
  cs *T, *C;
	
  n = A->n; Ai = A->i, Ax = A->x;
  mask  = calloc( n , sizeof(bool) );         /* create mask */
  uD = malloc( n * sizeof(double) );  
  if ( !mask ) return(0);                     /* check memory */
  for ( k = 0; k < nfixed; k++){
    Ax[fixed[k] ] = 1.0;
    mask[ fixed[k] ]  = 1;
    uD[ fixed[k] ] = b[ fixed[k] ];
  } 
  
  /* consider entries left and above the fixed diagonal entry */
  for ( k = 0; k < n; k++)
  {
    for (p = Ai[k] ; p < Ai[k+1]; p++)
    {
      if ( mask[ Ai[p] ] )            /* row i, col j, i = Ai[p]  >  j = k */
      {
        b[ k ] -= Ax[ p ] * uD[ Ai[p] ];
        Ax[p] = 0.0;
      }     
    }
  }
  /* consider entries right and below the fixed diagonal entry */
  for ( k = 0; k < nfixed; k++) 
  {
    for (p = Ai[fixed[k]] ; p < Ai[fixed[k]+1]; p++)
    {
      b[ Ai[p]] -= Ax[ p ] * uD[ fixed[k] ];
      Ax[p] = 0.0;
    }
    b[fixed[k]] = uD[fixed[k]];
  }
  free(mask); free(uD); 
  T = cs_sed(A); C = cs_compress(T); cs_free(T); cs_dropzeros(C);
  cs_cholsol(1,C,b);                    
  cs_free(C);
  return(1);
}