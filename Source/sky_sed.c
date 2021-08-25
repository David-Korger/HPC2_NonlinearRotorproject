#include "hpc.h"
/* S = sym sky matrix of a sym extrated diagonal matrix T */
sky *sky_sed (const sed *T)
{
  index k, p, n, *Sp, *Ti;
  double *Sx, *Sd, *Tx;
  sky *S ;
  
  if (!T) return (NULL) ;                        /* check inputs */
  n = T->n; Tx = T->x ; Ti = T->i ; 
  S = sky_alloc(n, 0) ;                          /* allocate result */
  if (!S) sky_free(S) ;                          /* check inputs */
  Sp = S->p; Sd = S->d;
  
  for (k = 0; k < n; k++) Sp[k] = k;             /* get envelope */
  
  for (k = 0; k < n; k++)
  {
    for (p = Ti[k]; p < Ti[k+1]; p++)  Sp[ Ti[p] ] = HPC_MIN( Sp[ Ti[p] ], k); 
  } 
  for (k = 1; k < n; k++) Sp[k] = Sp[k-1] + ( k - Sp[k] );
  S->x = Sx = calloc (Sp[n-1], sizeof (double)) ;
  if ( !Sx ) sky_free(S) ;                       /* out of memory */
  for (k = 0; k < n; k++) Sd[k] = Tx[k];
  for (k = 0; k < n; k++)                        /* put data into place, */
  {
    for (p = Ti[k]; p < Ti[k+1]; p++)            /* store strict lower part */ 
      Sx[ Sp[Ti[p]] - Ti[p] + k ] += Tx[p];
  }
  return (S) ;
}
