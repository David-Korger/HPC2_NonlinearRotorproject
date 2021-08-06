#include "hpc.h"

index hpc_mg(sed **A, double *b, double *x, double tol, index maxit,
             mesh **H, index nLevel, index pre, index post, index gamma)          
{
    
  index j, k, nIter, *p;
  double *wx, *wb, *wr, **xx, **bb, **rr, tmp, scal_tol ;
  
  if( !nLevel ) /* If no hierachy -> solve exactly */
  { 
    for ( k = 0; k < H[0]->nfixed; k++) b[H[0]->fixed[k]] = x[H[0]->fixed[k]];
    sed_cholsol_res(A[0], b, H[0]->fixed, H[0]->nfixed);
    for ( k = 0; k < A[0]->n ; k++) x[k] = b[k];
	}
  else
  {
    p = malloc ((nLevel+2) * sizeof(index)); 
    p[0] = 0;
    for (j = 0; j <=nLevel; j++) p[j+1] = p[j] + A[j]->n;
  
    wx = malloc( p[nLevel]   * sizeof(double) );
    wb = malloc( p[nLevel]   * sizeof(double) );
    wr = malloc( p[nLevel+1] * sizeof(double) );
  
    xx = malloc( (nLevel+1) * sizeof(double*) );
    bb = malloc( (nLevel+1) * sizeof(double*) );
    rr = malloc( (nLevel+1) * sizeof(double*) );
    
    for (j = 0; j < nLevel; j++){
      xx[j] = wx + p[j]; bb[j] = wb + p[j]; rr[j] = wr + p[j];
    }
    xx[nLevel] = x; bb[nLevel] = b; rr[nLevel] = wr + p[nLevel];
    free(p); 
    
    /* set tolerance, scal_tol =  || b(freenodes) ||^2 * tol */
    for ( k = 0 ; k < A[nLevel]->n ; k++) rr[nLevel][k] = bb[nLevel][k] ;
    for ( k = 0 ; k < H[nLevel]->nfixed; k++) 
                                   rr[nLevel][H[nLevel]->fixed[k]] = 0.0;      
    tmp = 0.0;
    for ( k=0 ; k < A[nLevel]->n ; k++) tmp += rr[nLevel][k] * rr[nLevel][k];
    scal_tol = tol * tmp;
    printf("scal_tol = %g\n",tmp);          
    /* iterate until convergence */
    for ( j = 0 ; j <= maxit; j++) 
    {
      /* Compute residual, r = b - A * x */
      for ( k = 0 ; k < A[nLevel]->n ; k++) rr[nLevel][k] = 0.0;
      sed_gaxpy(A[nLevel], xx[nLevel], rr[nLevel]);
      for ( k = 0 ; k < A[nLevel]->n ; k++) 
                                 rr[nLevel][k] = bb[nLevel][k] - rr[nLevel][k];
      /* Consider constrains */
      for ( k = 0 ; k < H[nLevel]->nfixed; k++) 
            rr[nLevel][H[nLevel]->fixed[k]] = 0.0;  
      tmp = 0.0;
      for ( k=0 ; k < A[nLevel]->n ; k++) tmp += rr[nLevel][k] * rr[nLevel][k];
      printf("Norm(resid)^2 = %g\n",tmp);         
      if ( tmp <= scal_tol ) 
      {
        free(wx); free(wb); free(wr);
        free(xx); free(bb); free(rr);
        printf("# Iter = %g\n",(double) j);         
        return (j);
      }
      hpc_mg_cycle(A, H, nLevel, bb, xx, rr, pre, post, gamma);
    }
    fprintf(stderr,
        "Max iterations reached: maxit = %ld, r = %g\n", maxit, sqrt(tmp));
    nIter = j;
    free(wx); free(wb); free(wr);
    free(xx); free(bb); free(rr);
	}
	return (nIter);
}