#include "hpc.h"

#include<time.h>
#include <sys/time.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

index hpc_solve_nonlin(mesh **H, sed **A, index N, double *x,         
         double (*f_uD)(double *),  double (*f_vol)(double *, index),
          void (*f_H)(double *, index, double *), double (*f_g)(double *, index),   
         double (*f_p)(double, index),     double (*f_Dp)(double, index), 
         index solver, double tol, index maxiter) 
{
  index k, n, ncoord, nbdry, nfixed, *fixed, nedges, *bdry, iter = 0;
  double x1[2], x2[2], m[2], scal_tol, *w, *b, *Coord;
  
  n = A[N]->n; nbdry = H[N]->nbdry ; bdry = H[N]->bdry; ncoord = H[N]->ncoord ; 

//   nelem = H[N]->nelem ; nedges = H[N]->nedges; 
//   H[N]->fixed = mesh_getFixed(ncoord, bdry, nbdry, &H[N]->nfixed);
  nfixed = H[N]->nfixed ; fixed = H[N]->fixed ; 
  
  Coord = H[N]->coord;
  
  w = calloc (n, sizeof(double));       /* get workspace for sol*/
  b = calloc (n, sizeof(double));       /* get workspace for sol*/

  /* fix Dirichlet nodes */
  for ( k = 0; k < nbdry; k++)
  {
    if (!bdry[4*k+3])
    {
      x1[0] = Coord[2 * bdry[4*k]];   x1[1] = Coord[2 * bdry[4*k]+1];
      x2[0] = Coord[2 * bdry[4*k+1]]; x2[1] = Coord[2 * bdry[4*k+1]+1];
      m[0] = ( x1[0] + x2[0] ) / 2.0 ; m[1] = ( x1[1] + x2[1] ) / 2.0 ;
      w[bdry[4*k]  ]          = x[bdry[4*k]  ]          = f_uD(x1);
      w[bdry[4*k+1]]          = x[bdry[4*k+1]]          = f_uD(x2);
      w[ncoord + bdry[4*k+2]] = x[ncoord + bdry[4*k+2]] = f_uD(m);
     }
  }
  mesh_buildRhs_nonlin(H[N], b, w, f_p, f_vol, f_H, f_g);
  
    printf("aa\n");
  hpc_print_double(b, A[N]->n, "b", 1);
    printf("bb\n");

  for ( k = 0; k < nfixed; k++) b[fixed[k]] = 0.0;
  scal_tol = HPC_MAX(1e-10,hpc_ddot(b, b, n) * tol);
  printf("scal_tol = %g\n",scal_tol);
  /* loop in Newton iteration */
    while ( iter < maxiter ){

    if (!iter) TIME_SAVE(10);   
    /* compute residuum */
    mesh_buildRhs_nonlin(H[N], b, x, f_p, f_vol, f_H, f_g);  
    if (!iter) TIME_SAVE(11);
    for ( k = 0; k < nfixed; k++) b[fixed[k]] = 0.0;
    printf("tol = %g\n",hpc_ddot(b, b, n));
    if (!iter) TIME_SAVE(12);

    if (hpc_ddot(b, b, n) < scal_tol){
      
    
    printf("\n");
    printf("build rhs                     = %9i ns\n", (int) TIME_ELAPSED(10,11));
    printf("fix & ddot                    = %9i ns\n", (int) TIME_ELAPSED(11,12));
    printf("buildS_nonlin                 = %9i ns\n", (int) TIME_ELAPSED(12,13));
    printf("buildS_byInterpolation        = %9i ns\n", (int) TIME_ELAPSED(13,14));
    printf("solve lse                     = %9i ns\n", (int) TIME_ELAPSED(14,15));
    printf("========================================\n\n");

      free(w); free(b);
      return(iter);
    }
    /* assemble coefficient matrix */
    if ( !sed_buildS_nonlin(H[N], A[N], x, f_p, f_Dp) ) return(0); 
    
    if (!iter) TIME_SAVE(13);
    for (k = N; k > 0; k--){
      hpc_buildS_byInterpolation(A[k], A[k-1], H[k]->nedges, H[k]->edge2no);
    }
    if (!iter) TIME_SAVE(14);

    for ( k = 0; k < n; k++) w[k] = 0.0;
    /* solve linear system */
    switch(solver)
    {
     case 1: hpc_fmg(A, b, w, 3, H, N, 2, 2, 1); break;
     case 2: hpc_mg(A, b, w, 1e-10, 100, H, N, 2, 2, 1); break;
     default: printf("Error: Only solver 1 possible \n");
    }
    if (!iter) TIME_SAVE(15);
    
    for ( k = 0; k < n; k++) x[k] -= w[k];
    iter++;
    }
  return(-iter);
}