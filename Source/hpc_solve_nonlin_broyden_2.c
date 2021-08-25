#include "hpc.h"

#include<time.h>
#include <sys/time.h>

struct timeval tv_solve_broyden_2[50];
#define TIME_SAVE(j)   (gettimeofday(&tv_solve_broyden_2[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv_solve_broyden_2[k].tv_sec - tv_solve_broyden_2[j].tv_sec)+(tv_solve_broyden_2[k].tv_usec - tv_solve_broyden_2[j].tv_usec))

index hpc_solve_nonlin_broyden_2(mesh **H, sed **B_0, index N, double *x,         
         double (*f_uD)(double *),  double (*f_vol)(double *, index),
         void (*f_H)(double *, index, double *), double (*f_g)(double *, index),   
         double (*f_p)(double, index),     double (*f_Dp)(double, index), 
         index solver, double tol, index maxiter, double damping) 
{
  index k, j, n, ncoord, nbdry, nfixed, *fixed, nedges, *bdry, iter = 0;
  double x1[2], x2[2], m[2], scal_tol, *y, *p, *p_k, *b_old, *b, *Coord;

  // create memory for all n  nxn-diagonal-sparseness-matrices:
  index **S = (index **)malloc(n * sizeof(index*)); 
  for (j = 0; j < n; j++) {
    S[j] = (index *)calloc (n, sizeof(index));       
  }
  
  n = B_0[N]->n; nbdry = H[N]->nbdry ; bdry = H[N]->bdry; ncoord = H[N]->ncoord ; 

//   nelem = H[N]->nelem ; nedges = H[N]->nedges; 
//   H[N]->fixed = mesh_getFixed(ncoord, bdry, nbdry, &H[N]->nfixed);
  nfixed = H[N]->nfixed ; fixed = H[N]->fixed ; 
  
  Coord = H[N]->coord;
  
  y = calloc (n, sizeof(double));      
  b = calloc (n, sizeof(double));       
  p = calloc (n, sizeof(double));       
  p_k = calloc (n, sizeof(double));       
  b_old = calloc (n, sizeof(double));       

  /* fix Dirichlet nodes */
  for ( k = 0; k < nbdry; k++)
  {
    if (!bdry[4*k+3])
    {
      x1[0] = Coord[2 * bdry[4*k]];   x1[1] = Coord[2 * bdry[4*k]+1];
      x2[0] = Coord[2 * bdry[4*k+1]]; x2[1] = Coord[2 * bdry[4*k+1]+1];
      m[0] = ( x1[0] + x2[0] ) / 2.0 ; m[1] = ( x1[1] + x2[1] ) / 2.0 ;
      y[bdry[4*k]  ]          = x[bdry[4*k]  ]          = f_uD(x1);
      y[bdry[4*k+1]]          = x[bdry[4*k+1]]          = f_uD(x2);
      y[ncoord + bdry[4*k+2]] = x[ncoord + bdry[4*k+2]] = f_uD(m);
     }
  }
  mesh_buildRhs_nonlin(H[N], b, y, f_p, f_vol, f_H, f_g);
  
  hpc_print_double(b, B_0[N]->n, "b", 1);
  
  for ( k = 0; k < nfixed; k++) b[fixed[k]] = 0.0;
  scal_tol = HPC_MAX(1e-10,hpc_ddot(b, b, n) * tol);
  printf("scal_tol = %10g\n",scal_tol);
  
  /* start broyden iteration */

  /* 1. compute first broyden matrix B_0 = jacobian of x_0 */
  /* jacobian of x_0 = stiffness matrix (I think) */
  if ( !sed_buildS_nonlin(H[N], B_0[N], x, f_p, f_Dp) ){
    printf("ERROR in computing S at start of Broyden iteration !\n");
    return(0); 
  } 

  /* 1b save sparseness pattern (S is zero when B is zero, else S is 1) */
  for(k = 0; k < n; k++) {if (B_0[N]->x[k] != 0) { S[k][k] = 1; }}

  for(k = 0; k < n; k++) {
    for(j = B_0[N]->i[k]; j < B_0[N]->i[k+1]; j++){
      /* When an index is listed within SED-Format, then it must be nonzero */
      /* so set S to one at every listed point in the matrix */
      S[k][j] = 1;
    }
  }

  /* 2. Solve B_0 p = -F(x_0),  for p */
  /* b <- F(x_0) */
  mesh_buildRhs_nonlin(H[N], b, x, f_p, f_vol, f_H, f_g);  
  /* empty solution vector p */
  for (k = 0; k < n; k++) p[k] = 0.0;
  /* p <- solve B_0 p = -F(x_0) */
  switch(solver)
  {
   case 1: hpc_fmg(B_0, b, p, 3, H, N, 2, 2, 1); break;
   case 2: hpc_mg(B_0, b, p, 1e-10, 100, H, N, 2, 2, 1); break;
   default: printf("Error: Only solver 1 possible \n");
  }

  /* 3. start loop, iterate: */
  iter = 0;
  double sum = 0; // contains sum_(j=1)^(j=iter-2) of ( s_j+1 * s_j / |s_j|^2 )
  double *p_k_square = calloc(n, sizeof(double)); // contains |p_k|^2 for all k = 1,...,n
  while ( iter < maxiter ) {
    /* #1: apply correction x <- x + p */
    for (k = 0; k < n; k++) { x[k] += damping * p[k]; }
    iter++;

    for ( k = 0; k < n; k++) b_old[k] = b[k]; // save old rhs-value
    /* #2: check convergence: */
    mesh_buildRhs_nonlin(H[N], b, x, f_p, f_vol, f_H, f_g);  
    for ( k = 0; k < nfixed; k++) b[fixed[k]] = 0.0;
    printf("%2i: tol = %11g\n", iter, hpc_ddot(b, b, n));
    
    if (hpc_ddot(b, b, n) < scal_tol){
      printf("broyden converged\n");
      
      free(b), free(b_old), free(p), free(p_k), free(y), free(S);
      return iter;
    }

    /* #4: calculate change of b, y = f_new - f */
    // not needed: for (k = 0; k < n; k++) { y[k] += b[k] - b_old[k]; }

    /* #5: calculate new B matrix */
    for (k = 0; k < n; k++) { 
      /* calculate p_k = S_k * p  and p_k^2 */
      for (j = 0; j < n; j++) {p_k[j] = S[k][j]*p[j]; }
      p_k_square[k] = hpc_ddot(p_k, p_k, n);            
    }

    /* p_k is not needed anymore, sparsity pattern is contained within sed format */
    for (k = 0; k < n; k++){
      for (j = B_0[N]->i[k]; j < B_0[N]->i[k+1]; j++){
        B_0[N]->x[j] -= b[j] * p[k] / p_k_square[k];
      }
    }
  }

  /* not solved within max. no. of iterations */
  printf("\n");
  printf("WARNING: Max. No. of iterations reached\n");
  //printf("build rhs                     = %9i ns\n", (int) TIME_ELAPSED(10,11));
  //printf("fix & ddot                    = %9i ns\n", (int) TIME_ELAPSED(11,12));
  //printf("buildS_nonlin                 = %9i ns\n", (int) TIME_ELAPSED(12,13));
  //printf("buildS_byInterpolation        = %9i ns\n", (int) TIME_ELAPSED(13,14));
  //printf("solve lse (1 cycle)           = %9i ns\n", (int) TIME_ELAPSED(14,15));
  printf("==================================================\n");
  return(-iter);
}