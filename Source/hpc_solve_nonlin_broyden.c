#include "hpc.h"

#include<time.h>
#include <sys/time.h>

struct timeval tv_solve_broyden[50];
#define TIME_SAVE(j)   (gettimeofday(&tv_solve_broyden[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv_solve_broyden[k].tv_sec - tv_solve_broyden[j].tv_sec)+(tv_solve_broyden[k].tv_usec - tv_solve_broyden[j].tv_usec))

index hpc_solve_nonlin_broyden(mesh **H, sed **B_0, index N, double *x,         
         double (*f_uD)(double *),  double (*f_vol)(double *, index),
         void (*f_H)(double *, index, double *), double (*f_g)(double *, index),   
         double (*f_p)(double, index),     double (*f_Dp)(double, index), 
         index solver, double tol, index maxiter, double damping) 
{
  index k, j, n, ncoord, nbdry, nfixed, *fixed, nedges, *bdry, iter = 0;
  double x1[2], x2[2], m[2], scal_tol, *w, *s, *s_old, *b, *Coord;
  
  n = B_0[N]->n; nbdry = H[N]->nbdry ; bdry = H[N]->bdry; ncoord = H[N]->ncoord ; 

//   nelem = H[N]->nelem ; nedges = H[N]->nedges; 
//   H[N]->fixed = mesh_getFixed(ncoord, bdry, nbdry, &H[N]->nfixed);
  nfixed = H[N]->nfixed ; fixed = H[N]->fixed ; 
  
  Coord = H[N]->coord;
  
  w = calloc (n, sizeof(double));      
  b = calloc (n, sizeof(double));       
  s = calloc (n, sizeof(double));       
  s_old = calloc (n, sizeof(double));       

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

  /* 2. Solve B_0 s = -F(x_0),  for s */
  /* b <- F(x_0) */
  mesh_buildRhs_nonlin(H[N], b, x, f_p, f_vol, f_H, f_g);  
  /* empty solution vector s */
  for (k = 0; k < n; k++) s[k] = 0.0;
  /* s <- solve B_0 s = -F(x_0) */
  switch(solver)
  {
   case 1: hpc_fmg(B_0, b, s, 3, H, N, 2, 2, 1); break;
   case 2: hpc_mg(B_0, b, s, 1e-10, 100, H, N, 2, 2, 1); break;
   default: printf("Error: Only solver 1 possible \n");
  }

  /* 3. start loop, iterate: */
  iter = 0;
  double sum = 0; // contains sum_(j=1)^(j=iter-2) of ( s_j+1 * s_j / |s_j|^2 )
  double s_square = hpc_ddot(s, s, n); // contains |s_iter|^2
  double temp = 0;
  while ( iter < maxiter ) {
    /* #1: apply correction x <- x + s_iter */
    for (k = 0; k < n; k++) { x[k] += damping * s[k]; }
    iter++;

    /* #2: check convergence: */
    mesh_buildRhs_nonlin(H[N], b, x, f_p, f_vol, f_H, f_g);  
    for ( k = 0; k < nfixed; k++) b[fixed[k]] = 0.0;
    printf("%2i: tol = %11g\n", iter, hpc_ddot(b, b, n));
    
    if (hpc_ddot(b, b, n) < scal_tol){
      printf("broyden converged\n");
      
      free(b), free(w);
      return iter;
    }
    
    /* #3: rhs (=b) changed, solve for new correction w: */
    for (k = 0; k < n; k++) w[k] = 0.0;
    switch(solver)
    {
     case 1: hpc_fmg(B_0, b, w, 3, H, N, 2, 2, 1); break;
     case 2: hpc_mg(B_0, b, w, 1e-10, 100, H, N, 2, 2, 1); break;
     default: printf("Error: Only solver 1 possible \n");
    }

    /* #4: correct correction w by sum */
    s_square = hpc_ddot(s, s, n);
    if(iter > 2) {
      sum += hpc_ddot(s_old, s, n) / s_square; // for iter = 0, s_old is empty (due to calloc)
      for (k = 0; k < n; k++){ w[k] += sum; }
    }

    /* #5: get new offset s_iter for x */
    for (k = 0; k < n; k++){ s_old[k] = s[k]; }  // save old value to calc. sum next time    

    temp = hpc_ddot(s_old, w, n);
    for (k = 0; k < n; k++){ s[k] = w[k] / (1 + temp/s_square); }

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