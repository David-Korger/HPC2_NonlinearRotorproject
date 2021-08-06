#include "hpc.h"

#include<time.h>
#include <sys/time.h>

struct timeval tv[50];
#define TIME_SAVE(j)   (gettimeofday(&tv[j], (struct timezone*)0))
#define TIME_ELAPSED(j,k)	(1.E+6*(tv[k].tv_sec-tv[j].tv_sec)+(tv[k].tv_usec-tv[j].tv_usec))

double f_f( double x[2], index typ )
{
  double absJ = 1.5e6;
  switch(typ)
  {
    case  7: return(-absJ) ; break;
    case  8: return(-absJ) ; break;
    case  9: return( absJ) ; break;
    case 10: return( absJ) ; break;
    default: return(0.0);
  }
}

double f_g( double x[2], index typ )
{
  return ( x[0] * x[1] );
}

void f_H( double x[2], index typ, double H[2] )
{
  double C = 883310.0;
  switch(typ)
  {
    case  2: H[0] = 0.0; H[1] =   C ; break;
    case  3: H[0] = 0.0; H[1] =  -C ; break;
    case 11: H[0] =  -C; H[1] = 0.0 ; break;
    case 12: H[0] =   C; H[1] = 0.0 ; break;
    default: H[0] = 0.0; H[1] = 0.0 ;
  }
  return;
}

double f_uD( double x[2] )
{
  return ( 0.0 );
}

double f_p_steel(double t) {
  int k, n = 8;
  double tms0;
  double a0[8] = { 547.945205479452,  652.173913043478,  761.904761904762,
                   869.565217391304, 1093.750000000000, 1408.450704225352, 
                  1960.784313725490, 2555.910543130990}, 
         a1[8] = {              0.0,   724.03198402383,   975.50842352355,
                   1358.29749960185,  1976.48832733799,  3800.93954995570,
                  14111.30505587862, 25560.22408963578}, 
         a2[8] = { 4850.9707957770,    835.9843597392,   -795.0068205266, 
                   3695.4742825839,   -7216.6406108357, -60450.3660030157,
                 -79200.4153732667, 244474.7867197961}, 
         a3[8] = {-10335.520945258,      672.988283383,   18059.681339454,
                   -6758.037193997,    65393.037307002,  650397.962164913,
                 4623931.458472351,                1.86264514923096e-08},
         s0[9] = {0.73, 0.92, 1.05, 1.15, 1.28, 1.42, 1.53, 1.565, 1.6};
  double k1 = 61613.0874074867, k2 = 38104.774433735, 
         nu_inf = 57142.8571428569 ;
  /* find span */
  if ( t <= s0[0]) return(  a0[0] );
  for ( k = 0; k < 8; k++ ){
    if (t < s0[k+1]) break;
  }
  /* eval spline */
  if ( k < 8 ) 
  {
    tms0 = t - s0[k]; 
    return( a0[k] + ( a1[k] + ( a2[k] + a3[k] * tms0 ) * tms0 ) * tms0 );
  }
  else
  {
    return (nu_inf - (k1 + k2/t)/t) ;
  }
}
double f_Dp_steel(double t) {
  int k;
  double tms0;
  double a0[8] = { 547.945205479452,  652.173913043478,  761.904761904762,
                   869.565217391304, 1093.750000000000, 1408.450704225352, 
                  1960.784313725490, 2555.910543130990}, 
         a1[8] = {              0.0,   724.03198402383,   975.50842352355,
                   1358.29749960185,  1976.48832733799,  3800.93954995570,
                  14111.30505587862, 25560.22408963578}, 
         a2[8] = { 4850.9707957770,    835.9843597392,   -795.0068205266, 
                   3695.4742825839,   -7216.6406108357, -60450.3660030157,
                 -79200.4153732667, 244474.7867197961}, 
         a3[8] = {-10335.520945258,      672.988283383,   18059.681339454,
                   -6758.037193997,    65393.037307002,  650397.962164913,
                 4623931.458472351,                1.86264514923096e-08},
         s0[9] = {0.73, 0.92, 1.05, 1.15, 1.28, 1.42, 1.53, 1.565, 1.6};
  double k1 = 61613.0874074867, 
         k2 = 38104.774433735, 
         nu_inf = 57142.8571428569 ;
  /* find span */
  if ( t <= s0[0]) return(  0.0 );
  for ( k = 0; k < 8; k++ ){
    if (t < s0[k+1]) break;
  }
  /* eval spline */
  if ( k < 8 ) 
  {
    tms0 = t - s0[k]; 
    return(  a1[k] + ( 2.0*a2[k] + 3.0*a3[k] * tms0 ) * tms0 ) ;
  }
  else
  {
    return ( (k1 + 2.0*k2/t)/(t*t) ) ;
  }
}

double f_p_aloy(double t) {
  int k, n = 8;
  double tms0;
  double a0[11] = { 105.263157894737,92.3076923076923,91.9540229885057,
                   91.3461538461538, 88.9830508474576, 96.7741935483871,
          110.062893081761, 123.076923076923, 136.363636363636, 
         149.253731343284, 176.470588235294}, 
         a1[11] = {           0,-4.6601534846584,-1.21740051225639,
          -6.2397514756575,50.2520107587619,315.995518828678,
          441.679820483837,581.190298171429,654.420206659012,
          939.518387695138,3972.95160157871}, 
         a2[11] = { -173.547855189937,25.9770183550912,-12.0738235995412,
            -631.503065772,-449.040088328816,5379.52259726901,
         -2506.77080409279,8810.309500256,-15742.2277033039,
         -88472.9921583701,149893.535234457}, 
         a3[11] = {244.177494155772,-55.0078175195464,-10.5795956119619,
          3967.90171174455,29595.2146916087,-71160.5288313992,
          119000.698434009,-232652.059602214,762322.740973568,
          5476960.75018197,-32971.1066419627},
         s0[12] = {0.19,0.65,0.87,1.04,1.18,1.24,1.272,1.3,1.32,1.34,1.36,1.45};
  double k1 = 4426.40163709085, k2 = 42752.8287373295, 
         nu_inf = 25111.1111111111 ;
  /* find span */
  if ( t <= s0[0]) return(  a0[0] );
  for ( k = 0; k < 11; k++ ){
    if (t < s0[k+1]) break;
  }
  /* eval spline */
  if ( k < 11 ) 
  {
    tms0 = t - s0[k]; 
    return( a0[k] + ( a1[k] + ( a2[k] + a3[k] * tms0 ) * tms0 ) * tms0 );
  }
  else
  {
    return (nu_inf - (k1 + k2/t)/t) ;
  }
}
double f_Dp_aloy(double t) {
  int k;
  double tms0;
  double a0[11] = { 105.263157894737,92.3076923076923,91.9540229885057,
                   91.3461538461538, 88.9830508474576, 96.7741935483871,
          110.062893081761, 123.076923076923, 136.363636363636, 149.253731343284,
          176.470588235294}, 
         a1[11] = {           0,-4.6601534846584,-1.21740051225639,
          -6.2397514756575,50.2520107587619,315.995518828678,
          441.679820483837,581.190298171429,654.420206659012,
          939.518387695138,3972.95160157871}, 
         a2[11] = { -173.547855189937,25.9770183550912,-12.0738235995412,
            -631.503065772,-449.040088328816,5379.52259726901,
         -2506.77080409279,8810.309500256,-15742.2277033039,
         -88472.9921583701,149893.535234457}, 
         a3[11] = {244.177494155772,-55.0078175195464,-10.5795956119619,
          3967.90171174455,29595.2146916087,-71160.5288313992,
          119000.698434009,-232652.059602214,762322.740973568,
          5476960.75018197,-32971.1066419627},
         s0[12] = {0.19,0.65,0.87,1.04,1.18,1.24,1.272,1.3,1.32,1.34,1.36,1.45};
  double k1 = 4426.40163709085, k2 = 42752.8287373295, 
         nu_inf = 25111.1111111111 ;
  /* find span */
  if ( t <= s0[0]) return(  0.0 );
  for ( k = 0; k < 11; k++ ){
    if (t < s0[k+1]) break;
  }
  /* eval spline */
  if ( k < 11 ) 
  {
    tms0 = t - s0[k]; 
    return(  a1[k] + ( 2.0*a2[k] + 3.0*a3[k] * tms0 ) * tms0 ) ;
  }
  else
  {
    return ( (k1 + 2.0*k2/t)/(t*t) ) ;
  }
}



double f_p(double t, index typ)
{
  double val, mu0 = 4.0e-7 * PI;
  switch(typ)
  {
    case  0: val = 1.0 / ( mu0 * 1.0 ) ; break; // Air
    case  1: val = f_p_aloy(t) ; break; // Cobalt Nickel Copper alloy
    case  2: val = 1.0 / ( mu0 *  1.045 ) ; break; // Magnet (NdFeB)
    case  3: val = 1.0 / ( mu0 *  1.045 ) ; break; // Magnet (NdFeB)
    case  4: val = f_p_steel(t) ; break; // Steel
    case 11: val = 1.0 / ( mu0 *  1.045 ) ; break; // Magnet (NdFeB)
    case 12: val = 1.0 / ( mu0 *  1.045 ) ; break; // Magnet (NdFeB)
    default: val = 1.0 / ( mu0 *  1.0 ) ;
  }
  return (val);
}

double f_Dp(double t, index typ)
{
  double val;
  switch(typ)
  {
    case  1: val = f_Dp_aloy(t) ; break; // Cobalt Nickel Copper alloy
    case  4: val = f_Dp_steel(t) ; break; // Steel
    default: val = 0.0 ;
  }
    return(val);
}

double kappa(double mid[2], index typ)
{
  return (1.0 );
}
 

int main (int argc, char **argv)
{
    index n, k, ncoord, nelem, nbdry, nfixed, nedges, total, *bdry,
          cnt = 0, N = 0, MAX_ITER = 100;
    char *Pdir = "../Problem/", fname[64];
    double *b, *x, *w, *Coord, x1[2], x2[2], m[2];
    mesh **H, *T ;
    sed **A;
    
    
    TIME_SAVE(0);
    printf("\n========================================\n");
    if (argc < 2 ){ printf("Problem not specified\n"); return(1); } 
    sprintf(fname,"%s%s",Pdir,argv[1]); /* get problem as parameter */
    if (argc>2){                        /* get no. of refinements  */ 
      if ( (atoi(argv[2]) >0) & (atoi(argv[2]) < 13)) N = atoi(argv[2]);
    } 
    printf("Load data form %s, no. refinements = %g\n", fname, (double) N);   
        
    /* Allocate memory for hierachy */
    H = malloc ( (N+1) * sizeof(mesh));
    A = malloc ( (N+1) * sizeof(sed*));
    
    /* Load geometry */
    H[0] = mesh_load (fname);             
    mesh_getEdge2no(H[0]->nelem, H[0]->elem, &H[0]->nedges, &H[0]->edge2no);
    H[0]->fixed = mesh_getFixed(H[0]->ncoord, H[0]->bdry, H[0]->nbdry, &H[0]->nfixed);
    printf("\nInit mesh  # dofs =  %10g\n", (double) H[0]->ncoord+H[0]->nedges);
    
    /* Build stiffness matrix, refine mesh and create hierachy  */ 
    TIME_SAVE(1);
    k = 0;
    while(1)
    {  
      A[k] = sed_nz_pattern(H[k]) ;            /* get pattern of matrix -> lookup in sed_buildS.c! */
      if (!A[k]) return(1);
      if (k >= N) break;
      H[k+1] = mesh_refine(H[k]);
      mesh_getEdge2no(H[k+1]->nelem, H[k+1]->elem,
                      &H[k+1]->nedges, &H[k+1]->edge2no);
      H[k+1]->fixed = mesh_getFixed(H[k+1]->ncoord, H[k+1]->bdry, 
                                    H[k+1]->nbdry, &H[k+1]->nfixed);
      k++;
    }
    TIME_SAVE(2);

    n = A[N]->n;
    x = calloc (n, sizeof(double));       /* get workspace for sol*/

//    cnt = hpc_solve_nonlin(H, A, N, x, f_uD, f_f, f_H, f_g, f_p, f_Dp, 2, 1e-10, 30);          
    cnt = hpc_solve_nonlin(H, A, N, x, f_uD, f_f, f_H, f_g, f_p, f_Dp, 1, 1e-10, 30);          
    TIME_SAVE(3);
    hpc_print_double(x, A[N]->n, "x", 1);
    sprintf(fname,"%s%s",Pdir,"motor.sol"); 
    hpc_save_double(fname,x,A[N]->n);
    
    printf("Final mesh # dofs =  %10g\n",(double)  H[N]->ncoord+H[N]->nedges);
    printf("# refinements     =  %10g\n",(double)  N);
                                   
    printf("\n");
    printf("Time load                     = %9i ns\n", (int) TIME_ELAPSED(0,1));
    printf("Time create hierarchy         = %9i ns\n", (int) TIME_ELAPSED(1,2));
    printf("Time solve nonlinear system   = %9i ns\n", (int) TIME_ELAPSED(2,3));
    printf("No. iterations                = %9g\n", (double) cnt); 
    printf("========================================\n\n");
     
    for (k=0; k<=N; k++) {mesh_free(H[k]); sed_free(A[k]);}
    free(H); free(A);

    return (0) ;
}