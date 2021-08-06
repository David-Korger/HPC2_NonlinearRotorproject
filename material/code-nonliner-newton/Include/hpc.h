#ifndef _HPC_H
#define _HPC_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

#include <errno.h>
#include <string.h>

#include <stdbool.h>
#include "cs.h"


#define index ptrdiff_t

#define PI 3.141592653589793

/* --- primary HPC routines and data structures ------------------------- */


typedef struct sed_sparse  /* matrix in sparse matrix in compressed col. */
{                          /* with extracted diagonal storage form      */
    index nzmax ;     /* maximum number of entries */
    index   n ;       /* number of rows/columns          */
    index  *i ;       /* col pointers and row indices    */
    double *x ;       /* numerical values, size i[n] */
} sed ;

typedef struct mesh_data  /* mesh */
{
    index ncoord ;    /* number of coordinates  */
    index nelem ;     /* number of elements   */
    index nedges ;    /* number of edges  */
    index nbdry ;     /* number of boundary elements  */
    index nfixed;     /* number of fixed nodes ????    */
    double *coord ;   /* coordinates (x1,y1,x2,y2, ... ,x_ncoord,y_ncoord) */
    index *elem ;     /* elements ([e1,e2,e3,m1,m2,m3,t1], ... ) */
    index *edge2no ;  /*  */
    index *bdry ;     /* bdry ([e1,e2,m1,t1], [e3,e4,m2,t2], ...) */
    index *fixed ;    /* bdry ([e1,e2,m1,t1], [e3,e4,m2,t2], ...) */
} mesh ;

typedef struct sky_pack  /* sym. matrix in sky storage form */
{
    index   n ;       /* number of rows/columns          */
    index  *p ;       /* col pointers (size n)         */
    double *d ;       /* diagonal entries (size n)       */
    double *x ;       /* off-diagonal entries, size p[n-1] */
} sky ;

/* utilities */
void *hpc_realloc (void *p, index n, size_t size, index *ok);
double hpc_cumsum (index *p, index *c, index n);
 
sed *sed_alloc (index n, index nzmax, index values);
index sed_realloc (sed *A, index nzmax);
sed *sed_free (sed *A);
sed *sed_done (sed *C, void *w, void *x, index ok);
// sed *sed_compress (const cs *A);
index sed_print (const sed *A, index brief);
index sed_gaxpy (const sed *A, const double *x, double *y);
index sed_dupl (sed *A);
index sed_gs_constr (const sed *A, const double *b, double *x, double *w, 
                     index *fixed, index nFixed, index forward);

mesh *mesh_alloc (index ncoord, index nelem, index nbdry);
mesh *mesh_free (mesh *M);
mesh *mesh_load (char *fname);
index *mesh_getFixed(const index nCoord, const index *bdry, 
                     const index nBdry, index *nFixed);
index mesh_print (const mesh *M, index brief);
mesh *mesh_refine(mesh *In);
index mesh_getEdge2no(const index nElem, const index *Elem, 
                      index *nEdges, index **edge2no);




void stima_laplace3(double p1[2], double p2[2], double p3[2],
                    index  typ, double dx[6], double ax[9]);

sed *sed_nz_pattern(mesh *M) ; 
index sed_buildS(mesh *M, sed *T);
index sed_buildS_nonlin(mesh *M, sed *T, double *u, double (*fp)(double, index),
                    double (*fDp)(double, index)) ;

void mesh_Neumann(double p1[2], double p2[2], 
                  index typ, double (*fc)(double *, index), double m[3]);  
void mesh_vol_elem(double p1[2], double p2[2], double p3[2], index typ, 
        double (*fc)(double *, index), void (*fdiv)(double *, index, double *), 
        double m[3]);  
void mesh_vol(double p1[2], double p2[2], double p3[2], index  typ, 
     double (*fc)(double *, index), void (*fdiv)(double *, index, double *), 
     double b[6]); 

void mesh_buildRhs(const mesh *M, double *b, double (*fV)(double *, index), 
                   void (*fDiv)(double *, index, double *), 
                   double (*g)(double *, index));
void mesh_buildRhs_nonlin(const mesh *M, double *b, double *u, 
                   double (*fp)(double , index), 
                   double (*fV)(double *, index), 
                   void (*fDiv)(double *, index, double *), 
                   double (*fN)(double *, index));


void hpc_fmg(sed **A, double *b, double *x, index nCycle,
             mesh **H, index nLevel, index pre, index post, index gamma);          
index hpc_mg(sed **A, double *b, double *x, double tol, index maxit,
             mesh **H, index nLevel, index pre, index post, index gamma);          
index hpc_mg_cycle(sed **A, mesh **H, index nLevel, 
                   double **b, double **x, double **r,
                   index pre, index post, index gamma);

index hpc_solve_nonlin(mesh **H, sed **A, index nlevel, double *x,         
         double (*f_uD)(double *),  double (*f_vol)(double *, index),
         void (*f_H)(double *, index, double *), double (*f_g)(double *, index),    
         double (*f_p)(double, index),     double (*f_Dp)(double, index), 
         index solver, double tol, index maxiter); 


void hpc_rest(double *x, index *edgeno, index nEdges, double *y, index ny);
void hpc_prol(double *x, index nx, index *edgeno, index nEdges, double *y);
void hpc_prol_quad(double *x, double *y, index *elem, index nC, index nT, index nE);

index hpc_buildS_byInterpolation(sed *Ain, sed *Bout, index nE, index *edge2no);

void hpc_print_double(double *v, index n, char *txt, index brief);
void hpc_print_index(index *v, index n, char *txt, index brief);
double hpc_ddot(double *u, double *v, index n);
index hpc_save_double(char *fname, double *v, index n);

/* skyline format */
sky *sky_alloc (index n, index nzmax);
sky *sky_free (sky *A);
sky *sky_done (sky *G, void *w, void *x, index ok);
index sky_print (const sky *A, index brief);
sky *sky_sed (const sed *T);
index sky_cholesky(sky *A); 
index sky_cholsol(sky *A, double *x); 


cs *cs_sed (const sed *T);
index sed_cholsol_res(sed *A, double *b, index *fixed, index nfixed);

double kappa( double x[2], index typ );
double F_vol( double x[2], index typ );

#define HPC_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define HPC_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define HPC_CSC(A) (A && (A->nz == -1))
#define HPC_CSR(A) (A && (A->nz == -2))
#define HPC_TRIPLET(A) (A && (A->nz >= 0))
#endif
