#include "hpc.h"

/* compute local stiffness matrix for nonlinear problem */

void stima_nonlin(double p1[2], double p2[2], double p3[2],  double u[3], 
                  index  typ, double (*fp)(double, index),
                   double (*fDp)(double, index), double m[3][3])  
{
  int i, j;
  double d[3][2], fac, area2, absDu, srot[2], Dval[3], tmp;
    
  for (i = 0 ; i < 2 ; i++ ){
     d[0][i] = p3[i]-p2[i];
     d[1][i] = p1[i]-p3[i];
     d[2][i] = p2[i]-p1[i]; 
     srot[i] = d[0][i] * u[0] + d[1][i]*u[1] + d[2][i]*u[2];
  }
  area2 = d[1][0] * d[2][1] - d[2][0] * d[1][1];
  absDu = sqrt( srot[0] * srot[0] + srot[1] * srot[1] ) / area2;
  Dval[0] = fp(absDu,typ); Dval[1] = Dval[0]; Dval[2] = 0.0 ;
  if ( absDu > 1e-12 ) 
  {
    tmp = fDp(absDu,typ) / (absDu * area2 * area2) ;
    Dval[0] += tmp * srot[0] * srot[0];
    Dval[1] += tmp * srot[1] * srot[1];
    Dval[2] += tmp * srot[0] * srot[1];
  }
  fac = 1.0 / ( 2.0 * area2 );  
  for ( i = 0 ; i < 3 ; i++ ){
    for ( j = 0 ; j < i ; j++ ){ 
      m[i][j] = fac * ( Dval[0] *  d[i][0]*d[j][0] + Dval[1]*d[i][1]*d[j][1]
                      + Dval[2] * (d[i][0]*d[j][1] +         d[i][1]*d[j][0])); 
      m[j][i] = m[i][j]; 
    }
    m[i][i] = fac * (     Dval[0]*d[i][0]*d[i][0] + Dval[1]*d[i][1]*d[i][1] 
                    + 2 * Dval[2]*d[i][0]*d[i][1] ); 
  }
  //printf("d1 = %f, d2 = %f, d3 = %f\n", m[0][0], m[1][1], m[2][2]); 
}

// void stima_laplace42_(double p1[2], double p2[2], double p3[2],  
//                    index  typ,  double m[3][3])  
// {
//   int i, j;
//   double d[3][2], fac, mid[2];
//   
//   for (i = 0 ; i < 2 ; i++ ){
//      d[0][i] = p3[i]-p2[i];
//      d[1][i] = p1[i]-p3[i];
//      d[2][i] = p2[i]-p1[i]; 
//      mid[i] = ( p1[i] + p2[i] + p3[i] ) / 3.0 ; 
//   }
//   fac = 1/(2.0*(d[1][0]*d[2][1]-d[2][0]*d[1][1])) * kappa(mid,typ);
//   for ( i = 0 ; i < 3 ; i++ ){
//     for ( j = 0 ; j < i ; j++ ){ 
//       m[i][j] = fac * (d[i][0]*d[j][0] + d[i][1]*d[j][1]); 
//       m[j][i] = m[i][j]; 
//     }
//     m[i][i] = fac * (d[i][0]*d[i][0] + d[i][1]*d[i][1]); 
//   }
// }



/* create stiffness matrix for uniform refined element
 * ordering w.r.t. [ p1, p2, p3, m1=(p1+p2)/2, m2=(p2+p3)/2, m3=(p3+p1)/2] */
void stima_macro_nonlin(double p1[2], double p2[2], double p3[2],
                        double u[6], index  typ, double (*fp)(double, index),
                    double (*fDp)(double, index), double dx[6], double ax[9]) 
{
  int i,j;
  double m[3][2], c[3], fac, uloc[3], s[3][3];
  
  for (j = 0 ; j < 2 ; j++ )
  {
     m[0][j] = ( p1[j] + p2[j] ) / 2.0 ;
     m[1][j] = ( p2[j] + p3[j] ) / 2.0 ;
     m[2][j] = ( p3[j] + p1[j] ) / 2.0 ; 
  }
  uloc[0] = u[0]; uloc[1] = u[3]; uloc[2] = u[5];
  stima_nonlin(p1, m[0], m[2], uloc, typ, fp, fDp, s);  
  dx[0]  = s[0][0]; dx[3]  = s[1][1]; dx[5]  = s[2][2];
  ax[0]  = s[0][1]; ax[1]  = s[0][2]; ax[7]  = s[1][2];
  uloc[0] = u[1]; uloc[1] = u[4]; uloc[2] = u[3];
  stima_nonlin(p2, m[1], m[0], uloc, typ, fp, fDp, s);  
  dx[1]  = s[0][0]; dx[4]  = s[1][1]; dx[3] += s[2][2];
  ax[3]  = s[0][1]; ax[2]  = s[0][2]; ax[6]  = s[1][2];
  uloc[0] = u[2]; uloc[1] = u[5]; uloc[2] = u[4];
  stima_nonlin(p3, m[2], m[1], uloc, typ, fp, fDp, s);  
  dx[2]  = s[0][0]; dx[5] += s[1][1]; dx[4] += s[2][2];
  ax[5]  = s[0][1]; ax[4]  = s[0][2]; ax[8]  = s[1][2];
  uloc[0] = u[4]; uloc[1] = u[5]; uloc[2] = u[3];
  stima_nonlin(m[1], m[2], m[0], uloc, typ, fp, fDp, s);  
  dx[4] += s[0][0]; dx[5] += s[1][1]; dx[3] += s[2][2];
  ax[8] += s[0][1]; ax[6] += s[0][2]; ax[7] += s[1][2];  
}

index sed_buildS_nonlin(mesh *M, sed *T, double *u, double (*fp)(double, index),
                    double (*fDp)(double, index))                         
{
  index j, k, n, p, nC, nT, nz, *Elem, ind[6], *Ti, *w, imin, imax;
  
  static int ai[9] = {0,0,1,1,2,2,3,3,4}, aj[9] = {3,5,3,4,4,5,4,5,5};
  double dx[6], ax[9], *Coord, *Tx, uloc[6];

  nT = M->nelem; nC = M->ncoord; Coord = M->coord; Elem = M->elem; 
  n = T->n ; Ti = T->i ; 
  
  if (!(T->x)) T->x = malloc( Ti[n] * sizeof (double)) ;
  if (!(T->x)) return(0);
  Tx = T->x;
  
  for ( k = 0 ; k < Ti[n]; k++) Tx[k] = 0.0;

  for ( k = 0 ; k < nT; k++)
  {
    for (j = 0 ; j < 3 ; j++) ind[j] =      Elem[7*k+j];                                   
    for (j = 3 ; j < 6 ; j++) ind[j] = nC + Elem[7*k+j];
    for (j = 0 ; j < 6 ; j++) uloc[j] = u[ ind[j] ];
    stima_macro_nonlin(Coord+2*ind[0],Coord+2*ind[1],Coord+2*ind[2],
                       uloc, Elem[7*k+6],fp,fDp,dx,ax);
    for (j = 0 ; j < 6 ; j++) Tx[ind[j]] += dx[j];
    for (j = 0 ; j < 9 ; j++)
    {
      imin = HPC_MIN( ind[ai[j]] , ind[aj[j]] );
      imax = HPC_MAX( ind[ai[j]] , ind[aj[j]] );

      for (p = Ti[imin] ; p < Ti[imin+1] ; p++)
      {
        if (Ti[p] == imax ) 
        {
          Tx[p] += ax[j];
          break;
        }
      }
    }
  }
  return(1);
}
