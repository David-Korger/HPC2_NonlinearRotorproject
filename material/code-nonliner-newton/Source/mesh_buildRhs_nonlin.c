#include "hpc.h"

void rhs_stima_nonlin(double p1[2], double p2[2], double p3[2],  double u[3], 
                  index  typ, double (*fp)(double, index), double m[3])  
{
  int i, j;
  double d[3][2], fac, area2, absDu, sDu[2], Dval, tmp;
  
  for (i = 0 ; i < 2 ; i++ ){
     d[0][i] = p3[i]-p2[i];
     d[1][i] = p1[i]-p3[i];
     d[2][i] = p2[i]-p1[i]; 
     sDu[i] = d[0][i] * u[0] + d[1][i]*u[1] + d[2][i]*u[2];
  }
  area2 = d[1][0] * d[2][1] - d[2][0] * d[1][1];
  absDu = sqrt( sDu[0] * sDu[0] + sDu[1] * sDu[1] ) / area2;
  fac = fp(absDu,typ) / ( 2.0 * area2 );  
  for ( i = 0 ; i < 3 ; i++ ) m[i] = 0.0;
  for ( i = 0 ; i < 3 ; i++ ){
    for ( j = 0 ; j < i ; j++ ){ 
      tmp = fac * (d[i][0]*d[j][0] + d[i][1]*d[j][1] ); 
      m[j] += tmp * u[i];
      m[i] += tmp * u[j]; 
    }
    m[i] += u[i] * fac * (d[i][0]*d[i][0] + d[i][1]*d[i][1] ); 
  }
}

/* create rhs - vector for uniform refined macro element
 * ordering w.r.t. [ p1, p2, p3, m1=(p1+p2)/2, m2=(p2+p3)/2, m3=(p3+p1)/2] */
void mesh_Fnon_macro(double p1[2], double p2[2], double p3[2], double u[6],
                    index  typ, double (*fc)(double, index), double b[6]) 
{
  int j;
  double m[3][2], c[3], fac, s[3], uloc[3];
  
  for (j = 0 ; j < 2 ; j++ )
  {
     m[0][j] = ( p1[j] + p2[j] ) / 2.0 ;
     m[1][j] = ( p2[j] + p3[j] ) / 2.0 ;
     m[2][j] = ( p3[j] + p1[j] ) / 2.0 ; 
  }
  uloc[0] = u[0]; uloc[1] = u[3]; uloc[2] = u[5];
  rhs_stima_nonlin(p1, m[0], m[2], uloc, typ, fc, s);
  b[0]  = s[0]; b[3]  = s[1]; b[5]  = s[2];

  uloc[0] = u[1]; uloc[1] = u[4]; uloc[2] = u[3];
  rhs_stima_nonlin(p2, m[1], m[0], uloc, typ, fc, s);
  b[1]  = s[0]; b[4]  = s[1]; b[3] += s[2];

  uloc[0] = u[2]; uloc[1] = u[5]; uloc[2] = u[4];
  rhs_stima_nonlin(p3, m[2], m[1], uloc, typ, fc, s);
  b[2]  = s[0]; b[5] += s[1]; b[4] += s[2];

  uloc[0] = u[4]; uloc[1] = u[5]; uloc[2] = u[3];
  rhs_stima_nonlin(m[1], m[2], m[0], uloc, typ, fc, s);
  b[4] += s[0]; b[5] += s[1]; b[3] += s[2];
}

/* compute right hand side using midpoint rule */
void mesh_buildRhs_nonlin(const mesh *M, double *b, double *u, 
                   double (*fp)(double , index), 
                   double (*fV)(double *, index), 
                   void (*fH)(double *, index, double *), 
                   double (*fN)(double *, index))
{  
  index j, k, nC, nT, nE, nB, *Elem, *Bdry, ind[6] ;
  double *Coord, bx[6], uloc[6];
  
  nT = M->nelem; nC = M->ncoord; nB = M->nbdry; nE = M->nedges;
  Coord = M->coord; Elem = M->elem; Bdry = M->bdry; 
  
  for ( k = 0 ; k < nC + nE; k++) b[k] = 0.0;
  for ( k = 0 ; k < nT; k++)
  {
    for (j = 0 ; j < 3 ; j++) ind[j]  =      Elem[7*k+j];                                   
    for (j = 3 ; j < 6 ; j++) ind[j]  = nC + Elem[7*k+j]; 
    for (j = 0 ; j < 6 ; j++) uloc[j] = u[ ind[j] ]; 
    
    mesh_vol(Coord+2*ind[0],Coord+2*ind[1],Coord+2*ind[2],Elem[7*k+6],fV,fH,bx);
    for (j = 0 ; j < 6 ; j++) b[ind[j]] -= bx[j] ;
    
    mesh_Fnon_macro(Coord+2*ind[0],Coord+2*ind[1],Coord+2*ind[2],
                    uloc,Elem[7*k+6],fp,bx);
    for (j = 0 ; j < 6 ; j++) b[ind[j]] += bx[j] ;  
  }                    

  for ( k = 0 ; k < nB; k++)
  {
    if (Bdry[4*k+3]) 
    {
      for (j = 0 ; j < 2 ; j++) ind[j] =      Bdry[4*k+j];                                   
      ind[2] = nC + Bdry[4*k+2];  
      mesh_Neumann(Coord+2*ind[0],Coord+2*ind[1],Bdry[4*k+3],fN,bx);
      for (j = 0 ; j < 3 ; j++) b[ind[j]] -= bx[j];  
    }
  }                    
}
