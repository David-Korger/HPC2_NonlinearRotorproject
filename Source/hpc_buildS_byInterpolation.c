#include "hpc.h"

index set_value(index *Ti, double *Tx, index first, index last, index idx, 
                double value)
{
  index p;
  for (p = first ; p < last ; p++)
  {
    if (Ti[p] == idx ) 
    {
      Tx[p] += value;
      return(1);
    }
  }
  return(0);
}

index hpc_buildS_byInterpolation(sed *Ain, sed *Bout, index nE, index *edge2no)
{
  index i, j, k, p, nC, ind, imin, imax, i0, i1, idx[2], *Ai, *Bi;
  double *Ax, *Bx;
  
  nC = Bout->n ;                               /* Small dimension */
  Ai = Ain->i;  Ax = Ain->x; Bi = Bout->i; 
  
  if (!(Bout->x)) Bout->x = malloc( Bi[nC] * sizeof (double)) ;
  
  if (!Bout->x) return(0);
  Bx = Bout->x;  

  for (k = nC; k < Bi[nC]; k++) Bx[k] = 0.0;   /* Init Bx */  
  for (k = 0; k < nC; k++) Bx[k] = Ax[k];      /* Copy D_cc */
  for (k = 0; k < nE; k++){                    /* Compute D_cc += D_ff * P_f */
    i0 = edge2no[ 2 * k ];       
    i1 = edge2no[ 2 * k + 1 ];
    imin = HPC_MIN(i0,i1);       
    imax = HPC_MAX(i0,i1);
    Bx[i0] += 0.25 * Ax[nC+k];    
    Bx[i1] += 0.25 * Ax[nC+k];  
    set_value(Bi, Bx, Bi[imin], Bi[imin+1], imax, 0.25*Ax[nC+k]);  
  } 
  for (k = 0; k < nC; k++)                     /* Compute A_cc + P_f * A_fc */
  {
    for (p = Ai[k]; p < Ai[k+1]; p++) 
    {
      if (Ai[p] < nC)
      {
        set_value(Bi, Bx, Bi[k], Bi[k+1], Ai[p], Ax[p]); 
      }
      else
      {
        for (j = 0; j < 2; j++)
        {
          ind = edge2no[ 2 * ( Ai[p] - nC ) + j ];
          if (ind == k)
          {
            Bx[ind] += Ax[p];                  /* 2.0 * 0.5 * Ax[p]  */ 
          }
          else
          {
            imin = HPC_MIN(ind,k);
            imax = HPC_MAX(ind,k);
            set_value(Bi, Bx, Bi[imin], Bi[imin+1], imax, 0.5*Ax[p]); 
          }
        }
      }         
    }
  }
  for (k = 0; k < nE; k++)             /* Compute A_cc += P_f * A_ff * P_f^T */
  {
    idx[0] = edge2no[ 2 * k ]; 
    idx[1] = edge2no[ 2 * k + 1 ];
    for (p = Ai[nC+k]; p < Ai[nC+k+1]; p++) 
    {
      for (j = 0; j < 2; j++)
      {
        ind = edge2no[ 2 * ( Ai[p] - nC ) + j ];
        for (i = 0; i < 2; i++)
        {
          if (ind == idx[i])
          {
            Bx[ind] += 0.5 * Ax[p];      /* 2.0 * 0.25 * Ax[p] */;
          }
          else
          {
            imin = HPC_MIN(ind,idx[i]);
            imax = HPC_MAX(ind,idx[i]);
            set_value(Bi, Bx, Bi[imin], Bi[imin+1], imax, 0.25*Ax[p]);  
          }
        }
      }         
    }
  }
  return(1);
}
  


