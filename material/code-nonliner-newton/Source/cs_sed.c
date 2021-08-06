#include "hpc.h"
/* C = cs matrix of a sym extrated diagonal matrix T */
cs *cs_sed (const sed *T)
{
  index i, j , k, p, n, *Ti;   
  double *Tx ;
  cs *C ;
  if (!T) return (NULL) ;                             /* check inputs */
  n = T->n; Tx = T->x ; Ti = T->i ; 
  C = cs_spalloc (0, 0, 1, 1, 1) ;                    /* allocate result */
  for (k = 0; k < n; k++)
  {
    if (!cs_entry (C, (csi) k, (csi) k, Tx[k])) return (cs_spfree (C)) ;
    for (p = Ti[k]; p < Ti[k+1]; p++)  
    {
      if (!cs_entry (C, (csi) Ti[p], (csi) k, Tx[p])) return (cs_spfree (C)) ;
      if (!cs_entry (C, (csi) k, (csi) Ti[p], Tx[p])) return (cs_spfree (C)) ;
    }
  }
  return (C) ;
}