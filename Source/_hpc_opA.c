#include "hpc.h"

void hpc_opA(double *absDu, index material, double (*p)(double *, index), 
             double (*Dp)(double *, index), double *val, double *Dval)
{
  double absDu;
  absDu = sqrt(Du[0]*Du[0]+Du[1]*Du[1]);
  val = p(absDu,material);
  if (!Dval)
  {
    Dval[0] = val; Dval[1] = 0.0; Dval[2] = val; 
    if ( absDu > 1e-12 ) 
    {
      tmp = Dp(absDu,material);
      Dpx[0] += tmp * Du[0] * Du[0];
      Dpx[1] += tmp * Du[0] * Du[1];
      Dpx[2] += tmp * Du[1] * Du[1];
    }
  }
}

