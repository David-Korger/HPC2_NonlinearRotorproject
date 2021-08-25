/* via x[0] < x[1] < ... < x[n] there are n 
 * intervalls defined, hpc_findspan determines k,
 * s.t. x[k] <= t < x[k+1] 
 */

int hpc_findspan(int n, double t, const double x[])
{
  int low, mid, high;
  
  if ( t >= x[0] && t < x[n] )
  { 
    low = 0; high = n;
    mid = (low+high)/2;
    while ( t < x[mid] || t >= x[mid+1] )
    {
      if (t < x[mid]) high = mid;
      else            low = mid;
      mid = (low+high)/2;
    }
    return(mid);
  } 
  else if( t >= ( 1 - 1e-14 ) * x[n] )
  {
    return(n+1);
  } 
  return -1;
}
