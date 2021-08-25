#include "hpc.h"


/* print (double) vector */
void hpc_print_double(double *v, index n, char *txt, index brief)
{
  index k;
  for (k = 0 ; k < n ; k++)
  {
     printf("%s[%g] = %g\n", txt, (double) k, v[k]);
     if (brief && k > 10) { printf ("  ... \n%s[%g] = %g\n\n", txt, (double) n, v[n]) ; break ; }
  }
}  

/* print (index) vector */
void hpc_print_index(index *v, index n, char *txt, index brief)
{
  index k;
  for (k = 0 ; k < n ; k++)
  {
     printf("%s[%g] = %g\n", txt, (double) k, (double) v[k]);
     if (brief && k > 10) { printf ("  ... \n%s[%g] = %g\n\n", txt, (double) n, (double) v[n]) ; break ; }
  }
}  

/* calculate inner vector product or scalar product */
double hpc_ddot(double *u, double *v, index n)
{
  index k;
  double tmp = 0.0;
  for (k = 0 ; k < n ; k++) tmp += u[k] * v[k];
  return (tmp);
}  

/* save double vector v of length n to file fname */
index hpc_save_double(char *fname, double *v, index n)
{
  FILE *file;
  index k;
  
  file = fopen(fname,"w");
  
  if ( !file ) return (0);
  for ( k = 0; k < n; k++) {
    fprintf(file,"%g \n", v[k]);
  }
  fclose(file);
  return (1);
}
