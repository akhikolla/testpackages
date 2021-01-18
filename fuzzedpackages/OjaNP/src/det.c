#include <stdio.h>
#include <math.h>
#include <R_ext/Lapack.h>
#define  TINY 0.0 /*1.0e-20;*/
#define NCOL 10
double det(double a[], int n)
{
 //  int i,j;  2009.05.08
   int i;
   double d=1.0;
   int ipiv[NCOL],info,neg=0; 
 //  double b[NCOL][NCOL];  2009.05.08
   
   F77_CALL(dgetrf)(&n,&n,a,&n,ipiv,&info);
  
   for(i=0;i<n;i++)
   {
     d *= a[i*n+i];
     if(ipiv[i]!=(i+1))
       neg=!neg;
   }
     
   return(neg?-d:d);
}

