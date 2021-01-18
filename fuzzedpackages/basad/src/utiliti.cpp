//
//
//  Created by Qingyan Xiang on 4/19/17.
//  Copyright © 2017 Qingyan Xiang. All rights reserved.
//
#include <Rcpp.h>
#include <math.h>
#include <RcppEigen.h>
#include "utiliti.h"


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define P1 10e-4
#define P2 10e-4

using namespace std;
using namespace Eigen;
using namespace Rcpp;


double mydnorm (double x, double mu, double sigmasq)
{ double pdf;
  pdf=(1/sqrt(2*PI*sigmasq))*exp(-(x-mu)*(x-mu)/(2*sigmasq));
  return(pdf);
}


/*unif(0,1) generator*/
double ran1(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1; 
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { 
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1; 
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 
  if (*idum < 0) *idum += IM1; 
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2; 
  iv[j] = *idum; 
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

/*N(0,1) generator */
double gasdev(long *idum){
  double ran1(long *idum);
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  if (*idum < 0) iset=0;
  if (iset == 0){ 
    do {
      v1=2.0*ran1(idum)-1.0; 
      v2=2.0*ran1(idum)-1.0;  
      rsq=v1*v1+v2*v2; 
    } while (rsq >= 1.0 || rsq == 0.0); 
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1; 
    return v2*fac;
  } else { 
    iset=0; 
    return gset;
  }
}

/*Gamma generator  */

double gamdev(double a, double b, long *idum)// a:alpha,shape parameter; b: beta, rate parameter. 
{double gasdev(long *idum),tempa=a;
  if(a<1) tempa=a+1;
  double ran1(long *idum);
  double d=tempa-1.0/3.0;
  double c=1.0/sqrt(9.0*d),u,z,x,v;
  int flag=1;
  do{
    u=ran1(idum);z=gasdev(idum);
    v=pow(1.0+c*z,3);
    if(z>-1.0/c){
      if((u<1.0-0.0331*z*z*z*z) | (log(u)<pow(z,2.0)/2.0+d-d*v+d*log(v) ) )  
      {flag=0;x=d*v;}
    }
    
  }while(flag);
  x=x/b;
  if(a<1) {u=ran1(idum);x*=pow(u,1/a);}
  return x;
}


double betadev( double a, double b, long *idum ){
    double gamdev( double a1, double b1, long *idum );
    
    double x, y, beta;
    x = gamdev( a, 1.0, idum );
    y = gamdev( b, 1.0, idum );
    
    beta = x /( x + y );
    
    return beta;
}



/*truncated norm, N(0,1), truncated at a, larger than a  */
double tndev(double a, long *idum)
{double gasdev(long *idum);
  double ran1(long *idum);
  int flag=1;
  double u, v, z;
  if(a<=0.45){
    do{
      //see++;
      u=gasdev(idum); if(u>a) flag=0;
    }while(flag);z=u;
  }
  else if(a>100){
      //cout<<"warning! large in truncated gamma"<<endl;
      return a;
  }
  else
  {
    do{//see++;
      u=ran1(idum); v=ran1(idum); z=a-log(1-u)/a;
      if(log(v)<-(z-a)*(z-a)/2) flag=0;
    }while(flag);
  }//cout<<see<<";a: "<<a<<endl;
  return z;
}


//inverse gaussian generator
double igasdev(double u, double l, long *idum){
    double gasdev(long *idum);
    double ran1(long *idum);
    
    u = ( u < 100000)? u:100000;
    
    double v, y, x, z, res;
    v = gasdev(idum);
    y = v * v;
    x = u + u * u * y /(2 * l) - u * sqrt( 4 * u * l * y + u * u * y * y )/ ( 2 * l );
    z = ran1(idum);
    
    if( z <= (  u/(u+x)  ))
        res = x;
    else
        res = (u*u)/x;

    return (res < 10000)? res:10000;
}


