#include "doe_Eval.h"
#include "doe_utility.h"
#include "doe_Matrix.h"
#include "doe_criteria.h"
#include <string.h>
#include <math.h>
#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>

using namespace Rcpp;
using namespace std;
#define DEF_PMM 20


// [[Rcpp::export]]
double CritEval(NumericMatrix X0, int q, int crit=0)
{
  switch(crit)
  {
  case 0:
    return MD2(X0,q);
  case 1:
    return CD2(X0,q);
  case 2:
    return WD2(X0,q);
  case 3:
    return maximin(X0,q);
  default:
    return MD2(X0,q);
  }
}

double MD2(NumericMatrix X0, int q)
{
  int i,j,k,nsamp=X0.nrow(),nv=X0.ncol();
  double result=1,part1=0,part2=0,mul=1;
  NumericMatrix X(clone(X0));
  for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) X(i,j) = (X(i,j)-0.5)/q;
  for(i=0;i<nv;i++) result *= 19/12.0;
  for(i=0;i<nsamp;i++)
  {
    mul=1;
    for(j=0;j<nv;j++) mul *= (5.0/3-0.25*ABS(X(i,j)-0.5)-0.25*ABS(X(i,j)-0.5)*ABS(X(i,j)-0.5));
    part1 += mul;
  }
  part1 *= (-2.0)/nsamp;
  for(i=0;i<nsamp;i++)
  {
    for(k=0;k<nsamp;k++)
    {
      mul=1;
      for(j=0;j<nv;j++) mul *=(15.0/8 - 0.25*ABS(X(i,j)-0.5) - 0.25*ABS(X(k,j)-0.5) -
          0.75*ABS(X(i,j)-X(k,j)) + 0.5*(X(i,j)-X(k,j))*(X(i,j)-X(k,j)));
      part2 += mul;
    }
  }
  part2 /= (nsamp*nsamp);
  result = result + part1 + part2;
  return result;
}

double CD2(NumericMatrix X0, int q)
{
  int i,j,k,nsamp=X0.nrow(),nv=X0.ncol();
  double result=1,part1=0,part2=0,mul=1;
  NumericMatrix X(clone(X0));
  for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) X(i,j) = (X(i,j)-0.5)/q;
  for(i=0;i<nv;i++) result *= 13/12.0;
  for(i=0;i<nsamp;i++)
  {
    mul = 1.0;
    for(k=0; k<nv; k++) mul *= (1.0+0.5*ABS(X(i,k)-0.5)-0.5*(X(i,k)-0.5)*(X(i,k)-0.5));
    part1 += mul;
  }
  part1 *= (-2.0/nsamp);
  for(i=0;i<nsamp;i++) for(j=0;j<nsamp;j++)
  {
    mul = 1.0;
    for(k=0; k<nv; k++)mul *= (1.0+0.5*ABS(X(i,k)-0.5)+0.5*ABS(X(j,k)-0.5)-0.5*ABS(X(i,k)-X(j,k)));
    part2 += mul;
  }
  part2 *= (1.0/(nsamp*nsamp));
  result = result + part1 + part2;
  
  return result;
}

double WD2(NumericMatrix X0, int q)
{
  int i,j,k,nsamp=X0.nrow(),nv=X0.ncol();
  double result=1,part1=0,mul=1;
  NumericMatrix X(clone(X0));
  for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) X(i,j) = (X(i,j)-0.5)/q;
  for(i=0;i<nv;i++) result *= 4/3.0;
  result *= (-1.0);
  for(i=0;i<nsamp;i++) for(j=0;j<nsamp;j++)
  {
    mul = 1.0;
    for(k=0; k<nv; k++) mul *= (1.5 - ABS(X(i,k)-X(j,k)) * (1-ABS(X(i,k) - X(j,k))));
    part1 += mul;
  }
  part1 *= (1.0/(nsamp*nsamp));
  result += part1;

  return result;
}

double maximin(NumericMatrix X0, int q)
{
  double d1,dt,maxmm,minmm,mmres1,mmres,**x,**D;
  int i,j,k,pmm,nsamp=X0.nrow(),nv=X0.ncol();
  pmm = DEF_PMM;
  minmm=pow(MINIDOUBLE,1.0/pmm);
  maxmm=1.0/minmm;
  x=NewDMatrix(nsamp,nv);
  D=NewDMatrix(nsamp,nsamp);

  pmm=(pmm+1)/2;

  for(i=0;i<nsamp; i++)
  {
      for(j=0;j<nv;j++)  x[i][j] = (X0(i,j)-1.0)/(q-1) ;
  }

  for(i=0;i<nsamp;i++)
  {
    for(j=i+1;j<nsamp;j++)
    {
      dt=0;
      for(k=0;k<nv;k++)
      {
        d1=x[i][k]-x[j][k];
        dt+=d1*d1;
      }
      D[j][i]=D[i][j]=dt;
      if(dt<minmm) D[j][i]=MAXDOUBLE;
      else
      {
        D[j][i]=mult(dt,pmm);//for(k=0;k<pmm-1;k++) D[j][i]*=dt;
        D[j][i]=1/D[j][i];
      }
    }
  }
  mmres1=0;
  for(i=0;i<nsamp;i++) for(j=0;j<i;j++) mmres1+=D[i][j];
  mmres1=MIN(MAXDOUBLE,mmres1);
  mmres=pow(mmres1,1.0/pmm/2.0);
  return(mmres);
}


