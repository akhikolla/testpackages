#include "curstat.h"
#include <R_ext/Random.h>

void convexmin(int n, double cumw[], double cs[], double y[])
{
  int	i, j, m;

  y[1] = cs[1]/cumw[1];
  for (i=2;i<=n;i++)
  {
    y[i] = (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1]);
    if (y[i-1]>y[i])
    {
      j = i;
      while (y[j-1] > y[i] && j>1)
      {
        j--;
        if (j>1)
          y[i] = (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1]);
        else
          y[i] = cs[i]/cumw[i];
        for (m=j;m<i;m++)	y[m] = y[i];
      }
    }
  }
}

double KK(double x)
{
  double u,y;

  u=x*x;

  if (u<=1)
    y = (16.0 + 35*x - 35*pow(x,3) + 21*pow(x,5) - 5*pow(x,7))/32.0;
  else
  {
    if (x>1)
      y=1;
    else
      y=0;

  }

  return y;
}


double bdf(double A, double B, int m, double t[], double p[], double u, double h)
{
  int			k;
  double		t1,t2,t3,sum;


  sum=0;
  for (k=1;k<=m;k++)
  {
    t1=(u-t[k])/h;
    t2=(u+t[k]-2*A)/h;
    t3=(2*B-u-t[k])/h;
    sum+= (KK(t1)+KK(t2)-KK(t3))*p[k];
    //sum+= KK(t1)*p[k];
  }
  return fmax(0,sum);
}


int compare(const void *a, const void *b)
{
  double x = *(double*)a;
  double y = *(double*)b;

  if (x < y)
    return -1;
  if (x > y)
    return 1;
  return 0;
}

int CompareTime(const void *a, const void *b)
{
  if ((*(SampleTime *) a).t < (*(SampleTime *) b).t)
    return -1;
  if ((*(SampleTime *) a).t > (*(SampleTime *) b).t)
    return 1;
  return 0;
}



void data_bootstrap(int N, int n, int *m, double x[], double x2[], double data2[], int **freq, int delta[], int delta2[])
{
  int	i,j,k;
  SampleTime2 *obs;

  obs = new SampleTime2[N+1];

  for (k=1;k<=n;k++)
    freq[0][k]=freq[1][k]=0;

  for (k=1;k<=N;k++)
  {
    j= 1+ N*unif_rand();

    x2[k]=x[j];
    delta2[k]=delta[j];
  }

  for (i=0;i<N;i++)
  {
    obs[i].t=x2[i+1];
    obs[i].delta=delta2[i+1];
  }

  qsort(obs,N,sizeof(SampleTime2),CompareTime);

  for (i=1;i<=N;i++)
  {
    x2[i]=obs[i-1].t;
    delta2[i]=obs[i-1].delta;
  }


  x2[0]=0;

  j=0;

  for (i=1;i<=N;i++)
  {
    if (x2[i]>x2[i-1])
    {
      j++;
      data2[j]=x2[i];
      freq[delta2[i]][j]=1;
    }
    else
    {
      data2[j]=x2[i];
      freq[delta2[i]][j]++;
    }
  }

  *m=j;


  delete[] obs;
}


void data_bootstrap2(int N, int nB, int n, int *m, double x[], double x2[], double data2[], int **freq, int delta[], int delta2[])
{
  int	i,j,k;
  SampleTime2 *obs;

  obs = new SampleTime2[nB+1];

  for (k=1;k<=nB;k++)
    freq[0][k]=freq[1][k]=0;

  for (k=1;k<=nB;k++)
  {
     j= 1+ N*unif_rand();

    x2[k]=x[j];
    delta2[k]=delta[j];
  }

  for (i=0;i<nB;i++)
  {
    obs[i].t=x2[i+1];
    obs[i].delta=delta2[i+1];
  }

  qsort(obs,nB,sizeof(SampleTime2),CompareTime);

  for (i=1;i<=nB;i++)
  {
    x2[i]=obs[i-1].t;
    delta2[i]=obs[i-1].delta;
  }


  x2[0]=0;

  j=0;

  for (i=1;i<=nB;i++)
  {
    if (x2[i]>x2[i-1])
    {
      j++;
      data2[j]=x2[i];
      freq[delta2[i]][j]=1;
    }
    else
    {
      data2[j]=x2[i];
      freq[delta2[i]][j]++;
    }
  }

  *m=j;


  delete[] obs;
}

double varF(int N, int n, int **freq, double *F, double A, double B, double t[], double h, double u)
{
  int			i;
  double		t1,t2,t3,sum;


  sum=0;

  for (i=1;i<=n;i++)
  {
    t1=(u-t[i])/h;
    t2=(u+t[i]-2*A)/h;
    t3=(2*B-u-t[i])/h;

    sum += SQR(K(t1)-K(t2)-K(t3))*(SQR(F[i]-1)*freq[1][i]+SQR(F[i])*freq[0][i]);
  }

  sum = sum/(N*N);

  return sum;
}



//////////////////////////////////////////////////////////////////////////////////
// Triweight kernel
//
// Export: No
// Date Created: 07.07.2017
// Created by: Kim Hendrickx
//////////////////////////////////////////////////////////////////////////////////
double K(double x)
{
    double u,y;

    u=x*x;

    if (u<=1)
        y=(35.0/32)*pow(1-u,3);
    else
        y=0.0;

    return y;
}
