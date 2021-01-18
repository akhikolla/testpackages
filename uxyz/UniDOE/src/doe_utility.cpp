#include <stdio.h>
#include <stdlib.h>
#include "doe_utility.h"
#include "doe_random.h"
#include "doe_Matrix.h"
#include "doe_Index.h"
#include <math.h>

double **range(int type,double **x,int nsamp,int nv);

double **normalize(int type,double **x,int nsamp,int nv)
{
	int i,j;
	double **rang,scale,ep;
	rang=range(type,x,nsamp,nv);
	for(j=0;j<nv;j++)
	{
		scale=rang[1][j]-rang[0][j];
		ep=EPS*(ABS(rang[1][j])+ABS(rang[0][j]));
		for(i=0;i<nsamp;i++)
		{
			if(scale>ep) x[i][j]=(x[i][j]-rang[0][j])/scale;
			else x[i][j]=0.5;
		}
	}
	return(rang);
}

void unnormalize(double **x,double **rang,int nsamp,int nv)
{
	int i,j;
	double scale;
	for(j=0;j<nv;j++)
	{
		scale=rang[1][j]-rang[0][j];
		for(i=0;i<nsamp;i++) x[i][j]=x[i][j]*scale+rang[0][j];
	}
}

double **range(int type,double **x,int nsamp,int nv)
{
	int i,j,*idx,level;
	double **rang,scale,*tempx;
	rang=NewDMatrix(2,nv);
	for(j=0;j<nv;j++)
	{
		rang[0][j]=rang[1][j]=x[0][j];
		for (i=1;i<nsamp;i++)
		{
			if(x[i][j]<rang[0][j]) rang[0][j]=x[i][j];
			else if(x[i][j]>rang[1][j]) rang[1][j]=x[i][j];
		}
	}
	if(type==2)
	{
		idx=NewIVector(nsamp);
		tempx=NewDVector(nsamp);
		for(j=0;j<nv;j++)
		{
			for(i=0;i<nsamp;i++) tempx[i]=x[i][j];
			indexx1(nsamp,tempx,(unsigned int*)idx);
            for(i=1,level=1;i<nsamp;i++)
                if(fabs(x[idx[i]][j]-x[idx[i-1]][j])>EPS2) level++;
			scale=0.5*(rang[1][j]-rang[0][j])/(level-1);
			rang[0][j]-=scale;
			rang[1][j]+=scale;
		}
		FreeVector(idx);
		FreeVector(tempx);
	}
	return(rang);
}

int mymin(int *v,int n)
{
	int i,minv;
	minv=v[0];
	for(i=1;i<n;i++) if(v[i]<minv) minv=v[i];
	return(minv);
}

int iCheckValue(int min,int max,int def,int value)
{
	if(def>max) def=max;
	else if(def<min) def=min;
	if(value<min||value>max) return(def);
	else return(value);
}

unsigned int uCheckValue(unsigned int min,unsigned int max,unsigned int def,unsigned int value)
{
	if(def>max) def=max;
	else if(def<min) def=min;
	if(value<min||value>max) return(def);
	else return(value);
}


char cCheckValue(char min,char max,char def,char value)
{
	if(def>max) def=max;
	else if(def<min) def=min;
	if(value<min||value>max) return(def);
	else return(value);
}

double dCheckValue(double min,double max,double def,double value)
{
	if(def>max) def=max;
	else if(def<min) def=min;
	if(ABS(value-min)<=EPS2) return(min);
	else if(ABS(value-max)<=EPS2) return(max);
	else if(value<min||value>max) return(def);
	else return(value);
}

//randome permutation
void permute(int *x,int n)
{
	double *temp;
	int *indx,*xt;
	int i;
	if(!n||!x) return;
	temp=NewDVector(n);
	indx=NewIVector(n);
	xt=NewIVector(n);
	for(i=0;i<n;i++) xt[i]=x[i];
	for(i=0;i<n;i++) temp[i]=Random();
	indexx1(n, temp, (unsigned int*)indx);
	for(i=0;i<n;i++) x[i]=xt[indx[i]];
	FreeVector(temp);
	FreeVector(indx);
	FreeVector(xt);
}

double mult(double x,int p)
{
	char k[10];//p<1024
	int t,i,t1,j;
	double y,y1;
	t=p; i=0; y=x;
	while(t>1)
	{
		t1=(t>>1);
		k[i]=t-(t1<<1);
		i++;
		t=t1;
	}
	for(j=i-1;j>=0;j--)
	{
		y1=y*y;
		if(k[j]) y=y1*x;
		else y=y1;
	}
	//printf("In utility mult: x = %f, p = %d; y = %f\n",x,p,y);
	return(y);
}
