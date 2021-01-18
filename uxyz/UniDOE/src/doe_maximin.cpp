#include "doe_Matrix.h"
#include "doe_criteria.h"
#include <math.h>
#include "doe_maximin.h"
#include "doe_utility.h"

#include <stdio.h>
#include <stdlib.h>
#define DEF_PMM 20
#define DEBUG 0

static double *scale;
static char scaled,pd; // 0: L1 otherwise:L2
static int nv,np,nnew,nsamp,pmm;
static char set_cnt;

static double **x,**D,mmres,mmres1,maxmm,minmm;
static double *xsnap,**Dsnap1,**Dsnap2,mmressnap,mmres1snap;
static double **xf, **Df1,**Df2,mmresfull,mmres1full;
static double *Dtemp;

double maximin_set_pinf(double **x0);
double maximin_cp1_pinf(int ncol,int i1,int i2);
double maximin_cp_pinf(int ncol,int npm,int *idx1,int *idx2);
double maximin_cp_set_pinf(int ncol,int npm,int *idx1,int *idx2);
double maximin_pm_pinf(int ncol, int npm,int *idx1,int *idx2);
double maximin_pm_set_pinf(int ncol, int npm,int *idx1,int *idx2);
double maximin_eval_pinf(double **x0);


// pars[1]:pmm, pars[0]: pd
//pars[2]: scale x'=scale*x  (not like entropy: x'=sqrt(theta)*x
/* The structure of D: upper matrix is L1 d or L2 d^2; lower matrix is d^(-p) or d^(-2p)
diagonal elements is used for temp mem*/
void create_maximin(double **x0, int nnew1, int np1,int nv1,CRITOPT *critopt)
{
	int i;
	nv=nv1;
	nnew=nnew1;
	np=np1;
	nsamp=nnew+np;
	scale=NewDVector(nv);
	pmm=DEF_PMM;
	pd=1;

	if(critopt->npars[0]>0) pmm=iCheckValue(0,100,DEF_PMM,(int)(critopt->pars[0][0]+0.01));
	if(critopt->npars[1]>0) pd=iCheckValue(0,1,1,(int)(critopt->pars[1][0]+0.01)-1);
//	for(i=0,maxmm=1.0;i<pmm;i++) maxmm*=10; //maxmm=1.0e6 if pmm=6
	minmm=pow(MINIDOUBLE,1.0/pmm);
	maxmm=1.0/minmm;

	#if DEBUG
		// printf("In maximin.cpp create_maximin: pmm = %d; pd = %d; scaled = %d\n",pmm,pd,scaled);
		// printf("In maximin.cpp create_maximin: minmm = %f; maxmm = %f\n",minmm,maxmm);
	#endif

	if(pd) pmm=(pmm+1)/2;

	scale=critopt->scale;
	scaled=0;
	for(i=0;i<nv;i++) if(ABS(scale[i]-1)>EPS) {scaled=1; break;}

	x=NewDMatrix(nsamp,nv);
	D=NewDMatrix(nsamp,nsamp);
	xsnap=NewDVector(nnew);
	Dsnap1=NewDMatrix(nnew,nsamp);
	if(pmm) Dsnap2=NewDMatrix(np,nnew);

	xf=NewDMatrix(nnew,nv);
	Df1=NewDMatrix(nnew,nsamp);
	if(pmm) Df2=NewDMatrix(np,nnew);

    Dtemp=NewDVector(2*nsamp);
	maximin_set(x0);
}

void free_maximin()
{
	FreeDMatrix(x);
	FreeDMatrix(D);
	FreeVector(xsnap);
	FreeDMatrix(Dsnap1);
	if (pmm) FreeDMatrix(Dsnap2);
	FreeDMatrix(xf);
	FreeDMatrix(Df1);
	if (pmm) FreeDMatrix(Df2);
	FreeVector(Dtemp);
}


double maximin_set(double **x0)
{
	double d1,dt;
	int i,j,k;
	set_cnt=0;

	if(!pmm) return(maximin_set_pinf(x0));
	if(x0!=NULL)
	{
		for(i=0;i<nsamp;i++)
		 {
			for(j=0;j<nv;j++)
				{
					if(scaled) x[i][j]=scale[j]*x0[i][j];
					else x[i][j]=x0[i][j];
				}
		 }
	}
	for(i=0;i<nsamp;i++)
	{
		for(j=i+1;j<nsamp;j++)
		{
			dt=0;
			for(k=0;k<nv;k++)
			{
				d1=x[i][k]-x[j][k];
				if(pd) dt+=d1*d1;
				else dt+=ABS(d1);
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
	if(pd) mmres=pow(mmres1,1.0/pmm/2.0);
	else mmres=pow(mmres1,1.0/pmm);
	return(mmres);
}

double maximin_cp_set(int ncol,int ncp,int *idx1,int *idx2)
{
    int i1,i2,i,j,id[2],ii1,ii2,m,counter=0;
	double dd,d1,d2,dp,dt,diff1,diff2,temp;

	if(!pmm) return(maximin_cp_set_pinf(ncol,ncp,idx1,idx2));

	diff1=diff2=0;
	for(j=0;j<ncp;j++)
	{
		i1=id[0]=idx1[j]; i2=id[1]=idx2[j];
		if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2)
		{
			temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
			counter++; continue;
		}
		for(i=0;i<nsamp;i++)
		{
			if(i != i1 && i != i2)
			{
				d1 = x[i][ncol]-x[i1][ncol];
				d2=  x[i][ncol]-x[i2][ncol];
				if(pd) dd=d2*d2-d1*d1;
				else dd=ABS(d2)-ABS(d1);
				for(m=0;m<2;m++)
				{
					if(i<id[m]) {ii1=i; ii2=id[m];}
					else {ii1=id[m]; ii2=i;}
					if(m) dd=-dd;
					dp=dt=MAX((D[ii1][ii2]+dd),0);
					D[ii1][ii2]=dp;  //change to new value
					if(dp<minmm) dp=MINIDOUBLE;
					else dp=mult(dt,pmm);//for(k=0;k<pmm-1;k++) dp*=dt;
					dp=1/dp;
					diff1+=dp;
					diff2+=D[ii2][ii1];
					D[ii2][ii1]=dp;    //change to new value
				}
			}
		}
		temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
	}

	if(counter==ncp) return(mmres);
	diff1=MIN(diff1,MAXDOUBLE);
    diff2=MIN(diff2,MAXDOUBLE);
	mmres1=(mmres1+diff1)-diff2;
	if(set_cnt==5||mmres1/diff2<1.0e-8) //avoid rounding errors,including mmres1<0
	{
		set_cnt=0;
		mmres1=0; for(i=0;i<nsamp;i++) for(j=0;j<i;j++) mmres1+=D[i][j];
	}
	else set_cnt++;
	mmres1=MIN(mmres1,MAXDOUBLE);

	if(pd) mmres=pow(mmres1,1.0/pmm/2.0);
	else mmres=pow(mmres1,1.0/pmm);
	return(mmres);
}

double maximin_cp(int ncol,int ncp,int *idx1,int *idx2)
{
	int i1,i2,i,j,ncp2;
	char tmp_set;
	double newmmres,**tempD1,**tempD2,tmpmmres,tmpmmres1,temp;
	if(!pmm) return(maximin_cp_pinf(ncol,ncp,idx1,idx2));
	ncp2=ncp*2;
	tempD1=NewDMatrix(ncp2,nsamp);
	tempD2=NewDMatrix(ncp2,nsamp);
	//save
	tmpmmres=mmres; tmpmmres1=mmres1;
	for(j=0;j<ncp;j++)
	{
		i1=idx1[j]; i2=idx2[j];
		if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2) continue;
		for(i=0;i<nsamp;i++)
		{
			tempD1[j][i]=D[i][i1]; tempD1[ncp2-j-1][i]=D[i1][i];
			tempD2[j][i]=D[i][i2]; tempD2[ncp2-j-1][i]=D[i2][i];
		}
	}

	tmp_set=set_cnt;
	newmmres=maximin_cp_set(ncol,ncp,idx1,idx2);
	set_cnt=tmp_set;

	//restore
	for(j=0;j<ncp;j++)
	{
		i1=idx1[j]; i2=idx2[j];
		if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2) continue;
		for(i=0;i<nsamp;i++)
		{
			D[i][i1]=tempD1[j][i]; D[i1][i]=tempD1[ncp2-j-1][i];
			D[i][i2]=tempD2[j][i]; D[i2][i]=tempD2[ncp2-j-1][i];
		}
		temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
	}
	mmres=tmpmmres; mmres1=tmpmmres1;
	FreeDMatrix(tempD1);
	FreeDMatrix(tempD2);
	return(newmmres);
}


double maximin_cp1(int ncol,int i1,int i2)
{
	int i,j,id[2],m,ii1,ii2;
	double dd,d1,d2,diff1,diff2,newmmres,newmmres1,dp,dt;
	if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2) return(mmres);
	if(!pmm) return(maximin_cp1_pinf(ncol,i1,i2));
	/*diff=0;*/
	diff1=diff2=0;

	id[0]=i1; id[1]=i2;
	for(i=0,j=0;i<nsamp;i++)
	{
		if(i != i1 && i != i2)
		{
			d1 = x[i][ncol]-x[i1][ncol];
			d2=  x[i][ncol]-x[i2][ncol];
			if(pd) dd=d2*d2-d1*d1;
			else dd=ABS(d2)-ABS(d1);
			for(m=0;m<2;m++)
			{
				if(i<id[m]) {ii1=i; ii2=id[m];}
				else {ii1=id[m]; ii2=i;}
				if(m) dd=-dd;
				dp=dt=MAX((D[ii1][ii2]+dd),0);
				if(dp<minmm) dp=MINIDOUBLE;
				else dp=mult(dt,pmm);//for(k=0;k<pmm-1;k++) dp*=dt;
				dp=1/dp;
				diff1+=dp;
				diff2+=D[ii2][ii1];
				Dtemp[j]=D[ii2][ii1]; j++; // save a copy
				D[ii2][ii1]=dp;  //change to new value
			}
		}
	}
	diff1=MIN(diff1,MAXDOUBLE);
    diff2=MIN(diff2,MAXDOUBLE);
	newmmres1=(mmres1+diff1)-diff2;
	if(newmmres1/diff2<1.0e-8) //avoid rounding errors,including mmwmmres1<0
	{
		newmmres1=0; for(i=0;i<nsamp;i++) for(j=0;j<i;j++) newmmres1+=D[i][j];
	}
	newmmres1=MIN(newmmres1,MAXDOUBLE);

	if(pd) newmmres=pow(newmmres1,1.0/pmm/2.0);
	else newmmres=pow(newmmres1,1.0/pmm);
	for(i=0,j=0;i<nsamp;i++)
	{
		if(i != i1 && i != i2)
		{
			for(m=0;m<2;m++)
			{
				if(i<id[m]) {ii1=i; ii2=id[m];}
				else {ii1=id[m]; ii2=i;}
				D[ii2][ii1]=Dtemp[j]; j++; //restore
			}
		}
	}
	return(newmmres);
}


double maximin()
{
	return(mmres);
}

double **maximin_x(double **xnew)
{
	int i,j;
	if(xnew==NULL) return(x);

	for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
	{
		xnew[i][j]=x[i][j];
		if(scaled) xnew[i][j]/=scale[j];
	}
	return(xnew);
}

void maximin_snap(int ncol)
{
	int i,j;
	mmressnap=mmres;
	mmres1snap=mmres1;
	if(pmm)
		for(i=np;i<nsamp;i++)
		{
			xsnap[i-np]=x[i][ncol];
			for(j=0;j<nsamp;j++) Dsnap1[i-np][j]=D[i][j];
			for(j=0;j<np;j++) Dsnap2[j][i-np]=D[j][i];
		}
	else
		for(i=np;i<nsamp;i++)
		{
			xsnap[i-np]=x[i][ncol];
			for(j=0;j<i;j++) Dsnap1[i-np][j]=D[j][i];
		}
}

void maximin_reset(int ncol)
{
	int i,j;
	mmres=mmressnap;
	mmres1=mmres1snap;
	if(pmm)
		for(i=np;i<nsamp;i++)
		{
			x[i][ncol]=xsnap[i-np];
			for(j=0;j<nsamp;j++) D[i][j]=Dsnap1[i-np][j];
			for(j=0;j<np;j++) D[j][i]=Dsnap2[j][i-np];
		}
	else
		for(i=np;i<nsamp;i++)
		{
			x[i][ncol]=xsnap[i-np];
			//this is necessary D[j][i]=D[i][j] for multiple pairs exchange
			for(j=0;j<i;j++) D[i][j]=D[j][i]=Dsnap1[i-np][j];
		}
}

//when pmm is inf
double maximin_set_pinf(double **x0)
{
	double d1,dt,newmmres;
	int i,j,k;
	if(x0!=NULL)
	{
		for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
		{
			if(scaled) x[i][j]=scale[j]*x0[i][j];
			else x[i][j]=x0[i][j];
		}
	}

	newmmres=MAXDOUBLE;

	for(i=0;i<nsamp;i++)
	{
		for(j=i+1;j<nsamp;j++)
		{
			dt=0;
			for(k=0;k<nv;k++)
			{
				d1=x[i][k]-x[j][k];
				if(pd) dt+=d1*d1;
				else dt+=ABS(d1);
			}
			D[j][i]=D[i][j]=dt;
		}
	}

	for(i=0;i<np;i++) for(j=i+1;j<np;j++) if(D[i][j]<newmmres) newmmres=D[i][j];
	mmres1=-newmmres;
	for(j=np;j<nsamp;j++)
		for(i=0;i<j;i++) if(D[i][j]<newmmres) newmmres=D[i][j];
	mmres=-newmmres;
	return(mmres);
}

//when pmm is inf
double maximin_cp_set_pinf(int ncol,int ncp,int *idx1,int *idx2)
{
	int i1,i2,i,j,counter=0;
	double dd,d1,d2,temp,newmmres;

	for(j=0;j<ncp;j++)
	{
		i1=idx1[j]; i2=idx2[j];
		if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2)
		{
			temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
			counter++; continue;
		}

		for(i=0;i<nsamp;i++)
		{
			if(i != i1 && i != i2)
			{
				d1 = x[i][ncol]-x[i1][ncol];
				d2=  x[i][ncol]-x[i2][ncol];
				if(pd) dd=d2*d2-d1*d1;
				else dd=ABS(d2)-ABS(d1);

				if(i<i1) {D[i][i1]+=dd; D[i1][i]=D[i][i1];}
				else {D[i1][i]+=dd; D[i][i1]=D[i1][i];}
				if(i<i2) {D[i][i2]-=dd; D[i2][i]=D[i][i2];}
				else {D[i2][i]-=dd; D[i][i2]=D[i2][i];}
			}
		}
		temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
	}
	if(counter==ncp) return(mmres);
	newmmres=-mmres1;
	for(j=np;j<nsamp;j++)
		for(i=0;i<j;i++) if(D[i][j]<newmmres) newmmres=D[i][j];
	mmres=-newmmres;
	return(mmres);
}

double maximin_cp_pinf(int ncol,int ncp,int *idx1,int *idx2)
{
	int i1,i2,i,j,counter=0;
	double dd,d1,d2,newmmres,temp;

	for(j=0;j<ncp;j++)
	{
		i1=idx1[j]; i2=idx2[j];
		if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2) {counter++; continue;}

		for(i=0;i<nsamp;i++)
		{
			if(i != i1 && i != i2)
			{
				d1 = x[i][ncol]-x[i1][ncol];
				d2=  x[i][ncol]-x[i2][ncol];
				if(pd) dd=d2*d2-d1*d1;
				else dd=ABS(d2)-ABS(d1);

				if(i<i1) D[i][i1]+=dd;
				else D[i1][i]+=dd;
				if(i<i2) D[i][i2]-=dd;
				else D[i2][i]-=dd;
			}
		}
		temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
	}

	if(counter==ncp) return(mmres);

	newmmres=-mmres1;
	for(j=np;j<nsamp;j++)
		for(i=0;i<j;i++) if(D[i][j]<newmmres) newmmres=D[i][j];

	for(j=0;j<ncp;j++)
	{
		i1=idx1[j]; i2=idx2[j];
		if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2) continue;
		for(i=0;i<nsamp;i++)
		{
			if(i != i1 && i != i2)
			{
				if(i<i1) D[i][i1]=D[i1][i];
				else D[i1][i]=D[i][i1];
				if(i<i2) D[i][i2]=D[i2][i];
				else D[i2][i]=D[i][i2];
			}
		}
		temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
	}
	return(-newmmres);
}


double maximin_cp1_pinf(int ncol,int i1,int i2)
{
	int i,j;
	double dd,d1,d2,newmmres;
	if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2) return(mmres);
	newmmres=MAXDOUBLE;
	for(i=0;i<nsamp;i++)
	{
		if(i != i1 && i != i2)
		{
			d1 = x[i][ncol]-x[i1][ncol];
			d2=  x[i][ncol]-x[i2][ncol];
			if(pd) dd=d2*d2-d1*d1;
			else dd=ABS(d2)-ABS(d1);

			if(i<i1) D[i][i1]+=dd;
			else D[i1][i]+=dd;

			if(i<i2) D[i][i2]-=dd;
			else D[i2][i]-=dd;
		}
	}
	newmmres=MAXDOUBLE;
	for(j=np;j<nsamp;j++)
		for(i=0;i<j;i++) if(D[i][j]<newmmres) newmmres=D[i][j];
	newmmres=(newmmres>-mmres1?-mmres1:newmmres);
	for(i=0;i<nsamp;i++)
	{
		if(i != i1 && i != i2)
		{
			if(i<i1) D[i][i1]=D[i1][i];
			else D[i1][i]=D[i][i1];

			if(i<i2) D[i][i2]=D[i2][i];
			else D[i2][i]=D[i][i2];
		}
	}

	return(-newmmres);
}


double maximin_eval(double **x0)
{
	double d1,dt,**x1,**D1,newmmres;
	int i,j,k;
	if(!pmm) return(maximin_eval_pinf(x0));

	if(x0==NULL) x0=x;
	x1=NewDMatrix(nsamp,nv);
	D1=NewDMatrix(nsamp,nsamp);
	for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
	{
		if(x0!=x && scaled) x1[i][j]=scale[j]*x0[i][j];
		else x1[i][j]=x0[i][j];
	}
	for(i=0;i<nsamp;i++)
	{
		for(j=i+1;j<nsamp;j++)
		{
			dt=0;
			for(k=0;k<nv;k++)
			{
				d1=x1[i][k]-x1[j][k];
				if(pd) dt+=d1*d1;
				else dt+=ABS(d1);
			}
			D1[i][j]=dt;
			if(dt<minmm) D1[j][i]=MAXDOUBLE;
			else
			{
				D1[j][i]=dt;
				D1[j][i]=mult(dt,pmm);//for(k=0;k<pmm-1;k++) D1[j][i]*=dt;
				D1[j][i]=1/D1[j][i];
			}
		}
	}
	newmmres=0;
	for(i=0;i<nsamp;i++)
		for(j=0;j<i;j++) newmmres+=D1[i][j];
	if(pd) newmmres=pow(newmmres,1.0/pmm/2.0);
	else newmmres=pow(newmmres,1.0/pmm);
	FreeDMatrix(D1);
	FreeDMatrix(x1);
	return(newmmres);
}


//when pmm is inf
double maximin_eval_pinf(double **x0)
{
	double d1,dt,newmmres,**x1,**D1;
	int i,j,k;

	x1=NewDMatrix(nsamp,nv);
	D1=NewDMatrix(nsamp,nsamp);
	if(x0==NULL) x0=x;

	for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
	{
		if(x0==NULL && scaled) x1[i][j]=scale[j]*x0[i][j];
		else x1[i][j]=x0[i][j];
	}

	newmmres=MAXDOUBLE;

	for(i=0;i<nsamp;i++)
	{
		for(j=i+1;j<nsamp;j++)
		{
			dt=0;
			for(k=0;k<nv;k++)
			{
				d1=x1[i][k]-x1[j][k];
				if(pd) dt+=d1*d1;
				else dt+=ABS(d1);
			}
			D1[j][i]=D1[i][j]=dt;
		}
	}

	for(i=0;i<np;i++) for(j=i+1;j<np;j++) if(D1[i][j]<newmmres) newmmres=D1[i][j];
	for(j=np;j<nsamp;j++) for(i=0;i<j;i++) if(D1[i][j]<newmmres) newmmres=D1[i][j];
	newmmres=-newmmres;
	FreeDMatrix(D1);
	FreeDMatrix(x1);
	return(newmmres);
}

double maximin_pm(int ncol, int npm,int *idx1,int *idx2)
{
	int i1,i,j,npm2;
	double newmmres,**tempD,tmpmmres,tmpmmres1;
	char tmp_set;

	if(!pmm) return(maximin_pm_pinf(ncol,npm,idx1,idx2));
	for(i=0;i<npm;i++) D[i][i]=x[idx2[i]][ncol]; //D[i][i] is used to save x

	npm2=npm*2;
	tempD=NewDMatrix(npm2,nsamp);
	//save
	tmpmmres=mmres; tmpmmres1=mmres1;
	for(j=0;j<npm;j++)
	{
		i1=idx1[j];
		for(i=0;i<nsamp;i++)
		{
			tempD[j][i]=D[i][i1]; tempD[npm2-j-1][i]=D[i1][i];
		}
	}

	tmp_set=set_cnt;
    newmmres=maximin_pm_set(ncol,npm,idx1,idx2);
    set_cnt=tmp_set;

	//restore
	mmres=tmpmmres; mmres1=tmpmmres1;
	for(i=0;i<npm;i++) x[idx2[i]][ncol]=D[i][i];
	for(j=0;j<npm;j++)
	{
		i1=idx1[j];
		for(i=0;i<nsamp;i++)
		{
			D[i][i1]=tempD[j][i]; D[i1][i]=tempD[npm2-j-1][i];
		}
	}
	FreeDMatrix(tempD);
	return(newmmres);
}


double maximin_pm_set(int ncol, int npm,int *idx1,int *idx2)
{
	int i1,i,j,ii1,ii2,counter=0;
	double diff1,diff2,xj,dp,dt,d1,d2,dd;

	if(!pmm) return(maximin_pm_set_pinf(ncol,npm,idx1,idx2));

	for(i=0;i<npm;i++) D[i][i]=x[idx2[i]][ncol]; //D[i][i] is used to save x

	diff1=diff2=0;
	for(j=0;j<npm;j++)
	{
		i1=idx1[j];  xj=D[j][j];
		if(ABS(x[i1][ncol]-xj)<EPS2) {x[i1][ncol]=xj; counter++; continue;}
		for(i=0;i<nsamp;i++)
		{
			if(i != i1)
			{
				d1 = x[i][ncol]-x[i1][ncol];
				d2=  x[i][ncol]-xj;
				if(pd) dd=d2*d2-d1*d1;
				else dd=ABS(d2)-ABS(d1);
				if(i<i1) {ii1=i; ii2=i1;}
				else {ii1=i1; ii2=i;}
				dp=dt=MAX((D[ii1][ii2]+dd),0);
				D[ii1][ii2]=dp;  //change to new value
				if(dp<minmm) dp=MINIDOUBLE;
				dp=mult(dt,pmm);//for(k=0;k<pmm-1;k++) dp*=dt;
				dp=1/dp;
				diff1+=dp;
				diff2+=D[ii2][ii1];
				D[ii2][ii1]=dp;    //change to new value
			}
		}
		x[i1][ncol]=xj;
	}

	if(counter==npm) return(mmres);

	diff1=MIN(diff1,MAXDOUBLE);
    diff2=MIN(diff2,MAXDOUBLE);
	mmres1=(mmres1+diff1)-diff2;
	if(set_cnt==5||mmres1/diff2<1.0e-8) //avoid rounding errors,including mmres1<0
	{
		mmres1=0; for(i=0;i<nsamp;i++) for(j=0;j<i;j++) mmres1+=D[i][j];
		set_cnt=0;
	}
	else set_cnt++;

	mmres1=MIN(mmres1,MAXDOUBLE);

	if(pd) mmres=pow(mmres1,1.0/pmm/2.0);
	else mmres=pow(mmres1,1.0/pmm);
	return(mmres);
}

//when pmm is inf
double maximin_pm_set_pinf(int ncol,int npm,int *idx1,int *idx2)
{
	int i1,i,j,counter=0;
	double dd,d1,d2,newmmres,xj;

 	for(i=0;i<npm;i++) D[i][i]=x[idx2[i]][ncol]; //D[i][i] is used to save x
	for(j=0;j<npm;j++)
	{
		i1=idx1[j];  xj=D[j][j];
		if(ABS(x[i1][ncol]-xj)<EPS2) {x[i1][ncol]=xj; counter++; continue;}

		for(i=0;i<nsamp;i++)
		{
			if(i != i1)
			{
				d1 = x[i][ncol]-x[i1][ncol];
				d2=  x[i][ncol]-xj;
				if(pd) dd=d2*d2-d1*d1;
				else dd=ABS(d2)-ABS(d1);

				if(i<i1) {D[i][i1]+=dd; D[i1][i]=D[i][i1];}
				else {D[i1][i]+=dd; D[i][i1]=D[i1][i];}
			}
		}
		x[i1][ncol]=xj;
	}

	if(counter==npm) return(mmres);
	newmmres=-mmres1;
	for(j=np;j<nsamp;j++)
		for(i=0;i<j;i++) if(D[i][j]<newmmres) newmmres=D[i][j];
	mmres=-newmmres;
	return(mmres);
}

double maximin_pm_pinf(int ncol,int npm,int *idx1,int *idx2)
{
	int i1,i,j,counter=0;
	double dd,d1,d2,newmmres,xj;

	for(i=0;i<npm;i++) D[i][i]=x[idx2[i]][ncol]; //D[i][i] is used to save x

	for(j=0;j<npm;j++)
	{
		i1=idx1[j]; xj=D[j][j];
		if(ABS(x[i1][ncol]-xj)<EPS2) {counter++; continue;}

		for(i=0;i<nsamp;i++)
		{
			if(i != i1)
			{
				d1 = x[i][ncol]-x[i1][ncol];
				d2=  x[i][ncol]-xj;
				if(pd) dd=d2*d2-d1*d1;
				else dd=ABS(d2)-ABS(d1);

				if(i<i1) D[i][i1]+=dd;
				else D[i1][i]+=dd;
			}
		}
		x[i1][ncol]=xj;
	}

	if(counter==npm) return(mmres);

	newmmres=-mmres1;
	for(j=np;j<nsamp;j++)
		for(i=0;i<j;i++) if(D[i][j]<newmmres) newmmres=D[i][j];

	for(i=0;i<npm;i++) x[idx2[i]][ncol]=D[i][i];

	for(j=0;j<npm;j++)
	{
		i1=idx1[j];
		if(ABS(x[i1][ncol]-D[j][j])<EPS2) continue;
		for(i=0;i<nsamp;i++)
		{
			if(i<i1) D[i][i1]=D[i1][i];
			else if(i>i1) D[i1][i]=D[i][i1];
		}
	}
	return(-newmmres);
}

void maximin_full_snap()
{
	int i,j;
	mmresfull=mmres;
	mmres1full=mmres1;
	if(pmm)
		for(i=np;i<nsamp;i++)
		{
			for(j=0;j<nv;j++) xf[i-np][j]=x[i][j];
			for(j=0;j<nsamp;j++) Df1[i-np][j]=D[i][j];
			for(j=0;j<np;j++) Df2[j][i-np]=D[j][i];
		}
	else
		for(i=np;i<nsamp;i++)
		{
			for(j=0;j<nv;j++) xf[i-np][j]=x[i][j];
			for(j=0;j<i;j++) Df1[i-np][j]=D[j][i];
		}
}

void maximin_full_reset()
{
	int i,j;
	mmres=mmresfull;
	mmres1=mmres1full;
	if(pmm)
		for(i=np;i<nsamp;i++)
		{
			for(j=0;j<nv;j++) x[i][j]=xf[i-np][j];
			for(j=0;j<nsamp;j++) D[i][j]=Df1[i-np][j];
			for(j=0;j<np;j++) D[j][i]=Df2[j][i-np];
		}
	else
		for(i=np;i<nsamp;i++)
		{
			for(j=0;j<nv;j++) x[i][j]=xf[i-np][j];
			//this is necessary D[j][i]=D[i][j] for multiple pairs exchange
			for(j=0;j<i;j++) D[i][j]=D[j][i]=Df1[i-np][j];
		}
}

void maximin_global_x(double **xnew)
{
	int i,j;
	for(j=0;j<nv;j++)
	{
		for(i=0;i<np;i++)  xnew[i][j]=x[i][j];
    	for(i=np;i<nsamp;i++)  xnew[i][j]=xf[i-np][j];
	}
	if(scaled) for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) xnew[i][j]/=scale[j];
}
