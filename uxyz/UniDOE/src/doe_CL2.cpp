#include "doe_Matrix.h"
#include <math.h>
#include "doe_criteria.h"
#include "doe_CL2.h"
#include "doe_utility.h"

#define DEF_SCALE 1.0
#define DEF_PMM 6
#define DEF_DIST 1

static double *scale;
static char scaled;
static int nv,nactive;
static int np;
static int nnew;
static int nsamp;

static double **x,**xc,**D,discr,**xf,*Df,discrfull;
static double *xsnap,*xcsnap,*Dsnap,discrsnap;


/* The structure of D (nsamp+1*nsamp+1):
  lower matrix is exactly a copy of upper matrix, except the nsamp column
*/
void create_discrcl2( double **x0,int nnew1, int np1,int nv1,CRITOPT *critopt)
{
	int i;
	nv=nv1;
	nnew=nnew1;
	np=np1;
	nsamp=nnew+np;

	scale=critopt->scale;
	scaled=0;
	for(i=0,nactive=0;i<nv;i++) if(scale[i]>EPS2) nactive++;
	for(i=0;i<nv;i++) if(ABS(scale[i]-1)>EPS) {scaled=1; break;}

	x=NewDMatrix(nsamp,nv);
	xc=NewDMatrix(nsamp,nv);
	D=NewDMatrix(nsamp+1,nsamp+1); //note: (nsamp+1)*(nsamp+1)
	xsnap=NewDVector(nnew);
	xcsnap=NewDVector(nnew);
	Dsnap=NewDVector((np+nsamp+2)*(nnew+1)/2);
	xf=NewDMatrix(nnew,nv);
	Df=NewDVector((np+nsamp+2)*(nnew+1)/2);
	discrcl2_set(x0);
}

void free_discrcl2()
{
	FreeDMatrix(x);
	FreeDMatrix(xc);
	FreeDMatrix(D);
	FreeVector(xsnap);
	FreeVector(xcsnap);
	FreeVector(Dsnap);
	FreeDMatrix(xf);
	FreeVector(Df);
}

double discrcl2_set(double **x0)
{
	double d1,d2;
	int i,j,k,nsamp2;
	if(x0!=NULL)
	{
		for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
		{
			if(scaled)
			{
				x[i][j]=scale[j]*x0[i][j]; 	xc[i][j]=scale[j]*(x0[i][j]-0.5);
			}
			else
			{
				x[i][j]=x0[i][j]; 	xc[i][j]=x0[i][j]-0.5;
			}
			xc[i][j]=ABS(xc[i][j]);
		}
	}
	nsamp2=nsamp*nsamp;
	D[nsamp][nsamp]=0;
	for(i=0;i<nsamp;i++)
	{
		for(k=0,d1=1,d2=1;k<nv;k++)
		{
			d1*=1+xc[i][k];
			d2*=(2+xc[i][k]-xc[i][k]*xc[i][k])/2;
		}
		D[i][i]=d1/nsamp2;
		D[nsamp][i]=D[i][nsamp]=-2*d2/nsamp; //D has nsamp+1 column
		for(j=i+1;j<nsamp;j++)
		{
			for(k=0,d1=1;k<nv;k++)
			{
				d2=x[i][k]-x[j][k];
				d1*=(2+xc[i][k]+xc[j][k]-ABS(d2))/2;
			}
			D[j][i]=D[i][j]=2*d1/nsamp2;
		}
	}
	for(i=0,discr=1;i<nactive;i++) discr*=(13.0/12.0);
	for(i=0;i<nsamp+1;i++) 	for(j=i;j<nsamp+1;j++)	discr+=D[i][j];
	return(discr);
}

double discrcl2_cp_set(int ncol,int ncp,int *idx1,int *idx2)
{
	int i1,i2,i,j;
	double d1,d2,dd,da1,da2,db1,db2,t1,t2,tt,temp,diff;
	diff=0;
	for(j=0;j<ncp;j++)
	{
		i1=idx1[j]; i2=idx2[j];
		if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2)
		{
			temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
			temp=xc[i1][ncol]; xc[i1][ncol]=xc[i2][ncol]; xc[i2][ncol]=temp;
			continue;
		}
		t1=1+xc[i1][ncol]; t2=1+xc[i2][ncol];
		tt=t2/t1; //alpha(i1,i2,k)
		d1=2-xc[i1][ncol]; d2=2-xc[i2][ncol];
		dd=d2/d1; //beta(i1,i2,k)
		da1=D[i1][i1]*tt; da2=D[i2][i2]/tt;
        db1=D[i1][nsamp]*dd*tt; db2=D[i2][nsamp]/dd/tt;

		diff+=da1-D[i1][i1]+da2-D[i2][i2]+db1-D[i1][nsamp]+db2-D[i2][nsamp];
		D[i1][i1]=da1; D[i2][i2]=da2;
        D[nsamp][i1]=D[i1][nsamp]=db1; D[nsamp][i2]=D[i2][nsamp]=db2;

		for(i=0;i<nsamp;i++)
		{
			if(i != i1 && i != i2)
			{
				t1=x[i][ncol]-x[i1][ncol];
				t2=x[i][ncol]-x[i2][ncol];
				d1 = 2+xc[i1][ncol]+xc[i][ncol]-ABS(t1);
				d2 = 2+xc[i2][ncol]+xc[i][ncol]-ABS(t2);
                dd=d2/d1;

				da1=D[i][i1]*dd;
				diff+=da1-D[i][i1];
				D[i1][i]=D[i][i1]=da1;  //change to new value

				da1=D[i][i2]/dd;
				diff+=da1-D[i][i2];
				D[i2][i]=D[i][i2]=da1;  //change to new value
			}
		}
		temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
		temp=xc[i1][ncol]; xc[i1][ncol]=xc[i2][ncol]; xc[i2][ncol]=temp;
	}

	discr+=diff;
	return(discr);
}


double discrcl2_cp(int ncol,int ncp,int *idx1,int *idx2)
{
	int i1,i2,i,j;
	double d1,d2,dd,da1,da2,db1,db2,t1,t2,tt,diff,newdiscr,temp;
	diff=0;
	for(j=0;j<ncp;j++)
	{
		i1=idx1[j]; i2=idx2[j];
		if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2) continue;

		t1=1+xc[i1][ncol]; t2=1+xc[i2][ncol];
		tt=t2/t1;
		d1=2-xc[i1][ncol]; d2=2-xc[i2][ncol];
		dd=d2/d1;
		da1=D[i1][i1]*tt; da2=D[i2][i2]/tt;
        db1=D[i1][nsamp]*dd*tt; db2=D[i2][nsamp]/dd/tt;

		diff+=da1-D[i1][i1]+da2-D[i2][i2]+db1-D[i1][nsamp]+db2-D[i2][nsamp];

		for(i=0;i<nsamp;i++)
		{
			if(i != i1 && i != i2)
			{
				t1=x[i][ncol]-x[i1][ncol];
				t2=x[i][ncol]-x[i2][ncol];
				d1 = 2+xc[i1][ncol]+xc[i][ncol]-ABS(t1);
				d2 = 2+xc[i2][ncol]+xc[i][ncol]-ABS(t2);
                dd=d2/d1;
				if(i<i1)
				{
					da1=D[i][i1]*dd; diff+=da1-D[i][i1];
					D[i][i1]=da1;  //change to new value
				}
				else
				{
					da1=D[i1][i]*dd; diff+=da1-D[i1][i];
					D[i1][i]=da1;  //change to new value
				}
				if(i<i2)
				{
					da1=D[i][i2]/dd; diff+=da1-D[i][i2];
					D[i][i2]=da1;  //change to new value
				}
				else
				{
					da1=D[i2][i]/dd; diff+=da1-D[i2][i];
					D[i2][i]=da1;  //change to new value
				}
			}
		}
		temp=x[i1][ncol]; x[i1][ncol]=x[i2][ncol]; x[i2][ncol]=temp;
		temp=xc[i1][ncol]; xc[i1][ncol]=xc[i2][ncol]; xc[i2][ncol]=temp;
	}

	//restore changed value
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
		temp=xc[i1][ncol]; xc[i1][ncol]=xc[i2][ncol]; xc[i2][ncol]=temp;
	}
	newdiscr=discr+diff;
	return(newdiscr);
}

double discrcl2_cp1(int ncol,int i1,int i2)
{
	int i;
	double d1,d2,dd,da1,da2,db1,db2,t1,t2,tt,diff,newdiscr;
    if(ABS(x[i1][ncol]-x[i2][ncol])<EPS2) return(discr);

	diff=0;
	t1=1+xc[i1][ncol]; t2=1+xc[i2][ncol];
	tt=t2/t1;
	d1=2-xc[i1][ncol]; d2=2-xc[i2][ncol];
	dd=d2/d1;
	da1=D[i1][i1]*tt; da2=D[i2][i2]/tt;
    db1=D[i1][nsamp]*dd*tt; db2=D[i2][nsamp]/dd/tt;

	diff+=da1-D[i1][i1]+da2-D[i2][i2]+db1-D[i1][nsamp]+db2-D[i2][nsamp];

	for(i=0;i<nsamp;i++)
	{
		if(i != i1 && i != i2)
		{
			t1=x[i][ncol]-x[i1][ncol];
			t2=x[i][ncol]-x[i2][ncol];
			d1 = 2+xc[i1][ncol]+xc[i][ncol]-ABS(t1);
			d2 = 2+xc[i2][ncol]+xc[i][ncol]-ABS(t2);
            dd=d2/d1;
			if(i<i1) diff+=D[i][i1]*dd-D[i][i1];
			else diff+=D[i1][i]*dd-D[i1][i];
			if(i<i2) diff+=D[i][i2]/dd-D[i][i2];
			else diff+=D[i2][i]/dd-D[i2][i];
		}
	}
	newdiscr=discr+diff;
	return(newdiscr);
}

double discrcl2_pm(int ncol, int npm,int *idx1,int *idx2)
{
	int i1,i2,i,j;
	double d1,d2,dd,da1,db1,t1,t2,tt,diff,newdiscr,*xt,*xct;
	xt=NewDVector(npm); xct=NewDVector(npm);
	for(i=0;i<npm;i++) {xt[i]=x[idx2[i]][ncol]; xct[i]=xc[idx2[i]][ncol];}
	diff=0;

	for(j=0;j<npm;j++)
	{
		i1=idx1[j];
		if(ABS(x[i1][ncol]-xt[j])<EPS2) continue;
		t1=1+xc[i1][ncol]; t2=1+xct[j]; tt=t2/t1;
		d1=2-xc[i1][ncol]; d2=2-xct[j]; dd=d2/d1;
		da1=D[i1][i1]*tt;  db1=D[i1][nsamp]*dd*tt;
		diff+=da1-D[i1][i1]+db1-D[i1][nsamp];

		for(i=0;i<nsamp;i++)
		{
			if(i != i1)
			{
				t1=x[i][ncol]-x[i1][ncol];
				t2=x[i][ncol]-xt[j];
				d1 = 2+xc[i1][ncol]+xc[i][ncol]-ABS(t1);
				d2 = 2+xct[j]+xc[i][ncol]-ABS(t2);
                dd=d2/d1;
				if(i<i1)
				{
					da1=D[i][i1]*dd;
					diff+=da1-D[i][i1];
					D[i][i1]=da1;  //change to new value
				}
				else
				{
					da1=D[i1][i]*dd;
					diff+=da1-D[i1][i];
					D[i1][i]=da1;  //change to new value
				}
			}
		}
		x[i1][ncol]=xt[j]; xc[i1][ncol]=xct[j];
	}

	//restore changed value
	for(j=0;j<npm;j++)
	{
		i2=idx2[j];
		if(ABS(x[i2][ncol]-xt[j])<EPS2) continue;
		for(i=0;i<nsamp;i++)
		{
			if(i<i2) D[i][i2]=D[i2][i];
			else if(i>i2) D[i2][i]=D[i][i2];
		}
		x[i2][ncol]=xt[j];  xc[i2][ncol]=xct[j];
	}
	newdiscr=discr+diff;
	FreeVector(xt);
	FreeVector(xct);

	return(newdiscr);
}


double discrcl2_pm_set(int ncol, int npm,int *idx1,int *idx2)
{
	int i1,i,j;
	double d1,d2,dd,da1,db1,t1,t2,tt,diff,*xt,*xct;
	xt=NewDVector(npm); xct=NewDVector(npm);
	for(i=0;i<npm;i++) {xt[i]=x[idx2[i]][ncol]; xct[i]=xc[idx2[i]][ncol];}
	diff=0;
	for(j=0;j<npm;j++)
	{
		i1=idx1[j];
		if(ABS(x[i1][ncol]-xt[j])<EPS2)
		{
			x[i1][ncol]=xt[j]; xc[i1][ncol]=xct[j];
			continue;
		}
		t1=1+xc[i1][ncol]; t2=1+xct[j]; tt=t2/t1;
		d1=2-xc[i1][ncol]; d2=2-xct[j]; dd=d2/d1;
		da1=D[i1][i1]*tt;  db1=D[i1][nsamp]*dd*tt;
		diff+=da1-D[i1][i1]+db1-D[i1][nsamp];
		D[i1][i1]=da1; D[nsamp][i1]=D[i1][nsamp]=db1;

		for(i=0;i<nsamp;i++)
		{
			if(i != i1)
			{
				t1=x[i][ncol]-x[i1][ncol];
				t2=x[i][ncol]-xt[j];
				d1 = 2+xc[i1][ncol]+xc[i][ncol]-ABS(t1);
				d2 = 2+xct[j]+xc[i][ncol]-ABS(t2);
                dd=d2/d1;
				da1=D[i][i1]*dd;
				diff+=da1-D[i][i1];
				D[i][i1]=D[i1][i]=da1;  //change to new value
			}
		}
		x[i1][ncol]=xt[j]; xc[i1][ncol]=xct[j];
	}
	discr+=diff;
	FreeVector(xt);
	FreeVector(xct);
	return(discr);
}


double discrcl2()
{
	return(discr);
}

double **discrcl2_x(double **xnew)
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

void discrcl2_snap(int ncol)
{
	int i,j,k;
	discrsnap=discr;
	for(i=np;i<nsamp;i++)
	{
		xsnap[i-np]=x[i][ncol];
		xcsnap[i-np]=xc[i][ncol];
	}
	for(k=0,i=np;i<nsamp+1;i++)
		for(j=0;j<=i;j++,k++) Dsnap[k]=D[i][j];

}

void discrcl2_reset(int ncol)
{
	int i,j,k;
	discr=discrsnap;
	for(i=np;i<nsamp;i++)
	{
		x[i][ncol]=xsnap[i-np];
		xc[i][ncol]=xcsnap[i-np];
	}
	for(k=0,i=np;i<nsamp+1;i++)
		//this is necessary D[j][i]=D[i][j] for multiple pairs exchange
		for(j=0;j<=i;j++,k++) D[j][i]=D[i][j]=Dsnap[k];
}


double discrcl2_eval(double **x0)
{
	double d1,d2,newdiscr,**x1,**xc1,**D1;
	int i,j,k,nsamp2;

	if(x0==NULL) x0=x;
	x1=NewDMatrix(nsamp,nv);
	xc1=NewDMatrix(nsamp,nv);
	D1=NewDMatrix(nsamp+1,nsamp+1);

	for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
	{
		if(scaled &&x0!=x)
		{
			x1[i][j]=scale[j]*x0[i][j];
			xc1[i][j]=scale[j]*(x0[i][j]-0.5);
		}
		else
		{
			x1[i][j]=x0[i][j];
			xc1[i][j]=x0[i][j]-0.5;
		}
		xc1[i][j]=ABS(xc1[i][j]);
	}
	nsamp2=nsamp*nsamp;
	D1[nsamp][nsamp]=0;
	for(i=0;i<nsamp;i++)
	{
		for(k=0,d1=1,d2=1;k<nv;k++)
		{
			d1*=1+xc1[i][k];
			d2*=(2+xc1[i][k]-xc1[i][k]*xc1[i][k])/2;
		}
		D1[i][i]=d1/nsamp2;
		D1[i][nsamp]=-2*d2/nsamp; //D1 has nsamp+1 column
		for(j=i+1;j<nsamp;j++)
		{
			for(k=0,d1=1;k<nv;k++)
			{
				d2=x1[i][k]-x1[j][k];
				d1*=(2+xc1[i][k]+xc1[j][k]-ABS(d2))/2;
			}
			D1[j][i]=D1[i][j]=2*d1/nsamp2;
		}
	}

	for(i=0,newdiscr=1;i<nactive;i++) newdiscr*=(13.0/12.0);
	for(i=0;i<nsamp+1;i++) 	for(j=i;j<nsamp+1;j++)	newdiscr+=D1[i][j];
	FreeDMatrix(D1);
	FreeDMatrix(x1);
	FreeDMatrix(xc1);
	return(newdiscr);
}

void discrcl2_full_snap()
{
	int i,j,k;
	discrfull=discr;
	for(i=np;i<nsamp;i++)
		for(j=0;j<nv;j++) xf[i-np][j]=x[i][j];
	for(k=0,i=np;i<nsamp+1;i++)
		for(j=0;j<=i;j++,k++) Df[k]=D[i][j];
}

void discrcl2_full_reset()
{
	int i,j,k;
	discr=discrfull;
	for(k=0,i=np;i<nsamp+1;i++)
		//this is necessary D[j][i]=D[i][j] for multiple pairs exchange
		for(j=0;j<=i;j++,k++) D[j][i]=D[i][j]=Df[k];

	for(i=np;i<nsamp;i++)
	{
		for(j=0;j<nv;j++)
		{
			x[i][j]=xf[i-np][j];
			if(scaled) xc[i][j]=xf[i-np][j]-0.5*scale[j];
			else xc[i][j]=xf[i-np][j]-0.5;
			xc[i][j]=ABS(xc[i][j]);
		}
	}
}

void discrcl2_global_x(double **xnew)
{
	int i,j;
	for(j=0;j<nv;j++)
	{
		for(i=0;i<np;i++)  xnew[i][j]=x[i][j];
    	for(i=np;i<nsamp;i++)  xnew[i][j]=xf[i-np][j];
	}
	if(scaled) for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) xnew[i][j]/=scale[j];
}
