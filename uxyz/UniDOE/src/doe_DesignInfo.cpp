#include "doe_DesignInfo.h"
#include "doe_Matrix.h"
#include "doe_Index.h"
#include "doe_utility.h"
#include "doe_criteria.h"

#define MINDIF 1.0e-15

//extern XINFO xinfo; // newly added from xinfo.h
//extern XLEVEL xlevel; //  newly added from xinfo.h

static unsigned  int *xidsnap,*xidtmp,**xidf;
static int nnew,nv,np,nsamp;

void create_xinfo(double **x,int nnew1,int np1,int nv1)
{
	int i,j,*idx,start,noxvl,ne;
	double *tempx;

	nnew=nnew1;
	nv=nv1;
	np=np1;
	nsamp=nnew+np;
	xinfo.xid=NewIMatrix(nnew,nv);
	xidf=(unsigned int**)NewIMatrix(nnew,nv);
	xinfo.xvl=NewIMatrix(nnew,nv);
	xinfo.nxvl=NewIVector(nv);
	xinfo.isbalanced=NewCVector(nv);


	xinfo.nle=NewIMatrix(nnew,nv);

	xidsnap=(unsigned int*)NewIVector(nnew);
	xidtmp=(unsigned int*)NewIVector(nnew);

	tempx=NewDVector(nnew);
	idx=NewIVector(nnew);

	for(j=0;j<nv;j++)
	{
		for(i=0;i<nnew;i++) tempx[i]=x[i+np][j];
		indexx1(nnew,tempx,(unsigned int*)idx);
		xinfo.isbalanced[j]=1;
		for(i=noxvl=ne=start=0;i<nnew;i++)
		{
			xinfo.xid[i][j]=(int)idx[i]+np;
			if(i==0) xinfo.xvl[i][j]=0;
			else if(tempx[idx[i]]-tempx[idx[i-1]]<MINDIF) xinfo.xvl[i][j]=noxvl;
			else
			{
				if(noxvl==0)
				{
					ne=i-start;
					if(ne==1) xinfo.isbalanced[j]=0;
				}
				else if(noxvl>0 && ne!=i-start) xinfo.isbalanced[j]=0;

				xinfo.nle[noxvl][j]=i-start;
				noxvl++;
				xinfo.xvl[i][j]=noxvl;
				start=i;
			}
		}
		if(ne!=nnew-start) xinfo.isbalanced[j]=0;

		xinfo.nle[noxvl][j]=nnew-start;
		xinfo.nxvl[j]=xinfo.xvl[nnew-1][j]+1;
	}

	FreeVector(tempx);
	FreeVector(idx);
}


void free_xinfo()
{
	FreeIMatrix(xinfo.xid);
	FreeIMatrix((int**)xidf);
	FreeIMatrix(xinfo.xvl);
	FreeVector(xinfo.nxvl);
	FreeVector(xinfo.isbalanced);
	FreeIMatrix(xinfo.nle);

	FreeVector(xidsnap);
	FreeVector(xidtmp);
}

void xinfo_snap(int ncol)
{
	int i;
	for(i=0;i<nnew;i++) xidsnap[i]=xinfo.xid[i][ncol];
}

void xinfo_reset(int ncol)
{
	int i;
	for(i=0;i<nnew;i++) xinfo.xid[i][ncol]=xidsnap[i];
}

void xinfo_full_snap()
{
	int i,j;
	for(i=0;i<nnew;i++) for(j=0;j<nv;j++) xidf[i][j]=xinfo.xid[i][j];
}

void xinfo_full_reset()
{
	int i,j;
	for(i=0;i<nnew;i++) for(j=0;j<nv;j++) xinfo.xid[i][j]=xidf[i][j];
}

void xinfo_cp(int ncol,int idx1,int idx2)
{
	int temp;
	temp=xinfo.xid[idx1][ncol]; xinfo.xid[idx1][ncol]=xinfo.xid[idx2][ncol];
	xinfo.xid[idx2][ncol]=temp;
}

void xinfo_symm_cp(int ncol,int idx1,int idx2)
{
	int temp,nt;
	temp=xinfo.xid[idx1][ncol]; xinfo.xid[idx1][ncol]=xinfo.xid[idx2][ncol];
	xinfo.xid[idx2][ncol]=temp;
	nt=nsamp+np-1;
	if(idx2!=nt-idx1)//;
	{
		idx1=nt-idx1;
		idx2=nt-idx2;
		temp=xinfo.xid[idx1][ncol]; xinfo.xid[idx1][ncol]=xinfo.xid[idx2][ncol];
		xinfo.xid[idx2][ncol]=temp;
	}
}

void xinfo_pm(int ncol,int npm,int *idx1,int *idx2)
{
	int i;
	for(i=0;i<npm;i++) xidtmp[i]=xinfo.xid[idx1[i]][ncol];
	for(i=0;i<npm;i++) xinfo.xid[idx2[i]][ncol]=xidtmp[i];
}

void xinfo_symm_pm(int ncol,int npm,int *idx1,int *idx2)
{
	int i;
	for(i=0;i<npm;i++) xidtmp[i]=xinfo.xid[idx1[i]][ncol];
	for(i=0;i<npm;i++) xinfo.xid[idx2[i]][ncol]=xidtmp[i];
	for(i=0;i<npm;i++) xidtmp[i]=xinfo.xid[nnew-idx1[i]-1][ncol];
	for(i=0;i<npm;i++) xinfo.xid[nnew-idx2[i]-1][ncol]=xidtmp[i];
}

/*
Here a level is not a fixed value, rather it is defined as a set of values, which
after sorting a column in x, are close to each other. Different levels have same
number of elements. And there are no same elements in different levels. However,
in the same levels, there may exist same elements.
For example: x=[x1 x2]. x1 (1, 2,2,3,4,4)' is two levels (1 2 2) and (3 4 4)
Another example: x1(1 1 1 2 2 2 3 3 3) three levels (1 1 1) (2 2 2) (3 3 3)

Specially, if number of levels is not provided or it violates the definition(no
elements in different levels has same value), we consider the variable as 1 level.

There are two kinds of exchange: one is between levels, the other is in a level.

 if nsamp%level[i]!=0, 1 entry per level is used,
in that case level[i] will not be used
Requirement: the elements in different levels are different. However, elements
in the same level may have same values.

       The errorness in level is detected based on whether the two elements
	   in the boundray of two  levels (after sorting) have same value.

If that's true, it will change the xinfo.levels and enties
and regard the variable as nnew levels.
*/


/*
This assume that each level has the same number of entries
The first row of levels_map is the number of entries
This function will mark duplicated elememts in levals_map with corresponding
negative value.
*/

void create_xlevel(double **x,int *levels)
{
	int i,j;

	xlevel.levels=NewIVector(nv);
	xlevel.entries=NewIVector(nv);
	xlevel.isequal=NewCVector(nv);

	for(i=0;i<nv;i++)
	{
		if(!(nnew%levels[i]) && levels[i]!=nnew)
		//if levels[i]=nnew, it is changed to levels[i]=1 for convenience
		// this will be used in pairs develop
		{
			xlevel.levels[i]=levels[i];
			xlevel.entries[i]=nnew/levels[i];
		}
		else {xlevel.levels[i]=1;xlevel.entries[i]=nnew;}
	}

	for(j=0;j<nv;j++)
	{
		for(i=xlevel.entries[j];i<nnew;i+=xlevel.entries[j])
	       if(xinfo.xvl[i][j]==xinfo.xvl[i-1][j])
		   {
				xlevel.levels[j]=1;
				xlevel.entries[j]=nnew;
				break;
		   }
	}
	for(j=0;j<nv;j++)
	{
		if(xinfo.isbalanced[j]&&xinfo.nxvl[j]==xlevel.levels[j]) xlevel.isequal[j]=1;
		else xlevel.isequal[j]=0;
	}
}


void free_xlevel()
{
	FreeVector(xlevel.levels);
	FreeVector(xlevel.entries);
	FreeVector(xlevel.isequal);
}

char checkSymm()
{
	int ne1,*tmp1,*tmp2,i,*indx1,*indx2,j;
	//there is a unfixed bug which will lead to error in updating criterion
	//induced by exchange, if the matrix is not a LHD. for example, a balanced
	//design. In this case, the symmetric exchange is disabled currently.
	for(i=0;i<nv;i++) if(xinfo.nxvl[i]!=nnew) return(0);
	ne1=nnew/2;
	tmp1=NewIVector(nnew);
	tmp2=NewIVector(nnew);
	indx1=NewIVector(ne1);
	indx2=NewIVector(ne1);
	for(i=0;i<ne1;i++)
		if(xinfo.xid[i][0]<xinfo.xid[nnew-i-1][0])
		{
			tmp1[i]=xinfo.xid[i][0]; tmp1[nnew-i-1]=xinfo.xid[nnew-i-1][0];
		}
		else
		{
			tmp1[i]=xinfo.xid[nnew-i-1][0]; tmp1[nnew-i-1]=xinfo.xid[i][0];
		}
	indexx2(ne1,tmp1,(unsigned int*)indx1);
	for(j=1;j<nv;j++)
	{
		for (i=0;i<ne1;i++)
			if(xinfo.xid[i][j]<xinfo.xid[nnew-i-1][j])
			{
				tmp2[i]=xinfo.xid[i][j]; tmp2[nnew-i-1]=xinfo.xid[nnew-i-1][j];
			}
			else
			{
				tmp2[i]=xinfo.xid[nnew-i-1][j]; tmp2[nnew-i-1]=xinfo.xid[i][j];
			}
		indexx2(ne1,tmp2,(unsigned int*)indx2);
		for (i=0;i<ne1;i++)
			if((tmp2[indx2[i]]!=tmp1[indx1[i]])||(tmp2[nnew-indx2[i]-1]!=tmp1[nnew-indx1[i]-1])) break;
		if(i<ne1) break;
	}
	FreeVector(tmp1);
	FreeVector(tmp2);
	FreeVector(indx1);
	FreeVector(indx2);
	if(j<nv) return(0);
	else return(1);
}
