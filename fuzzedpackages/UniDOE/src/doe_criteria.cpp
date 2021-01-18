#include <stdio.h>
#include <stdlib.h>
#include "doe_criteria.h"
#include "doe_maximin.h"
#include "doe_CL2.h"
#include "doe_MD2.h"
#include "doe_WD2.h"
#include "doe_utility.h"
#include "doe_Matrix.h"
#include "doe_DesignInfo.h"

static char criterion,ismax;
static int nsamp,nv;
static double goal;

void create_criteria(double **x,int nnew,int np,int nv1,CRITOPT *critopt)
{
	criterion=critopt->type;
	nsamp=nnew+np;
	nv=nv1;
	goal=critopt->goal;
	ismax=critopt->ismax;

	switch(criterion)
	{
	case 1:
		create_maximin(x,nnew,np,nv1,critopt);
		break;
	case 2:
		create_discrcl2(x,nnew,np,nv1,critopt);
		break;
	case 3:
		create_mxcl2(x,nnew,np,nv1,critopt);
		break;
	case 4:
		create_wdl2(x,nnew,np,nv1,critopt);
		break;
	default:
		create_mxcl2(x,nnew,np,nv1,critopt);
		break;
	}
}

void free_criteria(void)
{
	switch(criterion)
	{
	case 1:
		free_maximin();
		break;
	case 2:
		free_discrcl2();
		break;
	case 3:
		free_mxcl2();
		break;
	case 4:
		free_wdl2();
		break;
	default:
		free_mxcl2();
		break;
	}
}

double criteria_set(double **x)
{
	switch (criterion)
	{
	case 1:
		return(maximin_set(x));
	case 2:
		return(discrcl2_set(x));
	case 3:
		return(mxcl2_set(x));
	case 4:
		return(wdl2_set(x));
	default:
		return(mxcl2_set(x));
	}
}

double criteria_cp_set(int ncol,int ncp,int *idx1,int *idx2)
{
	switch (criterion)
	{
	case 1:
		return(maximin_cp_set(ncol,ncp,idx1,idx2));
	case 2:
		return(discrcl2_cp_set(ncol,ncp,idx1,idx2));
	case 3:
		return(mxcl2_cp_set(ncol,ncp,idx1,idx2));
	case 4:
		return(wdl2_cp_set(ncol,ncp,idx1,idx2));
	default:
		return(mxcl2_cp_set(ncol,ncp,idx1,idx2));
	}
}


double criteria_cp(int ncol,int ncp,int *idx1,int *idx2)
{
	switch (criterion)
	{
	case 1://Maximin
		return(maximin_cp(ncol,ncp,idx1,idx2));
	case 2:
		return(discrcl2_cp(ncol,ncp,idx1,idx2));
	case 3:
		return(mxcl2_cp(ncol,ncp,idx1,idx2));
	case 4:
		return(wdl2_cp(ncol,ncp,idx1,idx2));
	default:
		return(mxcl2_cp(ncol,ncp,idx1,idx2));
	}
}

double criteria_cp1(int ncol,int idx1,int idx2)
{
	switch (criterion)
	{
	case 1://Maximin
		return(maximin_cp1(ncol,idx1,idx2));
	case 2:
		return(discrcl2_cp1(ncol,idx1,idx2));
	case 3:
		return(mxcl2_cp1(ncol,idx1,idx2));
	case 4:
		return(wdl2_cp1(ncol,idx1,idx2));
	default:
		return(mxcl2_cp1(ncol,idx1,idx2));
	}
}

double criteria()
{
	switch (criterion)
	{
	case 1: //Maximin
		return(maximin());
	case 2:
		return(discrcl2());
	case 3:
		return(mxcl2());
	case 4:
		return(wdl2());
	default:
		return(mxcl2());
	}
}

double **criteria_x(double **x)
{
	switch (criterion)
	{
	case 1://Maximin
		return(maximin_x(x));
    case 2:
		return(discrcl2_x(x));
    case 3:
		return(mxcl2_x(x));
	case 4:
		return(wdl2_x(x));
	default:
		return(mxcl2_x(x));
	}
}

void criteria_snap(int ncol)
{
	switch (criterion)
	{
	case 1://Maximin
		maximin_snap(ncol); return;
    case 2:
		discrcl2_snap(ncol); return;
    case 3:
		mxcl2_snap(ncol); return;
	case 4:
		wdl2_snap(ncol); return;
	default:
		mxcl2_snap(ncol); return;
	}

}

void criteria_reset(int ncol)
{
	switch (criterion)
	{
	case 1://Maximin
		maximin_reset(ncol); return;
	case 2:
		discrcl2_reset(ncol); return;
	case 3:
		mxcl2_reset(ncol); return;
	case 4:
		wdl2_reset(ncol); return;
	default:
		mxcl2_reset(ncol); return;
	}
}

double criteria_pm(int ncol, int npm,int *idx,int *idxp)
{
	switch (criterion)
	{
	case 1://Maximin
		return(maximin_pm(ncol,npm,idx,idxp));
	case 2:
		return(discrcl2_pm(ncol,npm,idx,idxp));
	case 3:
		return(mxcl2_pm(ncol,npm,idx,idxp));
	case 4:
		return(wdl2_pm(ncol,npm,idx,idxp));
	default:
		return(mxcl2_pm(ncol,npm,idx,idxp));
	}
}

double criteria_pm_set(int ncol, int npm,int *idx,int *idxp)
{
	switch (criterion)
	{
	case 1://Maximin
		return(maximin_pm_set(ncol,npm,idx,idxp));
	case 2:
		return(discrcl2_pm_set(ncol,npm,idx,idxp));
	case 3:
		return(mxcl2_pm_set(ncol,npm,idx,idxp));
	case 4:
		return(wdl2_pm_set(ncol,npm,idx,idxp));
	default:
		return(mxcl2_pm_set(ncol,npm,idx,idxp));
	}
}

double criteria_eval(double **x)
{
	switch (criterion)
	{
	case 1:
		return(maximin_eval(x));
	case 2:
		return(discrcl2_eval(x));
	case 3:
		return(mxcl2_eval(x));
	case 4:
		return(wdl2_eval(x));
	default:
		return(mxcl2_eval(x));
	}
}

double criteria_min()
{
	return(goal);
}

void criteria_full_snap()
{
	switch (criterion)
	{
	case 1://Maximin
		maximin_full_snap(); return;
    case 2:
		discrcl2_full_snap(); return;
    case 3:
		mxcl2_full_snap(); return;
	case 4:
		wdl2_full_snap(); return;
	default:
		mxcl2_full_snap(); return;
	}
}

void criteria_full_reset()
{
	switch (criterion)
	{
	case 1://Maximin
		maximin_full_reset(); return;
	case 2:
		discrcl2_full_reset(); return;
	case 3:
		mxcl2_full_reset(); return;
	case 4:
		wdl2_full_reset(); return;
	default:
		mxcl2_full_reset(); return;
	}
}

void criteria_global_x(double **x)
{
	switch (criterion)
	{
	case 1://Maximin
		maximin_global_x(x);
		break;
    case 2:
		discrcl2_global_x(x);
		break;
    case 3:
		mxcl2_global_x(x);
		break;
	case 4:
		wdl2_global_x(x);
		break;
	default:
		mxcl2_global_x(x);
	}
}
