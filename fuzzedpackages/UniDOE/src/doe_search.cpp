//update: 20/11/2017,changes in soat_search
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "doe_search.h"
#include "doe_Matrix.h"
#include "doe_criteria.h"
#include "doe_Index.h"
#include "doe_random.h"
#include "doe_utility.h"
#include "doe_DesignInfo.h"
#include "doe_CL2.h"
using namespace std;
#define DEBUG 1
typedef struct
{
	char type;
	int size;
	int ncp; //number of exchanges in the same column,if type=1,ncp=1
	int *pi,*pj; //1*m
}PAIRSTABLE;

typedef struct
{
	char type;
	int size;
	int ne; //number of entries in each level,if type=1,ne=1
	int *pi,*pj;
}PERMUTETABLE;

//parameters changed by initialization
//Here, maxcol is M in paper, i.e. number of inner loops
static int maxiter,maxcol; //maximum exchanging columns in each iteration
static double maxtime,minchange,hits_ratio;
static int maxpairs,set_count=0;
static int nsamp,nnew,np,nv;
static char /*issymm,*/israndcol,israndpairs,isperm;

//threshold
static double tol,th0,factor,thmin,thmax,global_obj, obj,goal;
static double *ncol_prob;
static int *active_col,nactive_col,*nlpairs,*nepairs,niter;
static long ntotal;

static PERMUTETABLE pmtab;
static PAIRSTABLE ptab;

//local functions
void create_pairs_table(int npairs);
void create_permute_table();

void free_pairs_table();
void free_permute_table();

void save_before_change(int ncol);
void restore_after_change(int ncol);
void save_global();
void restore_global();
//void print_x(double **x);

int get_ncol(char isrnd);

int get_level_pairs(int ncol,int npairs);
int get_element_pairs(int ncol,int npairs);
int get_element_pairs_balanced(int ncol,int npairs);

double get_level_exchange(int ncol,int npairs);
double get_element_exchange(int ncol, int npairs,char flag);
double get_element_symm_exchange(int ncol, int npairs,char flag);
double find_exchange(int ncol,double cur_obj);

int get_level_set(int ncol,int npm);
int get_element_set(int ncol,int npm,int level);

double get_level_permute(int ncol,int npm, int nsets);
double get_element_permute(int ncol,int npm,int nsets);
double get_element_symm_permute(int ncol,int npm,int nsets);
double find_permute(int ncol,int npml,int npme,int nsets);
double full_permute();

//Update 2017.oct.13: automatically updating J in paper based on exploration
//and improving process.
void updateNepairs(double alpha);
std::vector<double> soat_search(double **x);


XINFO xinfo;
XLEVEL xlevel;

void create_search(double **x,int nnew1,int np1,int nv1,SEARCHOPT *options)
{
	int i,m,nl,ne,nexvl,totall,totale;

	double tweight,prob;
	////////////
	//TEST_INTEGER_VAR = 233;
	//printf("TEST_INTEGER_VAR = %d \n",TEST_INTEGER_VAR);
  ///////////
	srandomt();
	nnew=nnew1;
	np=np1;
	nv=nv1;
	nsamp=np+nnew;;
	maxtime=options->maxtime*CLOCKS_PER_SEC;
	maxpairs=options->maxpairs;
	//isprint=options->isprint;
  th0=options->th0;
	tol=options->tol;
	factor=options->factor;
	//issymm=options->issymm;
	israndcol=options->israndcol;
	israndpairs=options->israndpairs;
	//nochangeiter=options->nochangeiter;
	isperm=options->isperm;
	hits_ratio = options->hits_ratio;
	#if DEBUG
		//printf("In create_search of search.cpp before checkSymm, issymm = %d \n",issymm);
	#endif
	//if(issymm) issymm=checkSymm();
 	#if DEBUG
		//printf("In create_search of search.cpp after checkSymm, issymm = %d \n",issymm);
	#endif
	create_permute_table();
	create_xlevel(x,options->levels);
	//for(i=0;i<nv;i++) printf("\n In create_search function, options->levels[i] = %d",options->levels[i]);
	ncol_prob=NewDVector(nv);
	active_col=NewIVector(nv);

	for (i=0,tweight=EPS;i<nv;i++)
		if(options->colweight[i]>EPS) tweight+=options->colweight[i];
	for (i=0,m=0,prob=0;i<nv;i++)
	{
		if(options->colweight[i]>EPS)
		{
			ncol_prob[m]=options->colweight[i]/tweight+prob;
			active_col[m]=i;
			prob=ncol_prob[m];
			m++;
		}
	}
	nactive_col=m;
	ncol_prob[m-1]=1+EPS;
	nlpairs=NewIVector(nv);
	nepairs=NewIVector(nv);

	for(i=0,maxcol=0;i<nv;i++)
	{
		/*
			xlevel.entries[i] stores number of elements in column i
			xinfo.nxvl[i] stores number of levels in column i
		*/
		ne=xlevel.entries[i];
		nl=xlevel.levels[i];
		if(xinfo.isbalanced[i]) nexvl=nnew/xinfo.nxvl[i];
		else nexvl=1;
		totall=nl*(nl-1)/2;  totale=nl*ne*(ne-nexvl)/2;
		if(israndpairs)
		{
			if(nl>1)
			{
				//nlpairs is a vector that stores number of level permutations occuring in each column.
				//maximum value of it is approximately half of full permutation.
				nlpairs[i]=MIN((totall+1)/2,maxpairs); nepairs[i]=MIN((totale+1)/2,maxpairs);
			}
			else
			{
				nlpairs[i]=0;
				if(totale<=20) nepairs[i]=1;
				else nepairs[i]=(totale+1)/5;
				nepairs[i]=MIN(nepairs[i],maxpairs);
			}
		}
		else {nlpairs[i]=totall; nepairs[i]=totale;}
		// Below are procedures to decide number of inner loops.
		// If nlpairs and nepairs are not zero, maxcol is at least 4 * nb.of factors(i.e. columns)
		if(nlpairs[i]==0) maxcol=MAX(maxcol,(int)(2.0*nactive_col*totale/nepairs[i]));
		else if(nepairs[i]==0) maxcol=MAX(maxcol,(int)(2.0*nactive_col*totall/nlpairs[i]));
		else maxcol=MAX(maxcol,(int)(2* /*nactive_col*/nv*MAX(1.0*totall/nlpairs[i],1.0*totale/nepairs[i])));
	}
	maxcol=MAX(10,maxcol);
	maxcol=MIN(maxcol,100);
	maxcol=iCheckValue(1,100000000,maxcol,options->maxcol);

	maxiter=MAX(nactive_col*nnew*2,100);
	maxiter=iCheckValue(1,100000000,maxiter,options->maxiter);
	for(i=0,maxpairs=0;i<nv;i++)
	{
		maxpairs=MAX(maxpairs,nepairs[i]);
		maxpairs=MAX(maxpairs,nlpairs[i]);
	}

	create_pairs_table(maxpairs);
	#if DEBUG
		// printf("In search.cpp create_search: nepairs[0] = %d; nlpairs[0] = %d; maxpairs = %d; maxiter = %d\n",nepairs[0],nlpairs[0],maxpairs,maxiter);
		// printf("In search.cpp create_search: ptab.size = %d; pmtab.size = %d\n",ptab.size,pmtab.size);
		// printf("In search.cpp create_search: ptab.ncp = %d; pmtab.ne = %d\n",ptab.ncp,pmtab.ne);
		// for(i=0;i<nv;i++)
		// 	printf("In search.cpp create_search: options->levels[%d] = %d; xinfo.nxvl[i] = %d; xlevel.entries[i] = %d\n",i,options->levels[i],xinfo.nxvl[i],xlevel.entries[i]);
	#endif
	goal=criteria_min();
}


void get_new_options(SEARCHOPT *options)
{
	int i;
	options->factor=factor;
	options->th0=th0;
	options->tol=tol;
	//options->issymm=issymm;
  for(i=0;i<nv;i++) options->levels[i]=xlevel.levels[i];
	options->maxcol=maxcol;
	options->maxiter=maxiter;
	options->maxpairs=maxpairs;
}

void free_search()
{
	free_pairs_table();
	free_permute_table();
	free_xlevel();
	FreeVector(ncol_prob);
	FreeVector(active_col);
	FreeVector(nlpairs);
	FreeVector(nepairs);
}


void save_before_change(int ncol)
{
  criteria_snap(ncol);
	xinfo_snap(ncol);
}


void restore_after_change(int ncol)
{
	criteria_reset(ncol);
	xinfo_reset(ncol);
}


/* type 0: exchange between levels
   type 1: exchange between entries in the same level
*/
void create_pairs_table(int npairs)
{
	ptab.pi=NewIVector(npairs);
	ptab.pj=NewIVector(npairs);
}


void free_pairs_table()
{
	FreeVector(ptab.pi);
	FreeVector(ptab.pj);
}


void create_permute_table()
{
	pmtab.pi=NewIVector(nnew);
	pmtab.pj=NewIVector(nnew);
}


void free_permute_table()
{
	FreeVector(pmtab.pi);
	FreeVector(pmtab.pj);
}


void get_pairs_id(int tk,int ne,int *id1,int *id2)
{
	int p1,p2,i,j;
	p1=0; p2=ne-1;
	while(p2-p1>1)
	{
		i=(p1+p2)/2;
		if(tk>(2*ne-i-1)*i/2) p1=i;
		else p2=i;
	}
	if(tk>=(2*ne-p2-1)*p2/2) i=p2;
	else i=p1;
	j=tk-(2*ne-i-1)*i/2+i+1;
	*id1=i;
	*id2=j;
}


/*
ptab.ncp represents the number of elements in a level
ptab.type = 0 means doing level permutation between different levels, used in UDC
ptab.size represents number of level permutations happening in a column
xlevel.entries[*] means number of elements in a level
xlevel.levels[*] : number of levels in a column
*/
int get_level_pairs(int ncol,int npairs)
{
	int i,j,nl,/*ncp,*/m,total,tk,*indx,*spairs,maxsch,k;
  //unsigned int *indx;
	//char flag=0;
	ptab.type=0;
	/*ncp=*/ptab.ncp=xlevel.entries[ncol];
	nl=xlevel.levels[ncol];
  total=nl*(nl-1)/2;
	if(total==0) {ptab.size=0; return(0);}
	npairs=MIN(total,npairs);

  if(total>npairs)
	{
		maxsch=MAX(50,2*npairs); //maximum search
		spairs=NewIVector(maxsch);
		indx=(int*)NewLVector(maxsch);
		for(i=0;i<maxsch;i++) spairs[i]=(int) (Random()*total);
		indexx2(maxsch,spairs,(unsigned int*)indx);
		for(i=maxsch-1;i>0;i--) if(spairs[indx[i]]==spairs[indx[i-1]]) spairs[indx[i]]=-1;

		m=0; k=0;
		while(m<npairs && k<maxsch)
		{
			if(spairs[k]!=-1) {tk=spairs[k]; k++;}
			else {k++; continue;}
			get_pairs_id(tk,nl,&ptab.pi[m],&ptab.pj[m]);
			m++;
		}

		FreeVector(spairs);
		FreeVector(indx);
	}
	else for(i=0,m=0;i<nl;i++)
		for(j=i+1;j<nl;j++) {ptab.pi[m]= i; ptab.pj[m]= j; m++; }

	ptab.size=m;
	return(ptab.size);
}


//generate a set of exchange pairs. The two elements in each
//pair will be in the same group level. Different pairs may
//in different group level
int get_element_pairs(int ncol,int npairs)
{
	int i,j,t,nl,maxsch,ne1;
	int ne,total,total1,tk,k,m,*spairs,*indx;
	char isodd=0;

	if(xinfo.isbalanced[ncol]) return(get_element_pairs_balanced(ncol,npairs));

	ptab.type=1; ptab.ncp=1;
	ne=xlevel.entries[ncol];
	nl=xlevel.levels[ncol];
////////////
	//if(issymm) {ne1=nnew/2; total1=ne1*(2*nnew-ne1-1)/2; if(nnew%2) isodd=1;}
	//else {ne1=ne; total1=ne*(ne-1)/2;}
	ne1=ne;
	total1=ne*(ne-1)/2;
////////////
	total=nl*total1;
	npairs=MIN(npairs,total);
	if(npairs==0||xinfo.nxvl[ncol]==1) {ptab.size=0; return(0);}

	if(total>npairs)
	{
		maxsch=(int)(2.0*nnew/xinfo.nxvl[ncol]*npairs); //maximum search
		spairs=NewIVector(maxsch);
		indx=(int*)NewLVector(maxsch);
		for(i=0;i<maxsch;i++)
		{
			if(nl>1) t=(int)(nl*Random());
			else t=0;
			tk=(int) (Random()*total1);
			spairs[i]=t*total1+tk;
		}
		indexx2(maxsch,spairs,(unsigned int*)indx);
		for(i=maxsch-1;i>0;i--) if(spairs[indx[i]]==spairs[indx[i-1]]) spairs[indx[i]]=-1;

		m=0; k=0;
		while(m<npairs && k<maxsch)
		{
			if(spairs[k]!=-1) {t=spairs[k]/total1; tk=spairs[k]%total1; k++;}
			else {k++; continue;}

			get_pairs_id(tk,ne,&i,&j);

			if(isodd && j==ne1) continue;

			t=t*ne;
			if(xinfo.xvl[t+i][ncol]!=xinfo.xvl[t+j][ncol])
			{
				ptab.pi[m]= t+i; ptab.pj[m]= t+j; m++;
			}
		}
		FreeVector(indx);
		FreeVector(spairs);
	}
	else
	{
		for(k=0,m=0;k<nl;k++)
		{
			t=k*ne;
			for(i=0;i<ne1;i++) for(j=i+1;j<ne;j++)
			{
				if(isodd && j==ne1) continue;
				if(xinfo.xvl[t+i][ncol]!=xinfo.xvl[t+j][ncol])
				{
					ptab.pi[m]= t+i; ptab.pj[m]= t+j; m++;
				}
			}
		}
	}
	ptab.size=m;

	return(ptab.size);
}


int get_element_pairs_balanced(int ncol,int npairs)
{
	int i,j,t,nl,maxsch,tk1,ttk,i1,j1,ipdx=0;
	int ne,total,total1,tk,k,m,*spairs,*indx,nxvl,nexvl,nexvl2;

	ptab.type=1; ptab.ncp=1;
	ne=xlevel.entries[ncol];
	nl=xlevel.levels[ncol];

	nexvl=nnew/xinfo.nxvl[ncol];
	nexvl2=nexvl*nexvl;
	nxvl=ne/nexvl;

	total1=nxvl*(nxvl-1)/2;
	total=nl*nexvl2*total1;
	npairs=MIN(npairs,total);
	if(npairs==0) {ptab.size=0; return(0);}

	if(total>npairs)
	{
		maxsch=4; //maximum search
		spairs=NewIVector(2*npairs);
		indx=(int*)NewLVector(2*npairs);
		m=k=0;
		while(m<npairs && k<maxsch)
		{
			for(i=0;i<npairs;i++)
			{
				if(nl>1) t=(int)(nl*Random());
				else t=0;
				if(total1>1) tk=(int) (Random()*total1);
				else tk=0;
				tk1=(int) (Random()*nexvl2);
				spairs[i]=nexvl2*(t*total1+tk)+tk1;
			}
			indexx2(npairs+m,spairs,(unsigned int*)indx);
			for(i=0;i<npairs+m-1;i++) if(spairs[indx[i]]==spairs[indx[i+1]]) spairs[indx[i]]=-1;
			for(i=npairs+m-1; i>=0;i--)
			{
				if(indx[i]<npairs && spairs[indx[i]]!=-1) ipdx=indx[i];
				else if(indx[i]>=npairs && spairs[indx[i]]==-1)
				{
					spairs[indx[i]]=spairs[ipdx]; spairs[ipdx]=-1;
				}
			}

			for(i=0;i<npairs;i++)
				if(spairs[i]!=-1)
				{
					spairs[npairs+m]=spairs[i];
					m++;
					if(m==npairs) break;
				}
			k++;
		}

		for(k=0;k<m;k++)
		{
			ttk=spairs[npairs+k]/nexvl2; tk1=spairs[npairs+k]%nexvl2;
			t=ttk/total1;  tk=ttk%total1;

			get_pairs_id(tk,nxvl,&i,&j);
			i1=tk1/nexvl;  j1=tk1%nexvl;

			t=t*ne;
			ptab.pi[k]= t+i*nexvl+i1; ptab.pj[k]= t+j*nexvl+j1;
		}
		FreeVector(indx);
		FreeVector(spairs);
	}
	else
	{
		for(k=0,m=0;k<nl;k++)
		{
			t=k*ne;
			for(i=0;i<nxvl;i++) for(j=i+1;j<nxvl;j++)
			{
				for(i1=0;i1<nexvl;i1++) for(j1=0;j1<nexvl;j1++)
				{
					ptab.pi[m]= t+i*nexvl+i1; ptab.pj[m]= t+j*nexvl+j1; m++;
				}
			}
		}
	}
	ptab.size=m;
	return(ptab.size);
}


//R and x will be changed
/*
	get_level_exchange fucntion:
	ncol: the column that level permutation happens.\
	npairs: number of pair of levels that will exchange to find the best result during level permutation.
				e.g. If pairs of level = [(1,2),(2,3)], then npairs = 2.
*/
double get_level_exchange(int ncol,int npairs)
{

	int i,ncp,ti,tj,pair_id,*idx1,*idx2,j;
	unsigned int *rank;
	double *critobj;

	/*
		get_level_pairs function is used to generate level pairs, stored in ptab.
		ptab.pi and ptab.pj are two vectors to store level pairs to be exchanged.
		ptab.ncp is number of entries in a level.
		ptab.size is the actual number of pairs of levels to be exchanged.
	*/
	if(!(npairs=get_level_pairs(ncol,npairs))) return(criteria()) ;
	ncp=ptab.ncp;

	rank=(unsigned int*)NewLVector(npairs); critobj=NewDVector(npairs);
	idx1=NewIVector(ncp); 	idx2=NewIVector(ncp);

	//printf("\n In get_level_exchange function, ptab.size = %d ",ptab.size);

	for(i=0;i<ptab.size;i++)
	{
		ti=ncp*ptab.pi[i]; tj=ncp*ptab.pj[i];
		for(j=0;j<ncp;j++)
		{
			idx1[j]=xinfo.xid[ti+j][ncol]; 	idx2[j]=xinfo.xid[tj+j][ncol];
		}
		/*
			criteria_cp fucntion is used to change two elements in a column
			and calculate the discrepancy of design after the change.
		*/
		critobj[i]=criteria_cp(ncol,ncp,idx1,idx2);
	}

  /*
		indexx1 is used to rank the id according to corresponding criterion value.
	  In an assending rank
	*/
	indexx1(npairs,critobj,rank);
	/* select element for exchange */
	//pair_id is the id of level pair exchanged in ptab.pi and ptab.pj, which achieves lowest discrepancy.
	pair_id = rank[0];

	ncp=ptab.ncp;

	ti=ncp*ptab.pi[pair_id]; tj=ncp*ptab.pj[pair_id];
	for(i=0;i<ncp;i++)
	{
		idx1[i]=xinfo.xid[ti+i][ncol]; 	idx2[i]=xinfo.xid[tj+i][ncol];
		//change xinfo.xid correspondint to x
		xinfo_cp(ncol,ti+i,tj+i);
	}
	//criteria_cp_set will change x and R correspondingly
	criteria_cp_set(ncol,ptab.ncp,idx1,idx2);

	FreeVector(idx1);
	FreeVector(idx2);

	FreeVector(rank);
	FreeVector(critobj);
	return(criteria());
}


double get_element_exchange(int ncol,int npairs,char flag)
{
  //printf("\n In get_element_exchange function, ncol = %d; npairs = %d; flag=%d ",ncol,npairs,flag);

	int pair_id,idx1,idx2,ti,tj,i;
	unsigned int *rank;
	double *critobj;

	//if(issymm) return(get_element_symm_exchange(ncol,npairs,flag));
	if(!(npairs=get_element_pairs(ncol,npairs))) return(criteria());

	critobj = NewDVector(npairs);
	rank=(unsigned int *)NewLVector(npairs);

	for(i=0;i<ptab.size;i++)
	{
		idx1=xinfo.xid[ptab.pi[i]][ncol]; idx2=xinfo.xid[ptab.pj[i]][ncol];
		critobj[i]=criteria_cp1(ncol,idx1,idx2);
		//printf("\n In get_element_exchange function, critobj[%d]=%f ",i,critobj[i]);
	}
	   // acsending rank
	indexx1(npairs,critobj,rank);
	   //select element for exchange

	pair_id = rank[0];
	if(!flag||critobj[pair_id]<criteria())
	{
		ti=ptab.pi[pair_id];
		tj=ptab.pj[pair_id];


		idx1=xinfo.xid[ti][ncol]; idx2=xinfo.xid[tj][ncol];

		//criteria_cp_set will change x and R correspondingly
		criteria_cp_set(ncol,ptab.ncp,&idx1,&idx2);
		//change xid corresponding the change of x
		xinfo_cp(ncol,ti,tj);
	}

	FreeVector(rank);
	FreeVector(critobj);
	return(criteria());
}


double get_element_symm_exchange(int ncol,int npairs,char flag)
{
	int pair_id,idx1[2],idx2[2],ti,tj,i,ncp;
	unsigned int *rank;
	double *critobj;

	if(!(npairs=get_element_pairs(ncol,npairs))) return(criteria());

	critobj = NewDVector(npairs);
	rank=(unsigned int *)NewLVector(npairs);

	for(i=0;i<ptab.size;i++)
	{
		idx1[0]=xinfo.xid[ptab.pi[i]][ncol]; idx2[0]=xinfo.xid[ptab.pj[i]][ncol];
		ncp=1;
		if(ptab.pi[i]!=nnew-ptab.pj[i]-1)
		{
			idx1[1]=xinfo.xid[nnew-ptab.pi[i]-1][ncol]; idx2[1]=xinfo.xid[nnew-ptab.pj[i]-1][ncol];
			ncp=2;
		}
		critobj[i]=criteria_cp(ncol,ncp,idx1,idx2);
	}
	   /* assending rank*/
	indexx1(npairs,critobj,rank);
	    /* select element for exchange */

	pair_id =  rank[0];

	if(!flag||critobj[pair_id]<criteria())
	{
		ti=ptab.pi[pair_id];
		tj=ptab.pj[pair_id];

		idx1[0]=xinfo.xid[ti][ncol]; idx2[0]=xinfo.xid[tj][ncol];
		ncp=1;
		if(ti!=nnew-tj-1)
		{
			idx1[1]=xinfo.xid[nnew-ti-1][ncol]; idx2[1]=xinfo.xid[nnew-tj-1][ncol];
			ncp=2;
		}
		//criteria_cp_set will change x and R correspondingly
		criteria_cp_set(ncol,ncp,idx1,idx2);
		//change xid corresponding the change of x
		xinfo_symm_cp(ncol,ti,tj);
	}

	FreeVector(rank);
	FreeVector(critobj);
	return(criteria());
}


//return newobj and corresponding idx: this is not necessarily used in search
//x (ncol) and R will be changed, the caller should save a temp
double find_exchange(int ncol,double cur_obj)
{
	double ex_obj;
	set_count++;
	/*
		Update 2017/Oct/10 :
		When running element exchange in UniDOE, xlevel.isequal[ncol]=0 and xlevel.levels[ncol]=1
		When running level permutation for a general design instead of LHS, xlevel.isequal[ncol]=1
	*/
	if(xlevel.isequal[ncol]) {ex_obj=get_level_exchange(ncol,nlpairs[ncol]);}
	else
	{
		if(xlevel.levels[ncol]>1 && Random()>0.5)
		{
			ex_obj=get_level_exchange(ncol,nlpairs[ncol]);
			ex_obj=get_element_exchange(ncol,nepairs[ncol],(char)(ex_obj<cur_obj));
		}
		else ex_obj=get_element_exchange(ncol,nepairs[ncol],0);
	}
	if(set_count>40) {ex_obj=criteria_set(NULL); set_count=0;}
	return(ex_obj);
}


//for permutation, npm the number of levels/elements
int get_level_set(int ncol,int npm)
{
	int *tmp1,*tmp2;
	int i,nl,/*ne,*/k;
	nl=xlevel.levels[ncol];
	tmp1=NewIVector(nl);
	pmtab.type=0;
	/*ne=*/pmtab.ne=xlevel.entries[ncol];
    pmtab.size=MIN(npm,nl);
	//permutation should be no less than 3
	if(pmtab.size<=2)  {pmtab.size=0; return(0);}
	for(i=0;i<nl;i++) tmp1[i]=i;
	permute(tmp1,nl);
	tmp2=NewIVector(pmtab.size);
	for(i=0;i<pmtab.size;i++) tmp2[i]=tmp1[i];
	permute(tmp2,pmtab.size);
	for(i=0,k=0;i<pmtab.size;i++)
		if(tmp1[i]!=tmp2[i])
		{
			pmtab.pi[k]=tmp1[i]; pmtab.pj[k]=tmp2[i];
			k++;
		}
	pmtab.size=k;
	FreeVector(tmp1);
	FreeVector(tmp2);

	return(pmtab.size);
}


//get a permutation for a given level group
//if level<0, select a random level group
int get_element_set(int ncol,int npm,int level)
{
	int *tmp1,*tmp;
	int i,j,k,t,size,ne,nl,np_lev,ix=0;
	pmtab.type=1;
	pmtab.ne=1;
//////////////
	//if(issymm) ne=nnew/2;
	//else ne=xlevel.entries[ncol];
	ne=xlevel.entries[ncol];
//////////////
	size=MIN(ne,npm);
	//permutation should be no less than 3
	if(size<=2)  {pmtab.size=0; return(0);}

	nl=xlevel.levels[ncol];

	tmp1=NewIVector(ne);

	//one group level
	if(nl==1) t=0;
	//work on given group level
	else if(level>=0) t=MIN(level,nl-1)*ne;
	//work on a random group level
	else t=(int)(nl*Random())*ne;
	for(i=0;i<ne;i++) tmp1[i]=t+i;



	//with duplicates?
	if(nnew/xinfo.nxvl[ncol]>1)
	{
		tmp=NewIVector(xinfo.nxvl[ncol]);
    //search a set in which not all entries are the same (in the case size<ne)
		for(k=0;k<20;k++)
		{
		   permute(tmp1,ne);
		   for(i=0;i<size;i++) tmp[xinfo.xvl[tmp1[i]][ncol]]++;
		   np_lev=0;
		   for(i=0;i<xinfo.nxvl[ncol];i++)
			   if(tmp[i]) {np_lev++; ix=i;}
		   if(np_lev>1) break;
		   else tmp[ix]=0; //avoid accumulation
		}

		//this is to reduce the chance of duplicate permutations (due to duplicate enties)
		for(i=0;i<size;i++) pmtab.pi[i]=tmp1[i];
		for(i=0,j=0;i<=ix;i++)
		{
			//each iteration, pj take the first tmp[i] elements in tmp1
			permute(&tmp1[j],size-j);
			for(k=j;k<j+tmp[i];k++) pmtab.pj[k]=tmp1[k];
			j=k;
		}
		FreeVector(tmp);
	}
	else
	{
		// This is necessary when only part of tmp1 will be used (size<ne)
		permute(tmp1,ne);
		for(i=0;i<size;i++)	pmtab.pi[i]=pmtab.pj[i]=tmp1[i];
		permute(pmtab.pj,size);
	}

	for(i=0;i<size;i++)
		if(xinfo.xvl[pmtab.pi[i]][ncol]!=xinfo.xvl[pmtab.pj[i]][ncol]) break;
	if(i==size) pmtab.size=0;
	else pmtab.size=size;

	FreeVector(tmp1);

	return(pmtab.size);
}


double get_level_permute(int ncol,int npm,int nsets)
{
	int i,j,total,ne,k,m;
	int *idx1,*idx2,*pibst,*pjbst,ti,tj,nl;
	double critobj,critbst;

	ne=xlevel.entries[ncol];
	nl=xlevel.levels[ncol];
	//nl >=3
	nsets=MIN(nsets,nl*(nl-1)*(nl-2));
	if(nsets==0) return(criteria());

	total=ne*npm;
	idx1=NewIVector(total); 	idx2=NewIVector(total);
	pibst=NewIVector(npm);  	pjbst=NewIVector(npm);

	critbst=MAXDOUBLE;
	for(i=0;i<nsets;i++)
	{
		if(get_level_set(ncol,npm)==0) continue;
		for(j=0,k=0;j<pmtab.size;j++)
		{
			ti=ne*pmtab.pi[j]; tj=ne*pmtab.pj[j];
			for(m=0;m<pmtab.ne;m++,k++)
			{
				idx1[k]=xinfo.xid[ti+m][ncol];
				idx2[k]=xinfo.xid[tj+m][ncol];
			}
		}
		critobj=criteria_pm(ncol,k,idx1,idx2);

		if(critobj<critbst)
		{
			critbst=critobj;
			for(j=0;j<pmtab.size;j++) {pibst[j]=pmtab.pi[j]; pjbst[j]=pmtab.pj[j];}
		}
	}
	if(critbst==MAXDOUBLE) return(criteria());

	total=pmtab.size*ne;
	for(j=0,k=0;j<pmtab.size;j++)
	{
		ti=ne*pibst[j]; tj=ne*pjbst[j];
		for(m=0;m<pmtab.ne;m++,k++)
		{
			idx1[k]=xinfo.xid[ti+m][ncol];
			idx2[k]=xinfo.xid[tj+m][ncol];
		}
	}

	criteria_pm_set(ncol,total,idx1,idx2);
	for(j=0,k=0;j<pmtab.size;j++)
	{
		ti=ne*pibst[j]; tj=ne*pjbst[j];
		for(m=0;m<pmtab.ne;m++,k++)
		{
			idx1[k]=ti+m; 	idx2[k]=tj+m;
		}
	}
	xinfo_pm(ncol,total,idx1,idx2);

	FreeVector(idx1);
	FreeVector(idx2);
	FreeVector(pibst);
	FreeVector(pjbst);
	return(critbst);
}


double get_element_permute(int ncol,int npm,int nsets)
{
	int i,j,ne,nl;
	int *idx1,*idx2,*pibst,*pjbst;
	double critobj,critbst;
	//if(issymm) return(get_element_symm_permute(ncol,npm,nsets));
	ne=xlevel.entries[ncol];
	nl=xlevel.levels[ncol];

	//if ne=2,no permute; ne>=3
	nsets=MIN(nsets,nl*ne*(ne-1)*(ne-2));

	if(nsets==0) return(criteria());

	idx1=NewIVector(npm);
	idx2=NewIVector(npm);
	pibst=NewIVector(npm);
	pjbst=NewIVector(npm);

	critbst=MAXDOUBLE;
	for(i=0;i<nsets;i++)
	{
		if(get_element_set(ncol,npm,-1)==0) continue;
		for(j=0;j<pmtab.size;j++)
		{
			idx1[j]=xinfo.xid[pmtab.pi[j]][ncol];
			idx2[j]=xinfo.xid[pmtab.pj[j]][ncol];
		}
		critobj=criteria_pm(ncol,pmtab.size,idx1,idx2);
		if(critobj<critbst)
		{
			critbst=critobj;
			for(j=0;j<pmtab.size;j++) {pibst[j]=pmtab.pi[j]; pjbst[j]=pmtab.pj[j];}
		}
	}

	if(critbst==MAXDOUBLE) return(criteria());
	if(xlevel.levels[ncol]==1 ||criteria()>critbst)
	{
		for(j=0;j<pmtab.size;j++)
		{
			idx1[j]=xinfo.xid[pibst[j]][ncol];
			idx2[j]=xinfo.xid[pjbst[j]][ncol];
		}
		criteria_pm_set(ncol,pmtab.size,idx1,idx2);
		xinfo_pm(ncol,pmtab.size,pibst,pjbst);
	}
	FreeVector(idx1);
	FreeVector(idx2);
	FreeVector(pibst);
	FreeVector(pjbst);
	return(critbst);
}


double get_element_symm_permute(int ncol,int npm,int nsets)
{
	int i,j;
	int *idx1,*idx2,*pibst,*pjbst;
	double critobj,critbst;

	nsets=MIN(nsets,nnew*(nnew-1)*(nnew-2));
	if(nsets==0) return(criteria());

	idx1=NewIVector(2*npm);
	idx2=NewIVector(2*npm);
	pibst=NewIVector(npm);
	pjbst=NewIVector(npm);

	critbst=MAXDOUBLE;
	for(i=0;i<nsets;i++)
	{
		if(get_element_set(ncol,npm,-1)==0) continue;
		for(j=0;j<pmtab.size;j++)
		{
			idx1[j]=xinfo.xid[pmtab.pi[j]][ncol];
			idx2[j]=xinfo.xid[pmtab.pj[j]][ncol];
			idx1[pmtab.size+j]=xinfo.xid[nnew-pmtab.pi[j]-1][ncol];
			idx2[pmtab.size+j]=xinfo.xid[nnew-pmtab.pj[j]-1][ncol];
		}
		critobj=criteria_pm(ncol,pmtab.size*2,idx1,idx2);
		if(critobj<critbst)
		{
			critbst=critobj;
			for(j=0;j<pmtab.size;j++) {pibst[j]=pmtab.pi[j]; pjbst[j]=pmtab.pj[j];}
		}
	}

	if(critbst==MAXDOUBLE) return(criteria());
	if(xlevel.levels[ncol]==1 ||criteria()>critbst)
	{
		for(j=0;j<pmtab.size;j++)
		{
			idx1[j]=xinfo.xid[pibst[j]][ncol];
			idx2[j]=xinfo.xid[pjbst[j]][ncol];
			idx1[pmtab.size+j]=xinfo.xid[nnew-pibst[j]-1][ncol];
			idx2[pmtab.size+j]=xinfo.xid[nnew-pjbst[j]-1][ncol];
		}

		criteria_pm_set(ncol,pmtab.size*2,idx1,idx2);
		xinfo_symm_pm(ncol,pmtab.size,pibst,pjbst);
	}
	FreeVector(idx1);
	FreeVector(idx2);
	FreeVector(pibst);
	FreeVector(pjbst);
	return(critbst);
}


double find_permute(int ncol,int npml,int npme,int nsets)
{
	double perm_obj=0;
	if(xlevel.levels[ncol]>1) get_level_permute(ncol,npml,nsets);

	if(!xinfo.isbalanced[ncol]||xlevel.levels[ncol]!=xinfo.nxvl[ncol])
		perm_obj=get_element_permute(ncol,npme,nsets);
	perm_obj=criteria_set(NULL); set_count=0; 	//prevent accumulating errors
	return(perm_obj);
}


double full_permute(int nsets)
{
	int i,j;
	double perm_obj=0,bst_obj;
	bst_obj=criteria();
	for(i=0;i<nsets;i++)
	{
		save_global();
		for(j=0;j<nactive_col;j++) perm_obj=find_permute(active_col[j],xlevel.levels[j],xlevel.entries[j],1);
		if(perm_obj>bst_obj) restore_global();
		else bst_obj=perm_obj;
	}
	return(bst_obj);
}


int get_ncol(char isrnd)
{
	static int k=-1;
	double prob;
	int i;
	if(!isrnd) {
		#if DEBUG
			//printf("Current selected column is %d\n",active_col[(k+1)%nactive_col]);
		#endif
		return(active_col[(++k)%nactive_col]);
  }
	else
	{
		prob=Random();
		for(i=0;i<nv;i++) if(prob<ncol_prob[i]) break;
		#if DEBUG
			//printf("Current selected column is %d\n",active_col[i]);
		#endif
		return(active_col[i]);
	}
}


std::vector<double> soat_search(double **x)
{
	double old_global_obj,new_obj,th,prob_ratio,alpha1=1.2,alpha2=0.8,cprob=0,prob=0; //fixed
	int i,ncol,nochange=0;
	std::vector<double> return_list;
	clock_t start;

	start = clock();

	th=th0;
	new_obj=obj;
	factor=0.8;
	return_list.push_back(new_obj);
	#if DEBUG
		// for(i=0;i<nv;i++)
		// 	printf("In search.cpp soat_search: xlevel.levels[%d] = %d\n",i,xlevel.levels[i]);
		// printf("In search.cpp soat_search: maxcol = %d\n",maxcol);
	#endif
	while(global_obj-goal>goal*EPS2&&niter<maxiter&&(double)(clock()-start) < maxtime)
	{
		niter++;
		old_global_obj=global_obj;
		for(i=0,prob = 0;i<maxcol;i++)
		{
			#if DEBUG
				//printf("In search.cpp soat_search: nepairs[0] = %d, tol = %f\n",nepairs[0],tol);
			#endif
			ncol=get_ncol(israndcol);
			save_before_change(ncol);
			new_obj=find_exchange(ncol,obj);
			cprob = 1-MIN( 1, MAX( 0, (new_obj-obj)*1.0/th));
			prob = prob + cprob;
			if(Random() >= cprob) restore_after_change(ncol);
			else
			{
				return_list.push_back(new_obj);
				obj=new_obj;
				if(global_obj-obj>minchange)
				{
					global_obj=obj; criteria_x(x);
				}
			}

//////////////////////should be deleted when released////////////////////
/*			 if(return_list.size() >= 300000){
			   return_list.push_back(global_obj);
			   return(return_list);
			 }
*//////////////////////should be deleted when released////////////////////

		}
		if(old_global_obj-global_obj<tol) nochange=1;
		else nochange=0;


		prob_ratio = prob/MAX(10,maxcol);
		if(prob_ratio < hits_ratio) th = th*alpha1;
		else th = th*alpha2;

		/*
			xinfo.nxvl[0] is level.
		 	If latin hyper Cube Design, then adaptive element-wise selection
		*/
		if(xinfo.nxvl[0] == nsamp){
			if(nochange)	updateNepairs(alpha2);
			else 			updateNepairs(alpha1);
		}
	}

	return_list.push_back(global_obj);
	return(return_list);
}



std::vector<double> search(double **x)
{
  std::vector<double> final_obj;
  double objt;
	ntotal=0;
	niter=0;

	global_obj=obj=criteria();
	criteria_x(x);
	objt=ABS(obj)+EPS;

	//could extend to other algorithms here
	tol=dCheckValue(1.0e-15,0.1,objt*5.0e-6,tol);
	minchange=tol*0.01;
	th0=dCheckValue(0,MAXDOUBLE,0.005*objt,th0);
	thmin=th0/100; thmax=100*th0;
	factor=dCheckValue(0.3,0.999,0.8,factor);
	final_obj=soat_search(x);

	return(final_obj);
}

void save_global()
{
    criteria_full_snap();
	xinfo_full_snap();
}


void restore_global()
{
	criteria_full_reset();
	xinfo_full_reset();
}

long get_ntotal()
{
	return(maxpairs*maxcol*niter);
}


void updateNepairs(double alpha)
{
	int i,run = nnew,level = xinfo.nxvl[0],nexvl = nnew/xinfo.nxvl[0];
	int J = run*nexvl*(level-1)/2;
	#if DEBUG
		//printf("In search.cpp, updateNepairs : run = %d; level = %d; factor = %d,nexvl = %d; J = %d\n",run,level,factor,nexvl,J);
	#endif
	for(i=0; i<nv; i++)
	{
		nepairs[i] = (int)(nepairs[i]*alpha);
		nepairs[i] = MAX((int)(J*0.2),nepairs[i]);
		nepairs[i] = MIN((int)(J*0.3),nepairs[i]);
		nepairs[i] = MIN(maxpairs,nepairs[i]);
		/*nepairs[i] = MAX(40,nepairs[i]);
		nepairs[i] = MIN(50,nepairs[i]);*/
	}
}
