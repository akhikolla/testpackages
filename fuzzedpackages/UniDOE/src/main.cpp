#include "doe_criteria.h"
#include "doe_search.h"
#include "doe_Matrix.h"
#include "doe_Eval.h"
#include "doe_utility.h"
#include "doe_DesignInfo.h"
#include <string>
#include <Rcpp.h>
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>// std::random_shuffle
#include <vector>// std::vector
#include <ctime> // std::time
using namespace Rcpp;
using namespace std;
//CD2 related criteria. options.isnorm should be between 0~2. i.e. options.isnorm=cCheckValue(0,2,2,options.isnorm);

double **x;
static int RAND_STEP = 999;
SEARCHOPT options;
CRITOPT critopt;

void create_options(int nv);
void create_critopt(int nv);
void free_options();
void free_critopt();
void initialize_pars(int nv);
int check_pars(int nv,int nnew);
NumericMatrix Generate_init_matrix(StringVector opt, int n, int s, int q,NumericMatrix initX);
int myrandom (int i);


// [[Rcpp::export]]
List SATA_UD(int n, int s, int q, StringVector init, NumericMatrix initX, int crit, int maxiter, double hits_ratio)
{
  List lst;
  int nsamp,nv,np=0,i,j,orth;
  double critobj,**rang=NULL,critobj0,search_time;
  std::vector<double> critobj_vector;
  clock_t start;
  NumericMatrix InputX,return_matrix;
  InputX = Generate_init_matrix(init,n,s,q,initX);
	nsamp = InputX.nrow();
	nv = InputX.ncol();
	x = NewDMatrix(nsamp,nv);
  for(j=0;j<nv;j++)
  {
    for(i=0;i<nsamp;i++) x[i][j] = InputX(i,j);
  }

  create_options(nv);
	create_critopt(nv);
  initialize_pars(nv);
	critopt.type = crit;
  options.maxiter=maxiter;
  options.hits_ratio = hits_ratio;

  if(initX.ncol()>1&& as<string>(init) == "orth") orth=initX.ncol();
  else orth=0;
  for(i=0;i<nv;i++) options.colweight[i]=i<orth?0:options.colweight[i];//{ critopt.scale[i]=i<orth?0:critopt.scale[i];}
  check_pars(nv,nsamp);
  //Some input judgements are omitted;
  create_xinfo(x,nsamp,np,nv);
  return_matrix = NumericMatrix(nsamp,nv);//Matlaissymmb: plhs[0]=mxCreateDoubleMatrix(nnew,nv,mxREAL);
  if(options.isnorm) rang=normalize(options.isnorm,x,nsamp,nv);//,xinfo.nxvl);
	create_criteria(x,nsamp,np,nv,&critopt);
	critobj0=criteria();
	create_search(x,nsamp,np,nv,&options);
	start = clock();
	critobj_vector = search(x);
	search_time = (double)(clock()-start)/CLOCKS_PER_SEC;
	critobj=critobj_vector.back();
	if(options.isnorm)
  {
    unnormalize(x,rang,nsamp,nv);
    FreeDMatrix(rang);
  }
	for(i=0;i<nsamp; i++)
  {
      for(j=0;j<nv;j++)  return_matrix(i,j) = x[i][j];
  }
	//printf("\nIn main.cpp, StoUDC function, maxcol = %d\n",options.maxcol);

  FreeDMatrix(x);
	free_xinfo();
	free_criteria();
	free_search();
  free_critopt();
	free_options();

	lst["Init_Matrix"] = InputX;
	lst["UniDOE_Matrix"] = return_matrix;
	lst["obj0"] = critobj0;
	lst["obj"] = critobj;
	lst["time(s)"]= search_time;
	lst["obj_list"] = wrap(critobj_vector);
	return lst;
}

// [[Rcpp::export]]
List SATA_AUD( NumericMatrix XP,int n, int s, int q, StringVector init, NumericMatrix initX, int crit, int maxiter, double hits_ratio)
{
  int nv=s,nnew=n,nsamp,np=XP.nrow(),i,j,orth;
  double critobj,**rang=NULL,critobj0,search_time;
  std::vector<double> critobj_vector;
  clock_t start;
  NumericMatrix InputX(n,s),Init_matrix(n+np,s),return_matrix;
  List lst;
  InputX = Generate_init_matrix(init,n,s,q,initX);
  nsamp=np+nnew;
  x = NewDMatrix(nsamp,nv);
  for(j=0;j<nv;j++)
  {
    for(i=0;i<np;i++) {x[i][j] =Init_matrix(i,j)= XP(i,j);}
    for(i=0;i<nnew;i++) {x[i+np][j] = Init_matrix(i+np,j)=InputX(i,j);}
  }

  create_options(nv);
  create_critopt(nv);
  initialize_pars(nv);
  critopt.type = crit;
  options.israndcol=0;
  options.maxiter=maxiter;
  options.hits_ratio = hits_ratio;

  if(initX.ncol()>1&&as<string>(init) == "orth") orth=initX.ncol();
  else orth=0;
  for(i=0;i<nv;i++) options.colweight[i]=i<orth?0:options.colweight[i];
  check_pars(nv,nnew);
  //Some input judgements are omitted;
  create_xinfo(x,nnew,np,nv); //nnew = nsamp here
  return_matrix = NumericMatrix(nsamp,nv);//Matlab: plhs[0]=mxCreateDoubleMatrix(nnew,nv,mxREAL);
  if(options.isnorm) rang=normalize(options.isnorm,x,nsamp,nv);//,xinfo.nxvl);

  create_criteria(x,nnew,np,nv,&critopt);
  critobj0=criteria();
  create_search(x,nnew,np,nv,&options);
  start = clock();
  critobj_vector = search(x);
  search_time = (double)(clock()-start)/CLOCKS_PER_SEC;
  critobj=critobj_vector.back();
  if(options.isnorm) {unnormalize(x,rang,nsamp,nv); FreeDMatrix(rang);}
  for(i=0;i<nsamp; i++)
  {
    for(j=0;j<nv;j++)  return_matrix(i,j) = x[i][j];
  }

  FreeDMatrix(x);
  free_xinfo();
  free_criteria();
  free_search();
  free_critopt();
  free_options();

  lst["Init_Matrix"] = Init_matrix;
  lst["UniDOE_Matrix"] = return_matrix;
  lst["obj0"] = critobj0;
  lst["obj"] = critobj;
  lst["time(s)"]= search_time;
  lst["obj_list"] = wrap(critobj_vector);
  return lst;
}


// [[Rcpp::export]]
List SATA_LP(NumericMatrix X0, int q, int crit, int maxiter, double hits_ratio)
{
  List lst;
  clock_t start;
  int nsamp,nv,np=0,i,j;
  double critobj,**rang=NULL,critobj0,search_time;
  std::vector<double> critobj_vector;
  NumericMatrix InputX,return_matrix;
  InputX = X0;
  nsamp = InputX.nrow();
  nv = InputX.ncol();
  x = NewDMatrix(nsamp,nv);
  for(j=0;j<nv;j++)
  {
    for(i=0;i<nsamp;i++) x[i][j] = InputX(i,j);
  }

  create_options(nv);
  create_critopt(nv);
  initialize_pars(nv);
  critopt.type = crit;
  options.maxiter=maxiter;
  options.hits_ratio = hits_ratio;
  for(i=0;i<nv;i++) options.levels[i] = q;//level permutation;

  check_pars(nv,nsamp);
  create_xinfo(x,nsamp,np,nv); //nnew = nsamp here
  return_matrix = NumericMatrix(nsamp,nv);
  if(options.isnorm) rang=normalize(options.isnorm,x,nsamp,nv);
  create_criteria(x,nsamp,np,nv,&critopt);
  critobj0=criteria();
  create_search(x,nsamp,np,nv,&options);
  start = clock();
  critobj_vector = search(x);
  search_time = (double)(clock()-start)/CLOCKS_PER_SEC;
  critobj=critobj_vector.back();
  if(options.isnorm)
  {
    unnormalize(x,rang,nsamp,nv);
    FreeDMatrix(rang);
  }
  for(i=0;i<nsamp; i++)
  {
    for(j=0;j<nv;j++)  return_matrix(i,j) = x[i][j];
  }

  FreeDMatrix(x);
  free_xinfo();
  free_criteria();
  free_search();
  free_critopt();
  free_options();

  lst["Init_Matrix"] = InputX;
  lst["UniDOE_Matrix"] = return_matrix;
  lst["obj0"] = critobj0;
  lst["obj"] = critobj;
  lst["time(s)"]= search_time;
  lst["obj_list"] = wrap(critobj_vector);
  return lst;

}


void create_options(int nv)
{
	options.levels=NewIVector(nv);
	options.colweight=NewDVector(nv);
}

void create_critopt(int nv)
{
	critopt.scale=NewDVector(nv);
}


void initialize_pars(int nv)
{
	int i;
	critopt.type=-1;	    critopt.ismax=-1;
	critopt.npars[0]=-1;   critopt.npars[1]=-1;
	critopt.npars[2]=-1;
	critopt.pars[0]=NULL;  critopt.pars[1]=NULL;
	critopt.pars[2]=NULL;
	critopt.func[0]=-1;	   critopt.goal=-1.0e50;
	for(i=0;i<nv;i++) critopt.scale[i]=-1;

	options.maxcol=-1;     options.maxpairs=-1;
	options.maxtime=-1;    options.th0=-1;
	options.factor=-1;     options.maxiter=-1;
	options.isnorm=-1;		 /*options.issymm=-1;*/
	options.israndcol=1;   options.israndpairs=-1; //options.israndcol=-1; -> options.israndcol = 1; 2017/nov/1
	options.tol=-1;        options.hits_ratio = -1; // hits_ratio used in soat algorithm
	options.isperm=-1;
	for(i=0;i<nv;i++)
  {
    options.levels[i]=-1;
    options.colweight[i]=-1;
  }

}

void free_options()
{
	FreeVector(options.levels);
	FreeVector(options.colweight);
}

void free_critopt()
{
	FreeVector(critopt.scale);
}


int check_pars(int nv,int nnew)
{
	int i;
	critopt.type=cCheckValue(1,4,3,critopt.type); //HAVE BEEN CHANGED FROM cCheckValue(0,4,1,critopt.type)
	critopt.ismax=cCheckValue(0,1,0,critopt.ismax);
	critopt.npars[0]=iCheckValue(0,200,0,critopt.npars[0]);
	critopt.npars[1]=iCheckValue(0,200,0,critopt.npars[1]);
	critopt.npars[2]=iCheckValue(0,200,0,critopt.npars[2]);
	if(critopt.ismax) critopt.goal=dCheckValue(-1.0e20,1.0e20,-1.0e20,-critopt.goal);
	else critopt.goal=dCheckValue(-1.0e20,1.0e20,-1.0e20,critopt.goal);
	for(i=0;i<nv;i++) critopt.scale[i]=dCheckValue(1.0e-40,1,1,critopt.scale[i]);

	options.maxtime=dCheckValue(1.0e-4,10000000,100000,options.maxtime);//options.maxtime=dCheckValue(1.0e-4,1000000,10000,options.maxtime);
  options.hits_ratio = dCheckValue(0,1,0.1,options.hits_ratio);
  //options.maxpairs is J in paper, i.e. number of exchanges searched in one inner iteration.
	options.maxpairs=iCheckValue(1,nnew*(nnew-1)/2,50,options.maxpairs);

	if(critopt.type==2 || critopt.type==3 || critopt.type==4) options.isnorm=cCheckValue(0,2,2,options.isnorm);//HAVE BEEN CHANGED FROM critopt.type==2
  else options.isnorm=cCheckValue(0,2,1,options.isnorm);
	//options.issymm=cCheckValue(0,1,0,options.issymm);
	options.israndcol=cCheckValue(0,1,0,options.israndcol);
	options.israndpairs=cCheckValue(0,1,1,options.israndpairs);
	options.isperm=cCheckValue(0,1,1,options.isperm);


	for(i=0;i<nv;i++)
	{
		options.colweight[i]=dCheckValue(0,1,1,options.colweight[i]);
		if(critopt.scale[i]<EPS2) options.colweight[i]=0;
	}
  //if(options.issymm) for(i=0;i<nv;i++) options.levels[i]=1;
	//else for(i=0;i<nv;i++) options.levels[i]=iCheckValue(1,nnew,1,options.levels[i]);
  for(i=0;i<nv;i++) options.levels[i]=iCheckValue(1,nnew,1,options.levels[i]);

	return(0);
}

NumericMatrix Generate_init_matrix(StringVector opt, int n, int s, int q, NumericMatrix initX)
{
  int i,j;
  //std::srand( unsigned (std::time(0)) );
  std::vector<int> col;
  NumericMatrix return_matrix  = NumericMatrix(n,s);
  string option = as<string>(opt);

  for(i=1;i<=n;i++) col.push_back( (i%q)+1 );
  for(i=0;i<s;i++)
  {
    std::random_shuffle ( col.begin(), col.end(),myrandom );
    for(j=0;j<n;j++)  return_matrix(j,i) = col[j];
  }
  if(option == "orth" && initX.ncol()>1)
  {
    for(i=0;i<initX.ncol();i++)
    for(j=0;j<initX.nrow();j++) return_matrix(j,i) = initX(j,i);
  }
  else if(option == "input" && initX.nrow()>1)
  {
    return (initX);
  }

  return return_matrix;
}

int myrandom (int i)
{
  RAND_STEP = ((RAND_STEP+3)*7) % 1000000;
  return (RAND_STEP%i);
}
