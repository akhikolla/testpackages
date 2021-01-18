#include "run.h"

using namespace Rcpp ;
using namespace std ;
using namespace Eigen ;

RcppExport SEXP semR(SEXP X,SEXP m,SEXP K,SEXP Qsem,SEXP Bsem,SEXP Ql,SEXP Bl,SEXP RjSE,SEXP RjM,SEXP maxTry,SEXP run,SEXP detail)
{
  int g=as<int>(K),runC(as<int>(run));
  vector<int> mC=as<vector<int> >(m);
  SEMparameters param;
  param.nGibbsSE = as<vector<int> >(RjSE);
  param.nGibbsM = as<vector<int> >(RjM);
  param.maxIt = as<int>(Qsem);
  param.burnAlgo = as<int>(Bsem);
  param.nGibbsL = as<int>(Ql);
  param.burnL = as<int>(Bl);
  param.maxTry = as<int>(maxTry);
  param.detail = as<bool>(detail);
  
  NumericMatrix XR(X);
  int n(XR.nrow()),col(XR.ncol());
  vector<vector<int> > data(n,vector<int> (col));
  for(int i(0);i<n;i++)
  {
    for(int j(0);j<col;j++)
    {
      data[i][j]=XR[i+j*n];
    }
  }
  
  
  RankCluster semgibbs(data,g,mC,param);
  if(!semgibbs.dataOk())
  {
    vector<double> stock(1,2);
    return List::create(Named("stock")=wrap(stock),Named("indexPb")=wrap(semgibbs.indexPb()));  
  }
  
  semgibbs.run();
  
  //multiple run
  if(runC>1)
  {
    RankCluster semgibbsb(data,g,mC,param);
    double L(-numeric_limits<double>::max());
    if(semgibbs.convergence())
    L=semgibbs.L();
    
    for(int i(1);i<runC;i++)
    {
      semgibbsb.run();
      if(semgibbsb.convergence())
      {
        if(semgibbsb.L()>L)
        {
          L=semgibbsb.L();
          semgibbs=semgibbsb;
        }
      }
    }
  }
  
  if(semgibbs.convergence())
  {
    vector<double> stock(11);
    stock[0]=1;
    stock[1]=semgibbs.partial();
    stock[2]=semgibbs.L();
    stock[3]=semgibbs.bic();
    stock[4]=semgibbs.icl();
    int d=mC.size();
    int n=XR.nrow();
    vector<vector<vector<int> > > data(d,vector<vector<int> > (n));
    vector<vector<vector<double> > > scorePart(d,vector<vector<double> > (n));
    
    for(int i(0);i<n;i++)
    {
      for(int dim(0);dim<d;dim++)
      {
        data[dim][i]=semgibbs.rank(dim,i);
        scorePart[dim][i].resize(mC[dim],1);
      }
    }
    
    vector<vector<vector<int> > > dataInit(data);
    
    if(semgibbs.partial())
    {
      vector<vector<vector<int> > > dataInitSEM(semgibbs.initialPartialRank());
      vector<vector<int> > indexPartial=semgibbs.indexPartialData();
      vector<vector<vector<double> > > scorePartial(semgibbs.partialRankScore());
      
      for(int dim(0);dim<d;dim++)
      {
        int compteur(0);
        for(vector<int>::iterator it=indexPartial[dim].begin();it!=indexPartial[dim].end();it++)
        {
          dataInit[dim][*it]=dataInitSEM[dim][compteur];
          scorePart[dim][*it]=scorePartial[dim][compteur];        
          compteur++;
        }
      }
    }
    
    return List::create(
      //parameters
      Named("stock")=wrap(stock),
      Named("referenceRank")=wrap(semgibbs.mu()),
      Named("p")=wrap(semgibbs.p()),
      Named("proportion")=wrap(semgibbs.proportion()),
      Named("cluster")=wrap(semgibbs.z()),
      Named("tik")=convertMatrix<NumericMatrix,ArrayXXd>(semgibbs.tik()),
      Named("entropy")=convertvector<NumericVector,ArrayXd>(semgibbs.entropy()),
      Named("partialRank")=wrap(data),
      Named("probability")=convertMatrix<NumericMatrix,ArrayXXd>(semgibbs.probability()),
      //distance
      Named("distP")=wrap(semgibbs.distP()),
      Named("distMu")=wrap(semgibbs.distMu()),
      Named("distProp")=wrap(semgibbs.distProp()),
      Named("distZ")=wrap(semgibbs.distZ()),
      Named("distPartialRank")=wrap(semgibbs.distPartialRank()),
      //initialization
      Named("initMu")=wrap(semgibbs.initialMu()),
      Named("initZ")=wrap(semgibbs.initialZ()),
      Named("initPi")=wrap(semgibbs.initialP()),
      Named("initPartialRank")=wrap(dataInit),
      Named("initProportion")=wrap(semgibbs.initialProportion()),
      Named("scorePartial")=wrap(scorePart)//,
      //Named("indexPartialData")=wrap(semgibbs.indexPartialData())
      );
  }
  else
  {
    vector<double> stock(1,0);
    return List::create(Named("stock")=wrap(stock));
  }
}
