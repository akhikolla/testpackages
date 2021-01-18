#include "larsRmain.h"



#ifdef LARS_DEBUG
#include <iostream>
#endif

using namespace Rcpp;
using namespace STK;
using namespace std;
using namespace HD;

template<class Vector>
void convertToVector(SEXP const& rVector, Vector &output)
{
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering convertToVector")<<endl;
#endif
  STK::RVector<STK::Real> data(rVector);
#ifdef STK_BOUNDS_CHECK
  if (data.size() != output.size())
  {
    stk_cerr << _T("In convertToVector, sizes mismatch:\n");
    stk_cerr << _T("data.size() =") << data.size() << _T("\n");
    stk_cerr << _T("output.size() =") << output.size() << _T("\n");
  }
#endif
  for( int i =output.begin(), idata=data.begin(); i< output.end(); ++i, ++idata)
    output[i]=data[idata];
}

template<class Array>
void convertToArray(SEXP const& rMatrix, Array &output)
{
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering convertToArray")<<endl;
#endif
  STK::RMatrix<STK::Real> data(rMatrix);
#ifdef STK_BOUNDS_CHECK
  if (data.sizeRows() != output.sizeRows())
  {
    stk_cerr << _T("In convertToArray, sizes rows mismatch:\n");
    stk_cerr << _T("data.sizeRows() =") << data.sizeRows() << _T("\n");
    stk_cerr << _T("output.sizeRows() =") << output.sizeRows() << _T("\n");
  }
  if (data.sizeCols() != output.sizeCols())
  {
    stk_cerr << _T("In convertToArray, sizes cols mismatch:\n");
    stk_cerr << _T("data.sizeCols() =") << data.sizeCols() << _T("\n");
    stk_cerr << _T("output.sizeCols() =") << output.sizeCols() << _T("\n");
  }
#endif
  for(int i=output.beginRows(), iData = data.beginRows(); i<output.endRows(); ++i, ++iData)
    for(int j=output.beginCols(), jData = data.beginCols();j<output.endCols(); ++j, ++jData)
      output(i,j)=data(iData, jData);
}

RcppExport SEXP larsmain( SEXP data, SEXP response
                        , SEXP nbIndiv, SEXP nbVar
                        , SEXP maxStep, SEXP intercept, SEXP eps)
{
#ifdef LARS_DEBUG
  stk_cerr << _T("Entering larsmain")<<std::endl;
#endif
  //convert parameters
  int p = Rcpp::as<int>(nbVar), n = Rcpp::as<int>(nbIndiv), maxStepC = Rcpp::as<int>(maxStep);
  bool interceptC = Rcpp::as<bool>(intercept);
  STK::Real epsC  = Rcpp::as<STK::Real>(eps);

  STK::CArrayXX x(STK::Range(1,n), STK::Range(1,p));
  STK::CVectorX y(STK::Range(1,n));
  convertToArray(data,x);
  convertToVector(response,y);

#ifdef LARS_DEBUG
  stk_cerr << _T("larsmain. Creating Lars")<<endl;
#endif
  Lars lars(x,y,maxStepC,interceptC,epsC);
  lars.run();
#ifdef LARS_DEBUG
  stk_cerr << _T("larsmain. Lars.run() done")<<endl;
#endif

  int step=lars.step();
  vector<double> l1norm(step+1);
  vector<vector<int> > varIdx(step+1);
  vector<vector<double> > varCoeff(step+1);
  vector<vector<int> > evoIdxDrop(step);
  vector<vector<int> > evoIdxAdd(step);
  l1norm[0]=0;
  for(int i=1;i<=step;i++)
  {
    varIdx[i].resize(lars.path(i).size());
    varCoeff[i].resize(lars.path(i).size());
    for(int j=1;j<=lars.path(i).size();j++)
    {
      varCoeff[i][j-1]=lars.coefficient(i,j);
      varIdx[i][j-1]=lars.varIdx(i,j);
    }
    l1norm[i]=lars.l1norm(i);
    if(lars.evolution()[i-1].first.size()!=0)
        evoIdxAdd[i-1]=lars.evolution()[i-1].first;
    if(lars.evolution()[i-1].second.size()!=0)
        evoIdxDrop[i-1]=lars.evolution()[i-1].second;

  }
#ifdef LARS_DEBUG
  stk_cerr << _T("larsmain done")<<std::endl;
#endif

  return List::create( Named("l1norm")    =wrap(l1norm)
                     , Named("lambda")    =wrap(lars.lambda())
                     , Named("varIdx")    =wrap(varIdx)
                     , Named("varCoeff")  =wrap(varCoeff)
                     , Named("evoDropIdx")=wrap(evoIdxDrop)
                     , Named("evoAddIdx") =wrap(evoIdxAdd)
                     , Named("step")      =wrap(lars.step())
                     , Named("mu")        =wrap(lars.mu())
                     , Named("ignored")   =STK::wrap(lars.toIgnore().cast<int>())
                     , Named("error")     =wrap(lars.msg_error())
                     , Named("muX")       =STK::wrap(lars.muX())
                     );
}

RcppExport SEXP fusionmain(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps)
{
#ifdef FUSION_DEBUG
  stk_cerr << _T("Entering fusionmain")<<endl;
#endif
  //double t1,t2;
  //t1=clock();
  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv)), maxStepC(as<int>(maxStep));
  bool interceptC = as<bool>(intercept);
  Real epsC = as<STK::Real>(eps);

  STK::CArrayXX x(STK::Range(1,n),STK::Range(1,p));
  STK::CVectorX y(STK::Range(1,n));
  convertToArray(data, x);
  convertToVector(response,y);

  //t2=clock();
  //cout<<"Temps conversion des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //run algorithm
  //t1=clock();
#ifdef FUSION_DEBUG
  stk_cerr << _T("fusionmain. Creating Fusion")<<endl;
#endif
  Fusion fusion(x,y,maxStepC,interceptC,epsC);
  fusion.run();
#ifdef FUSION_DEBUG
  stk_cerr << _T("fusionmain. fusion.run() done")<<endl;
#endif
  //t2=clock();

  //cout<<"Temps fusion:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  //extract and convert results
  //t1=clock();
  int step=fusion.step();

  vector<double> l1norm(step+1);
  vector<vector<int> > varIdx(step+1);
  vector<vector<double> > varCoeff(step+1);
  vector<vector<int> > evoIdxDrop(step);
  vector<vector<int> > evoIdxAdd(step);

  l1norm[0]=0;
  for(int i = 1; i <= step; i++)
  {
    varIdx[i].resize(fusion.path(i).size());
    varCoeff[i].resize(fusion.path(i).size());
    for(int j = 1; j <= fusion.path(i).size(); j++)
    {
      varCoeff[i][j-1]=fusion.coefficient(i,j);
      varIdx[i][j-1]=fusion.varIdx(i,j);
    }
    l1norm[i]=fusion.l1norm(i);
    if(fusion.evolution()[i-1].first.size()!=0)
        evoIdxAdd[i-1]=fusion.evolution()[i-1].first;
    if(fusion.evolution()[i-1].second.size()!=0)
        evoIdxDrop[i-1]=fusion.evolution()[i-1].second;
  }

#ifdef FUSION_DEBUG
  stk_cerr << _T("fusionmain done")<<endl;
#endif
  //t2=clock();
  //cout<<"Temps extraction des données:"<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;

  return List::create( Named("l1norm")    =wrap(l1norm)
                     , Named("lambda")    =wrap(fusion.lambda())
                     , Named("varIdx")    =wrap(varIdx)
                     , Named("varCoeff")  =wrap(varCoeff)
                     , Named("step")      =wrap(step)
                     , Named("evoDropIdx")=wrap(evoIdxDrop)
                     , Named("evoAddIdx") =wrap(evoIdxAdd)
                     , Named("mu")        =wrap(fusion.mu())
                     , Named("ignored")   =STK::wrap(fusion.toIgnore().cast<int>())
                     , Named("error")     =wrap(fusion.msg_error())
                     , Named("muX")       =STK::wrap(fusion.muX())
                     );

}

RcppExport SEXP cvlarsmain( SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar
                          , SEXP maxStep, SEXP intercept, SEXP eps, SEXP nbFold
                          , SEXP partition, SEXP index, SEXP mode)
{
#ifdef CVLARS_DEBUG
  stk_cerr << _T("Entering cvlarsmain")<<endl;
#endif
  //double t1,t2;
  //t1=clock();
  //convert parameters
  int p(as<int>(nbVar)), n(as<int>(nbIndiv));
  int maxStepC(as<int>(maxStep)), nbFoldC(as<int>(nbFold));
  bool interceptC = as<bool>(intercept);
  bool modeLambda = as<bool>(mode);
  STK::Real epsC(as<STK::Real>(eps));

  vector<double> indexC=as<vector<double> >(index);
  vector<int> partitionC=as<vector<int> >(partition);

  STK::CArrayXX x(STK::Range(1,n), STK::Range(1,p));
  STK::CVectorX y(STK::Range(1,n));
  convertToArray(data,x);
  convertToVector(response,y);
  //t2=clock();
  //run algorithm
  //t1=clock();
#ifdef CVLARS_DEBUG
  stk_cerr << _T("cvlarsmain. Creating Cvlars")<<endl;
#endif
  Cvlars cvlars(x,y,nbFoldC,indexC,modeLambda,maxStepC,interceptC,epsC);
  if(partitionC[0]!=-1) { cvlars.setPartition(partitionC);}

#ifdef _OPENMP
  cvlars.run2();
#ifdef CVLARS_DEBUG
  stk_cerr << _T("cvlarsmain. Cvlars.run2() done")<<endl;
#endif
#else
  cvlars.run();
#ifdef CVLARS_DEBUG
  stk_cerr << _T("cvlarsmain. Cvlars.run() done")<<endl;
#endif
#endif
  //t2=clock();
#ifdef CVLARS_DEBUG
  stk_cerr << _T("cvlarsmain done")<<endl;
#endif
  return List::create( Named("cv")     =STK::wrap(cvlars.cv())
                     , Named("cvError")=STK::wrap(cvlars.cvError()));
}
