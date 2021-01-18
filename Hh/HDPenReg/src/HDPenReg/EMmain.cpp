#include <RTKpp.h>

#include "lassoModels/EM.h"

#include "lassoModels/Lasso.h"
#include "lassoModels/FusedLasso.h"
#include "lassoModels/LogisticLasso.h"
#include "lassoModels/LogisticFusedLasso.h"

#include <iostream>

using namespace Rcpp;
using namespace STK;
using namespace std;
using namespace HD;

RcppExport SEXP EMlassoMain( SEXP data, SEXP response
                           , SEXP lambda, SEXP intercept
                           , SEXP maxStep, SEXP burn
                           , SEXP threshold, SEXP eps, SEXP epsCG)
{
  //convert parameters
  int maxStepC(as<int>(maxStep)), burnC(as<int>(burn));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept);
  STK::RMatrix<double> dataC(data);
  STK::RVector<double> responseC(response);
  vector<STK::Real> lambdaC = as<vector<STK::Real> >(lambda);
  //
  int p = dataC.sizeCols(), n = dataC.sizeRows();
  ArrayXX x = dataC;
  VectorX y = responseC;
  //center the data if intercept
  STK::Real mu = 0.;
  if(interceptC)
  {
    mu = y.mean();
    y -= mu;
    x -= STK::Const::Vector<STK::Real>(x.rows()) * STK::Stat::mean(x);
  }

  // if lambdaC[0]=-1, we have to generate the lambda sequence in the same way
  // as the glmnet package
  if(lambdaC[0] == -1)
  {
    lambdaC.resize(100);
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax = (x.transpose() * y).abs().maxElt();
    STK::Real gapLambda = (log(corrMax)-log(corrMax * (n<p) ? 0.01 : 0.0001))/99;
    lambdaC[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambdaC[i-1] = exp(log(lambdaC[i]) - gapLambda);
  }

  //create EM
  EM algo(maxStepC,burnC,epsC);
  Lasso lasso( &x, &y, lambdaC[0], thresholdC, epsCGC);
  //run for all lambda
  vector<int> step;
  vector<double> logLikelihood;
  Rcpp::List pathCoefficients;
  Rcpp::List pathIndex;
  for(int i = 0; i < (int) lambdaC.size(); i++)
  {
#ifdef HD_DEBUG
    std::cout << "\nIn EMLasso. Set lambda =" << lambdaC[i] << "\n";
#endif
    //change the lambda
    lasso.setLambda(lambdaC[i]);
    //and initialize the solver
    lasso.initializeBeta();
    //run algorithm
    if (!algo.run(&lasso))
    {

#ifdef HD_DEBUG
      std::cout << "\nAn error occur in algo.run(&lasso).\nWhat: " << algo.error() << "\n";
#endif
    }
    // save results
    pathCoefficients.push_back(STK::wrap(lasso.currentBeta()));
    pathIndex.push_back(STK::wrap(lasso.currentSet()+1));
    step.push_back(algo.step());
    logLikelihood.push_back(lasso.lnLikelihood());
    //if there is 0 non-zeros coefficients, we stop at this values of lambda
    if( (lasso.currentBeta().size()==0))
    {
      lambdaC.erase(lambdaC.begin()+i+1,lambdaC.end());
      break;
    }
  }
  // return result
  return List::create( Named("variable")=wrap(pathIndex)
                     , Named("coefficient")=wrap(pathCoefficients)
                     , Named("lambda")=wrap(lambdaC)
                     , Named("mu")=wrap(mu)
                     , Named("logLikelihood") = wrap(logLikelihood)
                     , Named("step")=wrap(step)
                     );
}


RcppExport SEXP EMfusedLassoMain( SEXP data, SEXP response
                                , SEXP lambda1, SEXP lambda2
                                , SEXP intercept, SEXP maxStep
                                , SEXP burn
                                , SEXP eps, SEXP eps0, SEXP epsCG)
{
  //convert parameters
  int  maxStepC(as<int>(maxStep)), burnC(as<int>(burn));
  Real epsC(as<Real>(eps)), eps0C(as<Real>(eps0)), epsCGC(as<Real>(epsCG))
     , lambda1C(as<Real>(lambda1)), lambda2C(as<Real>(lambda2));
  bool interceptC=as<bool>(intercept);
  STK::RMatrix<double> dataC(data);
  STK::RVector<double> responseC(response);
  //
  //int p = dataC.sizeCols(), n = dataC.sizeRows();
  ArrayXX x = dataC;
  VectorX y = responseC;

  //center the data if intercept
  STK::Real mu = 0.;
  if(interceptC)
  {
    mu = y.mean();
    y -= mu;
    x -= STK::Const::Vector<STK::Real>(x.rows()) * STK::Stat::mean(x);
  }
  //create EM
  EM algo(maxStepC,burnC,epsC);
  //create fused lasso
#ifdef HD_DEBUG
  std::cout << "In EMfusedLassoMain\n"
            << "Creating fusedlasso class with\n"
            << "lambda1 =" << lambda1C << "\n"
            << "lambda2 =" << lambda2C << "\n"
            << "eps0 =" << eps0C << "\n"
            << "epsCGC =" << epsCGC << "\n";
#endif
  FusedLasso fusedlasso(&x,&y, lambda1C, lambda2C, eps0C, epsCGC);
  //run
  algo.run(&fusedlasso);
  //number of step made by the em
  return List::create( Named("coefficient")=wrap(fusedlasso.beta())
                     , Named("lambda1")=wrap(lambda1C)
                     , Named("lambda2")=wrap(lambda2C)
                     , Named("mu")=wrap(mu)
                     , Named("logLikelihood") = wrap(fusedlasso.lnLikelihood())
                     , Named("step")=wrap(algo.step())
                     );
}

RcppExport SEXP EMlogisticLassoMain( SEXP data, SEXP response
                                   , SEXP lambda, SEXP intercept
                                   , SEXP maxStep, SEXP burn
                                   , SEXP threshold, SEXP eps, SEXP epsCG)
{
  //convert parameters
  int maxStepC(as<int>(maxStep)), burnC(as<int>(burn));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  //bool interceptC=as<bool>(intercept);
  STK::RMatrix<double> dataC(data);
  STK::RVector<double> responseC(response);
  vector<STK::Real> lambdaC = as<vector<STK::Real> >(lambda);

  //
  int p = dataC.sizeCols(), n = dataC.sizeRows();
  ArrayXX x = dataC;
  VectorX y = responseC;
  //center the data if intercept
  STK::Real mu = 0.;
  // if lambdaC[0]=-1, we have to generate the lambda sequence in the same way
  // as the glmnet package
  if(lambdaC[0] == -1)
  {
    lambdaC.resize(100);
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax = (x.transpose() * y).abs().maxElt();
    STK::Real gapLambda = (log(corrMax)-log(corrMax * (n<p) ? 0.01 : 0.0001))/99;
    lambdaC[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambdaC[i-1] = exp(log(lambdaC[i]) - gapLambda);
  }
  //create EM
  EM algo(maxStepC,burnC,epsC);
#ifdef HD_DEBUG
  std::cout << "Creating LogisticLasso(x,y," << lambdaC[0] << ", " << thresholdC << ", " << epsCGC << ")" << std::endl;
#endif
  LogisticLasso lasso( &x, &y, lambdaC[0], thresholdC, epsCGC);
  //run for all lambda
  vector<int> step;
  vector<double> logLikelihood;
  Rcpp::List pathCoefficients;
  Rcpp::List pathIndex;
  for(int i = 0; i < (int) lambdaC.size(); i++)
  {
    //change the lambda
    lasso.setLambda(lambdaC[i]);
    //run algorithm
    algo.run(&lasso);
    // save results
    pathCoefficients.push_back(STK::wrap(lasso.currentBeta()));
    pathIndex.push_back(STK::wrap(lasso.currentSet()+1));
    step.push_back(algo.step());
    logLikelihood.push_back(lasso.lnLikelihood());
    //if there is 0 non-zeros coefficients, we stop at this values of lambda
    if( (lasso.currentBeta().size()==0))
    {
      lambdaC.erase(lambdaC.begin()+i+1,lambdaC.end());
      break;
    }

    //reinitialize the solver for next lambda
    lasso.initializeBeta();
  }
  // return result
  return List::create( Named("variable")=wrap(pathIndex)
                     , Named("coefficient")=wrap(pathCoefficients)
                     , Named("lambda")=wrap(lambdaC)
                     , Named("mu")=wrap(mu)
                     , Named("logLikelihood") = wrap(logLikelihood)
                     , Named("step")=wrap(step)
                     );
}

RcppExport SEXP EMlogisticFusedLassoMain( SEXP data, SEXP response
                                        , SEXP lambda1, SEXP lambda2
                                        , SEXP intercept, SEXP maxStep
                                        , SEXP burn
                                        , SEXP eps, SEXP eps0, SEXP epsCG)
{
  //convert parameters
  int  maxStepC(as<int>(maxStep)), burnC(as<int>(burn));
  Real epsC(as<Real>(eps)), eps0C(as<Real>(eps0)), epsCGC(as<Real>(epsCG))
     , lambda1C(as<Real>(lambda1)), lambda2C(as<Real>(lambda2));
  bool interceptC=as<bool>(intercept);
  STK::RMatrix<double> dataC(data);
  STK::RVector<double> responseC(response);
  //
  //int p = dataC.sizeCols(), n = dataC.sizeRows();
  ArrayXX x = dataC;
  VectorX y = responseC;

  //center the data if intercept
  STK::Real mu = 0.;
  //create EM
  EM algo(maxStepC,burnC,epsC);
#ifdef HD_DEBUG
  std::cout << "In EMlogisticFusedLassoMain\n"
            << "Creating LogisticFusedLasso class with\n"
            << "lambda1 =" << lambda1C << "\n"
            << "lambda2 =" << lambda2C << "\n"
            << "eps0 =" << eps0C << "\n"
            << "epsCGC =" << epsCGC << "\n";
#endif
  LogisticFusedLasso logisticfusedlasso(&x,&y, lambda1C, lambda2C, eps0C, epsCGC);
  //run
  algo.run(&logisticfusedlasso);
  //number of step made by the em
  return List::create( Named("coefficient")=wrap(logisticfusedlasso.beta())
                     , Named("lambda1")=wrap(lambda1C)
                     , Named("lambda2")=wrap(lambda2C)
                     , Named("mu")=wrap(mu)
                     , Named("logLikelihood") = wrap(logisticfusedlasso.lnLikelihood())
                     , Named("step")=wrap(algo.step())
                     );
}


