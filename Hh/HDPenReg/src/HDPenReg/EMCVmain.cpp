#include <RTKpp.h>

#include "lassoModels/CVLasso.h"
#include "lassoModels/CVFusedLasso.h"


#include <iostream>

using namespace Rcpp;
using namespace STK;
using namespace std;
using namespace HD;

RcppExport SEXP cvEMlassoMain( SEXP data, SEXP response
                             , SEXP lambda
                             , SEXP nbFolds
                             , SEXP intercept
                             , SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //convert parameters
  int maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
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

  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
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

  // Start CV method
  Residuals measure;

  //create cv for lasso
  CVLasso<Lasso> lassocv;
  //set data
  lassocv.setX(x);
  lassocv.setY(y);
  //set cv parameters
  lassocv.setNbFolds(nbFoldsC);
  lassocv.setIndex(lambdaC);
  //set em parameters
  lassocv.setBurn(burnC);
  lassocv.setMaxStep(maxStepC);
  lassocv.setEps(epsC);
  //set CG parameter
  lassocv.setEpsCG(epsCGC);
  //set threshold for lassosolver
  lassocv.setThreshold(thresholdC);
  //set type of measure
  lassocv.setTypeMeasure(&measure);
  //initialize the class
  lassocv.initialize();

  //run cv
  lassocv.run2();

  //find the position of the lambda with the smallest cv error

  int pos;
  STK::Real minCV = lassocv.cv().minElt(pos);

  // create list
  return List::create( Named("lambda") = wrap(lambdaC)
                     , Named("cv")     = STK::wrap(lassocv.cv())
                     , Named("cvError")= STK::wrap(lassocv.cvError())
                     , Named("minCV")  = wrap(minCV)
                     , Named("lambda.optimal")=wrap(lambdaC[pos])
                     );
}


RcppExport SEXP cvEMfusedLasso1DMain( SEXP data, SEXP response
                                    , SEXP lambda1, SEXP lambda2, SEXP optimL1
                                    , SEXP nbFolds
                                    , SEXP intercept
                                    , SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //convert parameters
  int  maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept), optimL1C=as<bool>(optimL1);
  //convert parameters
  STK::RMatrix<double> dataC(data);
  STK::RVector<double> responseC(response);
  //
  int p = dataC.sizeCols(), n = dataC.sizeRows();
  ArrayXX x = dataC;
  VectorX y = responseC;

  vector<STK::Real> lambda1C = as<vector<STK::Real> >(lambda1);
  vector<STK::Real> lambda2C = as<vector<STK::Real> >(lambda2);
  //center the data if intercept
  STK::Real mu = 0.;
  if(interceptC)
  {
    mu = y.mean();
    y -= mu;
    x -= STK::Const::Vector<STK::Real>(x.rows()) * STK::Stat::mean(x);
  }

  //if lambda1 has to be optimized, we can generate the lambda1 sequence
  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(optimL1C && (lambda1C[0] == -1) )
  {
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =(x.transpose() * y).abs().maxElt();
    STK::Real lambdaMinRatio =(n < p) ? 0.01 : 0.0001;
    STK::Real minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambda1C.resize(100);
    lambda1C[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambda1C[i-1] = exp(log(lambda1C[i]) - gapLambda);
  }
  // lambda sequence to optimize
  vector<STK::Real> lambdaC = (optimL1C) ? lambda1C : lambda2C;

  //create the 1Dcv for fused lasso
  CVFusedLasso1D<FusedLasso> fusedlassocv;
  Residuals measure;
  //set the data
  fusedlassocv.setX(x);
  fusedlassocv.setY(y);
  //set the cv parameter
  fusedlassocv.setNbFolds(nbFoldsC);
  fusedlassocv.setIndex(lambdaC);
  if(optimL1C)
    fusedlassocv.setLambda(lambda2C[0]);//if we have to optimize lambda1, we set lambda2 as lambda parameter
  else
    fusedlassocv.setLambda(lambda1C[0]);
  //set EM parameters
  fusedlassocv.setBurn(burnC);
  fusedlassocv.setEps(epsC);
  fusedlassocv.setMaxStep(maxStepC);
  //set fusedlasso solver parameter
  fusedlassocv.setThreshold(thresholdC);
  //set CG parameter
  fusedlassocv.setEpsCG(epsCGC);
  fusedlassocv.setTypeMeasure(&measure);

  //initialize the created class
  fusedlassocv.initialize();

  //run the CV
  fusedlassocv.run2();

  //find the position of the lambda with the smallest cv error
  int pos;
  STK::Real minCV = fusedlassocv.cv().minElt(pos);
  STK::Real lambdaMin = (optimL1C) ? lambda1C[pos] : lambda2C[pos];

  return List::create( Named("lambda")  = wrap(lambdaC)
                     , Named("cv")      = STK::wrap(fusedlassocv.cv())
                     , Named("cvError") = STK::wrap(fusedlassocv.cvError())
                     , Named("minCV")   = wrap(minCV)
                     , Named("lambda.optimal") = wrap(lambdaMin)
                     );
}

RcppExport SEXP cvEMfusedLasso2DMain( SEXP data, SEXP response
                                    , SEXP lambda1, SEXP lambda2
                                    , SEXP nbFolds
                                    , SEXP intercept
                                    , SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //convert parameters
  int  maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept);
  //convert parameters
  STK::RMatrix<double> dataC(data);
  STK::RVector<double> responseC(response);
  //
  //int p = dataC.sizeCols(), n = dataC.sizeRows();
  ArrayXX x = dataC;
  VectorX y = responseC;
  vector<STK::Real> lambda1C = as<vector<STK::Real> >(lambda1);
  vector<STK::Real> lambda2C = as<vector<STK::Real> >(lambda2);
  //center the data if intercept
  STK::Real mu = 0.;
  if(interceptC)
  {
    mu = y.mean();
    y -= mu;
    x -= STK::Const::Vector<STK::Real>(x.rows()) * STK::Stat::mean(x);
  }
//  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
//  if(lambda1C[0] == -1)
//  {
//    //corrMax is the first lambda value in the lars sequence
//    STK::Real corrMax =(x.transpose() * y).abs().maxElt();
//    STK::Real lambdaMinRatio = (n < p) ? 0.01 : 0.0001;
//    STK::Real minLambda = corrMax * lambdaMinRatio;
//    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
//    lambda1C.resize(100);
//    lambda1C[99] = corrMax;
//    for(int i = 99; i>0 ; i--)
//      lambda1C[i-1] = exp(log(lambda1C[i]) - gapLambda);
//  }

  //create the CV for fused lasso
  CVFusedLasso2D<FusedLasso> fusedlassocv;
  Residuals measure;
  //set the data
  fusedlassocv.setX(x);
  fusedlassocv.setY(y);
  //set CV parameters
  fusedlassocv.setNbFolds(nbFoldsC);
  fusedlassocv.setIndex(lambda1C);
  fusedlassocv.setIndexL2(lambda2C);
  //set EM parameters
  fusedlassocv.setEps(epsC);
  fusedlassocv.setMaxStep(maxStepC);
  fusedlassocv.setBurn(burnC);
  //set fusedsolver parameters
  fusedlassocv.setThreshold(thresholdC);
  //set CG parameter
  fusedlassocv.setEpsCG(epsCGC);
  fusedlassocv.setTypeMeasure(&measure);
  //initialize the class
  fusedlassocv.initialize();
  //run the algo
  fusedlassocv.run2();

  //find the position of the lambda with the smallest error
  int pos;
  STK::Real minCV = fusedlassocv.cv().minElt(pos);
  //convert the position in lambda1 and lambda2
  vector<STK::Real> lambdaMin(2);
  lambdaMin[0] = lambda1C[pos/lambda2C.size()];
  lambdaMin[1] = lambda2C[pos%lambda2C.size()];

  return List::create( Named("cv") = STK::wrap(fusedlassocv.cv())
                     , Named("cvError") = STK::wrap(fusedlassocv.cvError())
                     , Named("minCV") = wrap(minCV)
                     , Named("lambda.optimal")=wrap(lambdaMin)
                     );
}



// duplication pour cv logistic, trouver autre chose de mieux

RcppExport SEXP cvEMlogisticLassoMain( SEXP data, SEXP response
                                     , SEXP lambda, SEXP nbFolds
                                     , SEXP intercept
                                     , SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //convert parameters
  int maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept);
  STK::RMatrix<double> dataC(data);
  STK::RVector<double> responseC(response);
  vector<STK::Real> lambdaC = as<vector<STK::Real> >(lambda);
  //
  int p = dataC.sizeCols(), n = dataC.sizeRows();
  ArrayXX x = dataC;
  VectorX y = responseC;

  if(interceptC) interceptC=true;
  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(lambdaC[0] == -1)
  {
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =(x.transpose() * y).abs().maxElt();
    STK::Real lambdaMinRatio = (n < p) ? 0.01 : 0.0001;
    STK::Real minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambdaC.resize(100);
    lambdaC[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambdaC[i-1] = exp(log(lambdaC[i]) - gapLambda);
  }
  AUC measure;
  //create cv for lasso
  CVLasso<LogisticLasso> lassocv;
  //set data
  lassocv.setX(x);
  lassocv.setY(y);
  //set cv parameters
  lassocv.setNbFolds(nbFoldsC);
  lassocv.setIndex(lambdaC);
  //set em parameters
  lassocv.setBurn(burnC);
  lassocv.setMaxStep(maxStepC);
  lassocv.setEps(epsC);
  //set CG parameter
  lassocv.setEpsCG(epsCGC);
  //set threshold for lassosolver
  lassocv.setThreshold(thresholdC);
  //set type of measure
  lassocv.setTypeMeasure(&measure);
  //initialize the class
  lassocv.initialize();

  //run cv
  lassocv.run2();

  //find the position of the lambda with the smallest cv error
  int pos;
  STK::Real minCV = lassocv.cv().minElt(pos);

  return List::create( Named("lambda")  = wrap(lambdaC)
                     , Named("cv")      = STK::wrap(lassocv.cv())
                     , Named("cvError") = STK::wrap(lassocv.cvError())
                     , Named("minCV")   = wrap(minCV)
                     , Named("lambda.optimal") = wrap(lambdaC[pos]));
}



RcppExport SEXP cvEMlogisticFusedLasso1DMain( SEXP data, SEXP response
                                            , SEXP lambda1, SEXP lambda2, SEXP optimL1
                                            , SEXP nbFolds
                                            , SEXP intercept
                                            , SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //convert parameters
  int  maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  bool interceptC=as<bool>(intercept), optimL1C=as<bool>(optimL1);
  //convert parameters
  STK::RMatrix<double> dataC(data);
  STK::RVector<double> responseC(response);
  //
  int p = dataC.sizeCols(), n = dataC.sizeRows();
  ArrayXX x = dataC;
  VectorX y = responseC;

  vector<STK::Real> lambda1C = as<vector<STK::Real> >(lambda1);
  vector<STK::Real> lambda2C = as<vector<STK::Real> >(lambda2);

  if(interceptC) interceptC=true;
  //if lambda1 has to be optimized, we can generate the lambda1 sequence
  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
  if(optimL1C && (lambda1C[0] == -1) )
  {
    int pos;
    //corrMax is the first lambda value in the lars sequence
    STK::Real corrMax =(x.transpose() * y).abs().maxElt(pos);
    STK::Real lambdaMinRatio = (n < p) ? 0.01 : 0.0001;
    STK::Real minLambda = corrMax * lambdaMinRatio;
    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
    lambda1C.resize(100);
    lambda1C[99] = corrMax;
    for(int i = 99; i>0 ; i--)
      lambda1C[i-1] = exp(log(lambda1C[i]) - gapLambda);
  }
  // lambda sequence to optimize
  vector<STK::Real> lambdaC;
  lambdaC = (optimL1C) ? lambda1C : lambda2C;
  //create the 1Dcv for fused lasso
  CVFusedLasso1D<LogisticFusedLasso> fusedlassocv;
  AUC measure;
  //set the data
  fusedlassocv.setX(x);
  fusedlassocv.setY(y);
  //set the cv parameter
  fusedlassocv.setNbFolds(nbFoldsC);
  fusedlassocv.setIndex(lambdaC);
  if(optimL1C)
    fusedlassocv.setLambda(lambda2C[0]);//if we have to optimize lambda1, we set lambda2 as lambda parameter
  else
    fusedlassocv.setLambda(lambda1C[0]);
  //set EM parameters
  fusedlassocv.setBurn(burnC);
  fusedlassocv.setEps(epsC);
  fusedlassocv.setMaxStep(maxStepC);
  //set fusedlasso solver parameter
  fusedlassocv.setThreshold(thresholdC);
  fusedlassocv.setTypeMeasure(&measure);
  //set CG parameter
  fusedlassocv.setEpsCG(epsCGC);

  //initialize the created class
  fusedlassocv.initialize();

  //run the CV
  fusedlassocv.run2();

  //find the position of the lambda with the smallest cv error
  int pos;
  STK::Real minCV = fusedlassocv.cv().minElt(pos);
  STK::Real lambdaMin = (optimL1C) ? lambda1C[pos] : lambda2C[pos];

  return List::create( Named("lambda")  = wrap(lambdaC)
                     , Named("cv")      = STK::wrap(fusedlassocv.cv())
                     , Named("cvError") = STK::wrap(fusedlassocv.cvError())
                     , Named("minCV")   = wrap(minCV)
                     , Named("lambda.optimal") = wrap(lambdaMin)
                     );
}

RcppExport SEXP cvEMlogisticFusedLasso2DMain( SEXP data, SEXP response
                                            , SEXP lambda1, SEXP lambda2
                                            , SEXP nbFolds
                                            , SEXP intercept
                                            , SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
  //convert parameters
  int  maxStepC(as<int>(maxStep)), burnC(as<int>(burn)), nbFoldsC(as<int>(nbFolds));
  Real epsC(as<Real>(eps)), thresholdC(as<Real>(threshold)), epsCGC(as<Real>(epsCG));
  //bool interceptC=as<bool>(intercept);
  //convert parameters
  STK::RMatrix<double> dataC(data);
  STK::RVector<double> responseC(response);
  //
  //int p = dataC.sizeCols(), n = dataC.sizeRows();
  ArrayXX x = dataC;
  VectorX y = responseC;
  vector<STK::Real> lambda1C = as<vector<STK::Real> >(lambda1);
  vector<STK::Real> lambda2C = as<vector<STK::Real> >(lambda2);

//  if(interceptC) interceptC=true;
//  //if lambdaC[0]=-1, we have to generate the lambda sequence with the same way as the glmnet package
//  if(lambda1C[0] == -1)
//  {
//    //corrMax is the first lambda value in the lars sequence
//    STK::Real corrMax =(x.transpose() * y).abs().maxElt();
//    STK::Real lambdaMinRatio = (n < p) ? 0.01 : 0.0001;
//    STK::Real minLambda = corrMax * lambdaMinRatio;
//    STK::Real gapLambda = (log(corrMax)-log(minLambda))/99;
//    lambda1C.resize(100);
//    lambda1C[99] = corrMax;
//    for(int i = 99; i>0 ; i--)
//      lambda1C[i-1] = exp(log(lambda1C[i]) - gapLambda);
//  }

  //create the CV for fused lasso
  CVFusedLasso2D<LogisticFusedLasso> fusedlassocv;
  AUC measure;
  //set the data
  fusedlassocv.setX(x);
  fusedlassocv.setY(y);
  //set CV parameters
  fusedlassocv.setNbFolds(nbFoldsC);
  fusedlassocv.setIndex(lambda1C);
  fusedlassocv.setIndexL2(lambda2C);
  //set EM parameters
  fusedlassocv.setEps(epsC);
  fusedlassocv.setMaxStep(maxStepC);
  fusedlassocv.setBurn(burnC);
  //set fusedsolver parameters
  fusedlassocv.setThreshold(thresholdC);
  //set CG parameter
  fusedlassocv.setEpsCG(epsCGC);
  fusedlassocv.setTypeMeasure(&measure);
  //initialize the class
  fusedlassocv.initialize();

  //run the algo
  fusedlassocv.run2();

  //find the position of the lambda with the smallest error
  int pos;
  STK::Real minCV = fusedlassocv.cv().minElt(pos);
  //convert the position in lambda1 and lambda2
  vector<STK::Real> lambdaMin(2);
  lambdaMin[0] = lambda1C[pos/lambda2C.size()];
  lambdaMin[1] = lambda2C[pos%lambda2C.size()];

  return List::create( Named("cv")      = STK::wrap(fusedlassocv.cv())
                     , Named("cvError") = STK::wrap(fusedlassocv.cvError())
                     , Named("minCV")   = wrap(minCV)
                     , Named("lambda.optimal")=wrap(lambdaMin)
                     );
}

