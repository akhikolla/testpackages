#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <string>
#include "groupLasso.h"
#include "LUfit.h"
#include "pgLUfit.h"
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;
using Rcpp::as;

// // ' @export
//[[Rcpp::export]]
Rcpp::List LU_dense_cpp(Eigen::Map<Eigen::MatrixXd> X_, 
                        Eigen::VectorXd & z_, 
                        Eigen::VectorXd & icoef_,
                        Eigen::ArrayXd & gsize_,
                        Eigen::ArrayXd & pen_,
                        Eigen::ArrayXd & lambdaseq_,
                        bool user_lambdaseq_,
                        int pathLength_,
                        double lambdaMinRatio_,
                        double pi_,
                        int max_nUpdates_,
                        int maxit_,
                        Eigen::VectorXd & wei_,
                        bool weiOption_,
                        double tol_,
                        double inner_tol_,
                        bool useStrongSet_,
                        bool verbose_, 
                        double stepSize_,
                        double stepSizeAdj_,
                        int batchSize_,
                        int updateFreq_,
                        std::vector<double> samplingProbabilities_,
                        bool useLipschitz_,
                        std::string method_,
                        int trace_,
                        bool skipFitting_)
{
  
  try{
    
    if(method_=="CD"){
      LUfit<Eigen::Map<Eigen::MatrixXd> > lu(X_,z_,icoef_,gsize_,pen_,
                                             lambdaseq_,user_lambdaseq_,pathLength_,
                                             lambdaMinRatio_,pi_,max_nUpdates_,maxit_,wei_,weiOption_,tol_,
                                             inner_tol_,useStrongSet_,verbose_,trace_);
      if(!skipFitting_){lu.LUfit_main();}
      return Rcpp::List::create(Rcpp::Named("coef") = lu.getCoefficients(),
                                Rcpp::Named("std_coef") = lu.getStdCoefficients(),
                                Rcpp::Named("iters") = lu.getIters(),
                                Rcpp::Named("nUpdates") = lu.getnUpdates(),
                                Rcpp::Named("nullDev") = lu.getnullDev(),
                                Rcpp::Named("deviance") = lu.getDeviances(),
                                Rcpp::Named("lambda")=lu.getLambdaSequence(),
                                Rcpp::Named("convFlag")=lu.getconvFlag(),
                                Rcpp::Named("fVals") = lu.getfVals(),
                                Rcpp::Named("subgrads") = lu.getSubGradients(),
                                Rcpp::Named("fVals_all") = lu.getfVals_all(),
                                Rcpp::Named("beta_all") = lu.getbeta_all(),
                                Rcpp::Named("method") = method_);
      
    }else{
      pgLUfit<Eigen::Map<MatrixXd> > lu(X_,z_,icoef_,gsize_,pen_,lambdaseq_,user_lambdaseq_,pathLength_,
                                        lambdaMinRatio_,pi_,maxit_,tol_,verbose_,
                                        stepSize_,stepSizeAdj_,batchSize_,updateFreq_,samplingProbabilities_,useLipschitz_,method_,trace_);
      if(!skipFitting_){lu.pgLUfit_main();}
      return Rcpp::List::create(Rcpp::Named("coef") = lu.getCoefficients(),
                                Rcpp::Named("std_coef") = lu.getStdCoefficients(),
                                Rcpp::Named("iters") = lu.getIters(),
                                Rcpp::Named("nullDev") = lu.getnullDev(),
                                Rcpp::Named("deviance") = lu.getDeviances(),
                                Rcpp::Named("lambda")=lu.getLambdaSequence(),
                                Rcpp::Named("convFlag")=lu.getconvFlag(),
                                Rcpp::Named("fVals") = lu.getfVals(),
                                Rcpp::Named("subgrads") = lu.getSubGradients(),
                                Rcpp::Named("fVals_all") = lu.getfVals_all(),
                                Rcpp::Named("beta_all") = lu.getbeta_all(),
                                Rcpp::Named("stepSize") = lu.getStepSize(),
                                Rcpp::Named("samplingProbabilities")=lu.getSamplingProbabilities(),
                                Rcpp::Named("method") = method_);
      // return R_NilValue;
    }
  }
  catch(const std::invalid_argument& e ){
    throw std::range_error(e.what());
  }
  return R_NilValue;
}

// //' @export
//[[Rcpp::export]]
Rcpp::List LU_sparse_cpp(Eigen::SparseMatrix<double> & X_, 
                         Eigen::VectorXd & z_, 
                         Eigen::VectorXd & icoef_,
                         Eigen::ArrayXd & gsize_,
                         Eigen::ArrayXd & pen_,
                         Eigen::ArrayXd & lambdaseq_,
                         bool user_lambdaseq_,
                         int pathLength_,
                         double lambdaMinRatio_,
                         double pi_,
                         int max_nUpdates_,
                         int maxit_,
                         Eigen::VectorXd & wei_,
                         bool weiOption_,
                         double tol_,
                         double inner_tol_,
                         bool useStrongSet_, 
                         bool verbose_,
                         double stepSize_,
                         double stepSizeAdj_,
                         int batchSize_,
                         int updateFreq_,
                         std::vector<double> samplingProbabilities_,
                         bool useLipschitz_,std::string method_,
                         int trace_,
                         bool skipFitting_)
{
  try{
    
    if(method_=="CD"){
      LUfit<Eigen::SparseMatrix<double> > lu(X_,z_,icoef_,gsize_,pen_,
                                             lambdaseq_,user_lambdaseq_,pathLength_,
                                             lambdaMinRatio_,pi_,max_nUpdates_,maxit_,wei_,weiOption_,tol_,
                                             inner_tol_,useStrongSet_,verbose_,trace_);
      if(!skipFitting_){lu.LUfit_main();}
      return Rcpp::List::create(Rcpp::Named("coef") = lu.getCoefficients(),
                                Rcpp::Named("std_coef") = lu.getStdCoefficients(),
                                Rcpp::Named("iters") = lu.getIters(),
                                Rcpp::Named("nUpdates") = lu.getnUpdates(),
                                Rcpp::Named("nullDev") = lu.getnullDev(),
                                Rcpp::Named("deviance") = lu.getDeviances(),
                                Rcpp::Named("lambda")=lu.getLambdaSequence(),
                                Rcpp::Named("convFlag")=lu.getconvFlag(),
                                Rcpp::Named("fVals") = lu.getfVals(),
                                Rcpp::Named("subgrads") = lu.getSubGradients(),
                                Rcpp::Named("fVals_all") = lu.getfVals_all(),
                                Rcpp::Named("beta_all") = lu.getbeta_all(),
                                Rcpp::Named("method") = method_);
      
    }else{
      pgLUfit<Eigen::SparseMatrix<double> > lu(X_,z_,icoef_,gsize_,pen_,lambdaseq_,user_lambdaseq_,pathLength_,
                                        lambdaMinRatio_,pi_,maxit_,tol_,verbose_,
                                        stepSize_,stepSizeAdj_,batchSize_,updateFreq_,samplingProbabilities_,useLipschitz_,method_,trace_);
      if(!skipFitting_){lu.pgLUfit_main();}
      return Rcpp::List::create(Rcpp::Named("coef") = lu.getCoefficients(),
                                Rcpp::Named("std_coef") = lu.getStdCoefficients(),
                                Rcpp::Named("iters") = lu.getIters(),
                                Rcpp::Named("nullDev") = lu.getnullDev(),
                                Rcpp::Named("deviance") = lu.getDeviances(),
                                Rcpp::Named("lambda")=lu.getLambdaSequence(),
                                Rcpp::Named("convFlag")=lu.getconvFlag(),
                                Rcpp::Named("fVals") = lu.getfVals(),
                                Rcpp::Named("subgrads") = lu.getSubGradients(),
                                Rcpp::Named("fVals_all") = lu.getfVals_all(),
                                Rcpp::Named("beta_all") = lu.getbeta_all(),
                                Rcpp::Named("stepSize") = lu.getStepSize(),
                                Rcpp::Named("samplingProbabilities")=lu.getSamplingProbabilities(),
                                Rcpp::Named("method") = method_);
      // return R_NilValue;
    }
  }
  catch(const std::invalid_argument& e ){
    throw std::range_error(e.what());
  }
  return R_NilValue;
}

//[[Rcpp::export]]
Eigen::MatrixXd deviances_dense_cpp(Eigen::MatrixXd & coefMat_, 
                                    Eigen::Map<Eigen::MatrixXd> & X_,
                                    Eigen::VectorXd & z_, 
                                    double pi_,
                                    const Eigen::VectorXd & wei_, 
                                    bool weiOption_)
{
  int K = coefMat_.cols();
  VectorXd deviances(K);
  
  for (int j=0;j<K;j++)
  {
    deviances(j) = evalDeviance(X_,z_,pi_,coefMat_.middleCols(j,1),wei_,weiOption_);
  }
  return deviances;
}


//[[Rcpp::export]]
Eigen::MatrixXd deviances_sparse_cpp(Eigen::MatrixXd & coefMat_,
                                     Eigen::SparseMatrix<double> & X_, 
                                     Eigen::VectorXd & z_, 
                                     double pi_,
                                     const Eigen::VectorXd & wei_,
                                     bool weiOption_)
{
  int K = coefMat_.cols();
  VectorXd deviances(K);
  
  for (int j=0;j<K;j++)
  {
    deviances(j) = evalDeviance(X_,z_,pi_,coefMat_.middleCols(j,1),wei_,weiOption_);
  }
  return deviances;
}


