#include <Rcpp.h>
#include "robetsTargetFunction.h"
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List calc_out(std::vector<double> y,int m,std::vector<double> initstate,int errortype,int trendtype,int seasontype,
bool damped,double alpha,double beta,double gamma,double phi,int nmse, double k){

  RobetsTargetFunction* sp = new RobetsTargetFunction();
  
  sp->oneEval(y, errortype, trendtype, seasontype, damped, nmse, m,
		 alpha, beta, gamma, phi,initstate,k);
  
  return List::create(_["lik"] = sp->getLik(),_["amse"] = sp->getAmse(),_["e"] = sp->getE(),_["states"] = sp->getState(),
  _["roblik"] = sp->getRoblik(), _["tau2"] = sp->getTau2());
}

// [[Rcpp::export]]
double tau2(std::vector<double> x){
  
  RobetsTargetFunction* sp = new RobetsTargetFunction();
    
  return sp->computeTau2(x);
}
