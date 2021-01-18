#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace Rcpp;
#include "fct_kv_EvaluesQ.h"
#include "grplassoTEMP.h"
#include "penaliseTEMP.h"
SEXP calc_q(NumericVector Y, List Kv,List k_v, int qmax, double mumin,double muMax, int Num, int maxIter, double eps,
               bool verbose){
  int n; n = Y.size();

  // Error control
  if(mumin>muMax){
    Rcpp::stop("Error : mumin could not be greater than mumax, smaller mumin should be choseen.");
  }
  if(mumin<(muMax/1000)){
    Rcpp::stop("Error : mumin is too small, greater mumin should be choseen.");
  }
  List l;
  List L1;

  int q;
  NumericVector muv(2);
  muv[0] = muMax;
  muv[1] = mumin;
  double mu_new;
  int i=0;
  NumericVector mus;
  NumericVector qs;
  do{
    i+=1;
    mu_new = (muv[0]+muv[1])/2;
    mus.push_back(mu_new);
    L1 = grplasso(Y,Kv,k_v,mu_new,maxIter,eps,verbose);
    NumericVector supp; supp = L1["supp"];
    q = supp.size();
    qs.push_back(q);
    if(q==qmax){
      l = List::create(Named("mus",mus),Named("qs",qs),Named("mu",mu_new),Named("res",L1));
    }
    if(q>qmax){muv[0]=muv[0]; muv[1]=mu_new;}
    if(q<qmax){muv[0]=mu_new; muv[1]=muv[1];}
  } while (q!=qmax&& i<=Num);
  if (i>Num) {
    L1 = grplasso(Y,Kv,k_v,mumin,maxIter,eps,verbose);
    NumericVector supp; supp = L1["supp"];
    q = supp.size();
    mus.push_back(mumin);
    qs.push_back(q);
    if(q==qmax){
      l = List::create(Named("mus",mus),Named("qs",qs),Named("mu",mumin),Named("res",L1));
    }
  }
  if(q==qmax){
    return l;
  }
  if(q>qmax){
    Rcpp::stop("Error : Greater value for Num should be choseen.");
  }
  if(q<qmax){
    Rcpp::stop("Error : Smaller value for mumin should be choseen.");
  }
  return 0;
}//End calc_q

// [[Rcpp::export]]
SEXP RKHSMetMod_qmax(NumericVector Y, NumericMatrix X, String kernel,
                     int Dmax, NumericVector gamma, int qmax, double rat, int Num,
                     bool verbose=false){
  List Zpdme; Zpdme = calc_Kv(X,kernel, Dmax, true,verbose, 1e-08);
  List matZ; matZ = Zpdme[0];
  List resK;
  StringVector namG; namG = Zpdme[1];
  int vMax; vMax = namG.size();
  NumericVector zerosv(vMax);
  List k_v(vMax);
  int n; n = Y.size();
  double meany; meany = mean(Y);
  NumericVector y; y = Y-meany;
  VectorXd def; def = as<VectorXd>(y);
  NumericVector muTEMP(vMax);
  for(int v=0;v<vMax;v++){
    resK = matZ[v];
    NumericVector d; d = resK["Evalues"];
    MatrixXd Q; Q = resK["Q"];
    NumericMatrix D; D = diag(d);
    MatrixXd kv; kv = Q * as<MatrixXd>(D) * Q.transpose();
    k_v[v] = kv;
    double nrm; nrm = def.transpose()*kv*def;
    muTEMP(v) = 2*(sqrt(nrm/n));
  }

  double mumax = max(muTEMP);

  List qres; qres = calc_q(Y,Zpdme,k_v,qmax,mumax/rat,mumax, Num, 1000,1e-4,verbose);
  double mu_q; mu_q = qres["mu"];
  List grpres(1); grpres[0] = qres["res"];
  NumericVector mu(1); mu[0] = mu_q;
  List penMeMod; penMeMod = penMetaMod_cpp(Y,matZ,k_v,namG,grpres,gamma,mu,zerosv,zerosv,1000,
                                           verbose,false);
  NumericVector mus; mus = qres["mus"];
  NumericVector qs; qs = qres["qs"];
  return List::create(Named("mus",(mus)/sqrt(n)),Named("qs",qs),Named("MetaModel",penMeMod));
}//End RKHSMetMod_qmax
