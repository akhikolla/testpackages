#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace Rcpp;
#include "grplassoTEMP.h"
// [[Rcpp::export]]
SEXP grplasso_q(NumericVector Y, List Kv, int q, double rat, int Num){
  double eps=1e-4;
  int n; n = Y.size();
  List matZ; matZ = Kv[0];
  List resK(2);
  StringVector namesGrp; namesGrp = Kv[1];
  int vMax = namesGrp.size();
  List k_v(vMax);
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
    k_v[v]=kv;
    double nrm; nrm = def.transpose()*kv*def;
    muTEMP(v) = 2*(sqrt(nrm/n));
  }
  double mumax = max(muTEMP);
  // Error control
  double mumin; mumin = mumax/rat;
  if(mumin>mumax){
    Rcpp::stop("Error : mumin could not be greater than mumax, smaller mumin should be choseen.");
  }
  if(mumin<(mumax/1000)){
    Rcpp::stop("Error : mumin is too small, greater mumin should be choseen.");
  }
  List l;
  List L1;

    int qg;
    NumericVector muv(2);
      muv[0] = mumax;
      muv[1] = mumin;
      double mu_new;
      int i=0;
      NumericVector mus;
      NumericVector qs;
      do{
        i+=1;
        mu_new = (muv[0]+muv[1])/2;
        mus.push_back(mu_new);
        L1 = grplasso(Y,Kv,k_v,mu_new,1000,eps,false);
        NumericVector supp; supp = L1["supp"];
        qg = supp.size();
        qs.push_back(qg);
        if(qg==q){
          l = List::create(Named("mus",mus),Named("qs",qs),Named("mu",mu_new),Named("res",L1));
        }
        if(qg>q){muv[0]=muv[0]; muv[1]=mu_new;}
        if(qg<q){muv[0]=mu_new; muv[1]=muv[1];}
      } while (qg!=q&& i<=Num);
      if (i>Num) {
        L1 = grplasso(Y,Kv,k_v,mumin,1000,eps,false);
        NumericVector supp; supp = L1["supp"];
        qg = supp.size();
        mus.push_back(mumin);
        qs.push_back(qg);
        if(qg==q){
          l = List::create(Named("mus",mus),Named("qs",qs),Named("mu",mumin),Named("res",L1));
        }
      }
      if(qg==q){
        return l;
      }
      if(qg>q){
        Rcpp::stop("Error : Greater value for Num should be choseen.");
      }
      if(qg<q){
        Rcpp::stop("Error : Smaller value for mumin should be choseen.");
      }
  return 0;
}//End grplasso_q
