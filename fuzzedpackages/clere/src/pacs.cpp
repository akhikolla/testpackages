#include <cstdio>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <RcppEigen.h>

#include "Eigen/Dense"
#include "Eigen/SVD"
#include "Eigen/Core"
#include "conversion.h"

using namespace std;
using namespace Rcpp;
using namespace Eigen;

RcppExport SEXP pacs(SEXP R_IOPACSobj){  
  BEGIN_RCPP  
    Rcpp::S4 obj(R_IOPACSobj);
  // We need to access the following quantities
  // n, p, nItMax, eps, lambda
  // the data x, y and starting beta: bstart

  // Get model dimensions
  int n         = Rcpp::as<int> (obj.slot("n"));
  int p         = Rcpp::as<int> (obj.slot("p"));
  int nItMax    = Rcpp::as<int> (obj.slot("nItMax"));
  double eps    = Rcpp::as<double> (obj.slot("epsPACS"));
  double lambda = Rcpp::as<double> (obj.slot("lambda"));

  // Scalar parameters
  double littleeps = 1e-7;

  // Data
  MatrixXd x = MatrixXd(n, p);
  VectorXd y = VectorXd(n);

  // Primary Parameters
  VectorXd xTy   = VectorXd(p);
  MatrixXd xTx   = MatrixXd(p, p);

  VectorXd betal     = VectorXd(p);
  VectorXd betatilde = VectorXd(p);
  
  // Read the data
  NumericMatrix Rx(SEXP(obj.slot("x")));
  NumericVector Ry(SEXP(obj.slot("y")));
  convertMatrix<NumericMatrix,MatrixXd>(Rx,x);
  convertVector<NumericVector,VectorXd>(Ry,y);

  xTx = x.transpose()*x;
  xTy = x.transpose()*y;

  // Matrix dm and dp
  int q = p*(p-1)/2;
  MatrixXd dm   = MatrixXd::Zero(q, p);
  MatrixXd dp   = MatrixXd::Zero(q, p);
    
  int nrow_min = 0;
  int nrow_max = p-1;
  int l,k,ik;
  l=0; k=1;
  
  while(l<q){
    ik = k + 1;
    if(l<nrow_max and l>=nrow_min){
      dp(l,k-1)  = +1.0;
      dp(l,ik-1) = +1.0;
      dm(l,k-1)  = +1.0;
      dm(l,ik-1) = -1.0;
      l++;
      ik++;
    }
    if(l>=nrow_max){
      nrow_min  = nrow_max;
      nrow_max += p-k-1;
      k++;
    }
  }

  // Read starting parameters
  NumericVector Rbetal(SEXP(obj.slot("betaInput")));
  convertVector<NumericVector,VectorXd>(Rbetal,betal);
  
  int j;
  for(j = 0; j < p; j++){
    betatilde(j) = betal(j);
  }
  
  VectorXd dmb    = VectorXd(q);
  VectorXd dpb    = VectorXd(q);
  dmb             = dm * betatilde;
  dpb             = dp * betatilde;
  
  // Loic rewrites...  
  int p2q = p+2*q;
  int ppq = p+q;
  int l_p, l_ppq;
  VectorXd ascvec = VectorXd(p2q); 
  for(l=0;l<p;l++){
    ascvec(l) = 1.0 / (littleeps+abs(betatilde(l)));
  }
  for(l=p;l<ppq;l++){
    l_p = l-p;
    ascvec(l) = 1.0 / (littleeps+abs( dmb(l_p) ));
  }
  for(l=ppq;l<p2q;l++){
    l_ppq = l-ppq;
    ascvec(l) = 1.0 / (littleeps+abs( dpb(l_ppq) ));
  }
  
  VectorXd mm1 = VectorXd(p);
  VectorXd mm2 = VectorXd(q);
  VectorXd mm3 = VectorXd(q);
  MatrixXd dmd = MatrixXd(p, p);
  MatrixXd dpd = MatrixXd(p, p);
  MatrixXd D   = MatrixXd(p, p);
  VectorXd db  = VectorXd(p);
  int px2 = p * 2;
  int it  = 0;
  double crit = 1.0;
  double tmp;
  
  while(it<nItMax and crit>eps){
    for(j=0;j<p;j++){
      mm1(j) = ascvec(j)/(littleeps+abs(betatilde(j)));
    }
    for(l=0;l<q;l++){
      mm2(l) = ascvec(p+l)/(littleeps+abs(dmb(l)));
      mm3(l) = ascvec(px2+l)/(littleeps+abs(dpb(l)));
    }

    dmd   = lambda * ((dm.transpose() * mm2.asDiagonal()) * dm);
    dpd   = lambda * ((dp.transpose() * mm3.asDiagonal()) * dp);
    
    D     = dmd + dpd;
    
    for(j=0;j<p;j++){
      D(j,j) += lambda * mm1(j);
    }
    D.noalias() += xTx;
    betal = D.colPivHouseholderQr().solve(xTy);
    crit  = 0.0;
    for(j=0;j<p;j++){
      tmp = abs(betal(j)-betatilde(j));
      if(tmp>crit){
	crit = tmp;
      }
    }
    betatilde = betal;
    // cout<<"\tIteration n°"<<it<<" - eps = "<<crit<<endl; => couls be recorded in next version of of the code
    dmb   = dm * betatilde;
    dpb   = dp * betatilde;
    it++;
  }

  // Write output => betatilde
  obj.slot("betaOutput")  = outVector<NumericVector,VectorXd>(betatilde);
  END_RCPP
}

