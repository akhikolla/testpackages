// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
List bess_lm(Eigen::MatrixXd X, Eigen::VectorXd y, int T0, int max_steps, Eigen::VectorXd beta0) {
  int n=X.rows();
  int p=X.cols();
  double max_T=0.0;
  vector<int>E(p);
  for(int k=0;k<=p-1;k++) {
    E[k]=k;
  }
  vector<int>I(p-T0);
  vector<int>A(T0);
  vector<int>J(p-T0);
  vector<int>B(T0);
  Eigen::MatrixXd X_A = Eigen::MatrixXd::Zero(n, T0);
  Eigen::MatrixXd X_I = Eigen::MatrixXd::Zero(n, p-T0);
  Eigen::VectorXd beta_A = Eigen::VectorXd::Zero(T0);
  Eigen::VectorXd d_I = Eigen::VectorXd::Zero(p-T0);
  Eigen::VectorXd beta=beta0;
  Eigen::VectorXd d=(X.transpose()*(y-X*beta)) /double(n);
  Eigen::VectorXd bd=beta+d;
  bd=bd.cwiseAbs();
  for(int k=0;k<=T0-1;k++) {
    max_T=bd.maxCoeff(&A[k]);
    bd(A[k])=0.0;
  }
  sort (A.begin(),A.end());
  set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  for(int l=1;l<=max_steps;l++) {
    for(int mm=0;mm<=T0-1;mm++) {
      X_A.col(mm)=X.col(A[mm]);
    }
    for(int mm=0;mm<=p-T0-1;mm++) {
      X_I.col(mm)=X.col(I[mm]);
    }
    beta_A=X_A.colPivHouseholderQr().solve(y);
    for(int mm=0;mm<=T0-1;mm++) {
      d(A[mm])=0.0;
      beta(A[mm])=beta_A(mm);
    }
    d_I=(X_I.transpose()*(y-X_A*beta_A))/double(n);
    for(int mm=0;mm<=p-T0-1;mm++) {
      beta(I[mm])=0.0;
      d(I[mm])=d_I(mm);
    }
    bd=beta+d;
    bd=bd.cwiseAbs();
    for(int k=0;k<=T0-1;k++) {
      max_T=bd.maxCoeff(&B[k]);
      bd(B[k])=0;
    }
    sort (B.begin(),B.end());
    set_difference(E.begin(),E.end(), B.begin(),B.end(),J.begin());
    if(A==B) { break;} else {
      A=B;
      I=J;
    }
   }
   double mse=(y-X*beta).array().square().sum();
   mse=mse/(2*n);
   //return beta;
   return List::create(Named("beta")=beta,Named("max_T")=max_T,Named("mse")=mse);
}
// [[Rcpp::export]]
List get_A(Eigen::MatrixXd& X, Eigen::VectorXd& y, Eigen::VectorXd& beta, double& coef0,int& T0, Eigen::VectorXi& B, Eigen::VectorXd& weights) {
  double max_T=0.0;
  int n=X.rows();
  int p=X.cols();
  Eigen::VectorXd one_xbeta_exp = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd bd = Eigen::VectorXd::Zero(p);
  Eigen::VectorXi A_out = Eigen::VectorXi::Zero(T0);
  Eigen::VectorXi I_out = Eigen::VectorXi::Zero(p-T0);
  Eigen::VectorXd one = Eigen::VectorXd::Ones(n);
  vector<int>E(p);
  for(int k=0;k<=p-1;k++) {
    E[k]=k;
  }
  vector<int>A(T0);
  vector<int>I(p-T0);
  Eigen::VectorXd coef(n);
  for(int i=0;i<=n-1;i++) {
    coef(i)=coef0;
  }
  Eigen::VectorXd xbeta_exp = X*beta+coef;
  for(int i=0;i<=n-1;i++) {
    if(xbeta_exp(i)>25.0) xbeta_exp(i) = 25.0;
    if(xbeta_exp(i)<-25.0) xbeta_exp(i) = -25.0;
  }
  xbeta_exp = xbeta_exp.array().exp();
  Eigen::VectorXd pr = xbeta_exp.array()/(xbeta_exp+one).array();
  Eigen::VectorXd l1=-X.adjoint()*((y-pr).cwiseProduct(weights));
  //Eigen::MatrixXd X2=X.adjoint().array().square();
  X=X.array().square();
  Eigen::VectorXd l2=(X.adjoint())*((pr.cwiseProduct(one-pr)).cwiseProduct(weights));
  Eigen::VectorXd d=-l1.cwiseQuotient(l2);
  if(B.size()<p) {
    for(int k=0;k<=B.size()-1;k++) {
      d(B(k)-1)=0.0;
    }
  }
  bd = (beta+d).cwiseAbs().cwiseProduct(l2.cwiseSqrt());
  //bd=bd.cwiseAbs();
  //bd=bd.cwiseProduct(l2.cwiseSqrt());
  for(int k=0;k<=T0-1;k++) {
    max_T=bd.maxCoeff(&A[k]);
    bd(A[k])=0.0;
  }
  sort (A.begin(),A.end());
  set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  for(int i=0;i<T0;i++)
	  A_out(i) = A[i]+1;
  for(int i=0;i<p-T0;i++)
	  I_out(i) = I[i]+1;
  return List::create(Named("p")=pr,Named("A")=A_out,Named("I")=I_out,Named("max_T")=max_T);
}
// [[Rcpp::export]]
List getcox_A(Eigen::MatrixXd& X, Eigen::VectorXd& beta, int& T0, Eigen::VectorXi& B, Eigen::VectorXd& status, Eigen::VectorXd& weights) {
  double max_T=0.0;
  int n=X.rows();
  int p=X.cols();
  Eigen::VectorXi A_out = Eigen::VectorXi::Zero(T0);
  Eigen::VectorXi I_out = Eigen::VectorXi::Zero(p-T0);
  Eigen::VectorXd l1 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd l2 = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd cum_theta=Eigen::VectorXd::Zero(n);
  Eigen::VectorXd d = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd bd = Eigen::VectorXd::Zero(p);
  Eigen::MatrixXd xtheta(n,p);
  Eigen::MatrixXd x2theta(n,p);
  vector<int>E(p);
  for(int k=0;k<=p-1;k++) {
    E[k]=k;
  }
  vector<int>A(T0);
  vector<int>I(p-T0);
  Eigen::VectorXd theta=X*beta;
    for(int i=0;i<=n-1;i++) {
    if(theta(i)>25.0) theta(i) = 25.0;
    if(theta(i)<-25.0) theta(i) = -25.0;
  }
  theta=weights.array()*theta.array().exp();
  cum_theta(n-1)=theta(n-1);
  for(int k=n-2;k>=0;k--) {
    cum_theta(k)=cum_theta(k+1)+theta(k);
  }
  for(int k=0;k<=p-1;k++) {
    xtheta.col(k)=theta.cwiseProduct(X.col(k));
  }
  for(int k=0;k<=p-1;k++) {
    x2theta.col(k)=X.col(k).cwiseProduct(xtheta.col(k));
  }
  for(int k=n-2;k>=0;k--) {
    xtheta.row(k)=xtheta.row(k+1)+xtheta.row(k);
  }
  for(int k=n-2;k>=0;k--) {
    x2theta.row(k)=x2theta.row(k+1)+x2theta.row(k);
  }
  for(int k=0;k<=p-1;k++) {
    xtheta.col(k)=xtheta.col(k).cwiseQuotient(cum_theta);
  }
  for(int k=0;k<=p-1;k++) {
    x2theta.col(k)=x2theta.col(k).cwiseQuotient(cum_theta);
  }
  x2theta=x2theta.array()-xtheta.array().square().array();
  xtheta=X.array()-xtheta.array();
  if(status.size()>0) {
    for(int k=0;k<=status.size()-1;k++) {
      xtheta.row(status(k)-1)=Eigen::VectorXd::Zero(p);
      x2theta.row(status(k)-1)=Eigen::VectorXd::Zero(p);
    }
  }
  l1=-xtheta.adjoint()*weights;
  l2=x2theta.adjoint()*weights;
  d=-l1.cwiseQuotient(l2);
  if(B.size()<p) {
    for(int k=0;k<=B.size()-1;k++) {
      d(B(k)-1)=0.0;
    }
  }
  bd=beta+d;
  bd=bd.cwiseAbs();
  bd=bd.cwiseProduct(l2.cwiseSqrt());
  for(int k=0;k<=T0-1;k++) {
    max_T=bd.maxCoeff(&A[k]);
    bd(A[k])=0.0;
  }
  sort (A.begin(),A.end());
  set_difference(E.begin(),E.end(), A.begin(),A.end(),I.begin());
  for(int i=0;i<T0;i++)
	  A_out(i) = A[i]+1;
  for(int i=0;i<p-T0;i++)
	  I_out(i) = I[i]+1;
  return List::create(Named("A")=A_out,Named("I")=I_out,Named("max_T")=max_T);
}
