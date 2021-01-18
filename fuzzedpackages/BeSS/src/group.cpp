// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;
using namespace std;
vector<int> uniqueR(vector<int> x){
  vector<int>::iterator it;
  sort (x.begin(),x.end());
  it = unique(x.begin(), x.end());
  x.resize(distance(x.begin(),it));
  return x;
}
// [[Rcpp::export]]
Eigen::MatrixXd EigenR(Eigen::MatrixXd X){
  Eigen::EigenSolver<Eigen::MatrixXd> eig(X);
  Eigen::MatrixXd value = eig.pseudoEigenvalueMatrix();
  Eigen::MatrixXd vectors = eig.pseudoEigenvectors();
  Eigen::MatrixXd y = vectors*value.cwiseSqrt()*vectors.adjoint();
  return y;
}
// [[Rcpp::export]]
List gbess_lm(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd G, Eigen::VectorXd index, List PhiG, List invPhiG, int T0, int max_steps, Eigen::VectorXd beta0, int n, int p, int N){
  double max_T=0.0;
  Eigen::VectorXi B0 = Eigen::VectorXi::Zero(p);
  Eigen::VectorXd d0 = (X.transpose()*(y-X*beta0)) /double(n);
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd betabar = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd dbar = Eigen::VectorXd::Zero(p);
  vector<int>A(T0);
  vector<int>A0(T0);
  Eigen::VectorXi A_out = Eigen::VectorXi::Zero(T0);
  for(int i=0;i<T0;i++){
    A0[i] = 1;
  }
  vector<int>gr(T0);
  int mark = 0;
  int gr_size = 0;
  double mse = 0.0;
  for(int l=1;l<=max_steps;l++) {
    for(int i=0;i<N-1;i++){
      Eigen::MatrixXd phiG = PhiG[i];
      Eigen::MatrixXd invphiG = invPhiG[i];
      betabar.segment(index(i)-1,index(i+1)-index(i)) = phiG*beta0.segment(index(i)-1,index(i+1)-index(i));
      dbar.segment(index(i)-1,index(i+1)-index(i)) = invphiG*d0.segment(index(i)-1,index(i+1)-index(i));
    }
    Eigen::MatrixXd phiG = PhiG[N-1];
    Eigen::MatrixXd invphiG = invPhiG[N-1];
    betabar.segment(index(N-1)-1,p-index(N-1)+1) = phiG*beta0.segment(index(N-1)-1,p-index(N-1)+1);
    dbar.segment(index(N-1)-1,p-index(N-1)+1) = invphiG*d0.segment(index(N-1)-1,p-index(N-1)+1);
    Eigen::VectorXd bd = (betabar+dbar).array().abs();
    for(int k=0;k<T0;k++) {
      max_T=bd.maxCoeff(&A[k]);
      //cout<<bd(A[k])<<" ";
      bd(A[k])=0;
    }
    sort (A.begin(),A.end());
    if(A==A0) break;
    for(int i=0;i<T0;i++){
      gr[i] = G(A[i]);
    }
    vector<int>gr_unique = uniqueR(gr);
    mark = 0;
    gr_size = int(gr_unique.size());
    for(int i=0;i<gr_size-1;i++){
      B0.segment(mark, index(gr_unique[i])-index(gr_unique[i]-1)) = Eigen::VectorXi::LinSpaced(index(gr_unique[i])-index(gr_unique[i]-1), index(gr_unique[i]-1)-1, index(gr_unique[i])-2);
      mark = mark+index(gr_unique[i])-index(gr_unique[i]-1);
    }
    if(G(p-1)==gr_unique[gr_size-1]){
      B0.segment(mark, p-index(gr_unique[gr_size-1]-1)+1) = Eigen::VectorXi::LinSpaced(p-index(gr_unique[gr_size-1]-1)+1, index(gr_unique[gr_size-1]-1)-1, p-1);
      mark = mark + p-index(gr_unique[gr_size-1]-1)+1;
    }else{
      B0.segment(mark, index(gr_unique[gr_size-1])-index(gr_unique[gr_size-1]-1)) = Eigen::VectorXi::LinSpaced(index(gr_unique[gr_size-1])-index(gr_unique[gr_size-1]-1), index(gr_unique[gr_size-1]-1)-1, index(gr_unique[gr_size-1])-2);
      mark = mark + index(gr_unique[gr_size-1])-index(gr_unique[gr_size-1]-1);
    }

    Eigen::VectorXi B = B0.head(mark);
    Eigen::MatrixXd X_B(n, B.size());

    for(int i=0;i<B.size();i++){
      X_B.col(i) = X.col(B(i));
    }
    Eigen::VectorXd beta_B=X_B.colPivHouseholderQr().solve(y);
    Eigen::VectorXd d = (X.transpose()*(y-X_B*beta_B));
    for(int i=0;i<B.size();i++){
      beta(B(i)) = beta_B(i);
      d(B(i)) = 0;
    }
    beta0 = beta;
    d0 = d;
    A0 = A;
  }
  Eigen::VectorXi B = B0.head(mark);
  mse=(y-X*beta).array().square().sum()/double(n);
  for(int i=0;i<T0;i++)
    A_out(i) = A[i];
  return List::create(Named("beta")=beta,Named("max_T")=max_T,Named("mse")=mse, Named("A")=A_out, Named("B")=B, Named("gr_size")=gr_size);
}
// [[Rcpp::export]]
List gget_A(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXi G, Eigen::VectorXi index, int T0, Eigen::VectorXd beta0, double coef0, int n, int p, int N, Eigen::VectorXd weights, Eigen::VectorXi B00){
  double max_T=0.0;
  Eigen::VectorXi A_out = Eigen::VectorXi::Zero(T0);
  Eigen::VectorXd betabar = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd dbar = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd one = Eigen::VectorXd::Ones(n);
  Eigen::VectorXd bd = Eigen::VectorXd::Zero(p);
  Eigen::VectorXi B0 = Eigen::VectorXi::Zero(p);
  vector<int> A(T0);
  vector<int> gr(T0);
  Eigen::VectorXd coef(n);
  for(int i=0;i<=n-1;i++) {
    coef(i)=coef0;
  }
  int mark=0;
  Eigen::VectorXd eta = X*beta0+coef;
  for(int i=0;i<=n-1;i++) {
    if(eta(i)<-25.0) eta(i) = -25.0;
    if(eta(i)>25.0) eta(i) = 25.0;
  }
  Eigen::VectorXd exppr = eta.array().exp();
  Eigen::VectorXd pr = exppr.array()/(exppr+one).array();
  Eigen::VectorXd g = weights.array()*(y-pr).array();
  Eigen::VectorXd h = weights.array()*pr.array()*(one-pr).array();
  Eigen::VectorXd d0 = X.adjoint()*g;
  if(B00.size()<p){
    for(int i=0;i<B00.size();i++){
      d0(B00(i)-1)=0.0;
    }
  }
  for(int i=0;i<N-1;i++){
    Eigen::MatrixXd XG =X.middleCols(index(i)-1, index(i+1)-index(i));
    Eigen::MatrixXd XGbar = XG.adjoint()*h.asDiagonal()*XG;
    Eigen::MatrixXd phiG = -EigenR(XGbar);
    Eigen::MatrixXd invphiG = phiG.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(index(i+1)-index(i),index(i+1)-index(i)));
    betabar.segment(index(i)-1,index(i+1)-index(i)) = phiG*beta0.segment(index(i)-1,index(i+1)-index(i));
    dbar.segment(index(i)-1,index(i+1)-index(i)) = invphiG*d0.segment(index(i)-1,index(i+1)-index(i));
  }
  Eigen::MatrixXd XG =X.middleCols(index(N-1)-1,p-index(N-1)+1);
  Eigen::MatrixXd XGbar = XG.adjoint()*h.asDiagonal()*XG;
  Eigen::MatrixXd phiG = -EigenR(XGbar);
  Eigen::MatrixXd invphiG = phiG.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(p-index(N-1)+1,p-index(N-1)+1));
  betabar.segment(index(N-1)-1,p-index(N-1)+1) = phiG*beta0.segment(index(N-1)-1,p-index(N-1)+1);
  dbar.segment(index(N-1)-1,p-index(N-1)+1) = invphiG*d0.segment(index(N-1)-1,p-index(N-1)+1);
  bd = (betabar+dbar).array().abs();
  for(int k=0;k<T0;k++) {
    max_T=bd.maxCoeff(&A[k]);
    //cout<<bd(A[k])<<" ";
    bd(A[k])=0.0;
  }
  sort(A.begin(),A.end());
  for(int i=0;i<T0;i++){
    gr[i] = G(A[i]);
  }
  vector<int>gr_unique = uniqueR(gr);
  mark = 0;
  int gr_size = int(gr_unique.size());
  for(int i=0;i<gr_size-1;i++){
    B0.segment(mark, index(gr_unique[i])-index(gr_unique[i]-1)) = Eigen::VectorXi::LinSpaced(index(gr_unique[i])-index(gr_unique[i]-1), index(gr_unique[i]-1)-1, index(gr_unique[i])-2);
    mark = mark+index(gr_unique[i])-index(gr_unique[i]-1);
  }
  if(G(p-1)==gr_unique[gr_size-1]){
    B0.segment(mark, p-index(gr_unique[gr_size-1]-1)+1) = Eigen::VectorXi::LinSpaced(p-index(gr_unique[gr_size-1]-1)+1, index(gr_unique[gr_size-1]-1)-1, p-1);
    mark = mark+p-index(gr_unique[gr_size-1]-1)+1;
  }else{
    B0.segment(mark, index(gr_unique[gr_size-1])-index(gr_unique[gr_size-1]-1)) = Eigen::VectorXi::LinSpaced(index(gr_unique[gr_size-1])-index(gr_unique[gr_size-1]-1), index(gr_unique[gr_size-1]-1)-1, index(gr_unique[gr_size-1])-2);
    mark = mark+index(gr_unique[gr_size-1])-index(gr_unique[gr_size-1]-1);
  }
  Eigen::VectorXi B = B0.head(mark);
  for(int i=0;i<T0;i++)
    A_out(i) = A[i];
  return List::create(Named("p")=pr,Named("A")=A_out,Named("B")=B,Named("max_T")=max_T,Named("gr_size")=gr_size);
}
// [[Rcpp::export]]
List ggetcox_A(Eigen::MatrixXd X, Eigen::VectorXi G, Eigen::VectorXi index, int T0, Eigen::VectorXd beta0, int n, int p, int N, Eigen::VectorXd status, Eigen::VectorXd weights, Eigen::VectorXi B00){
  double max_T=0.0;
  Eigen::VectorXi A_out = Eigen::VectorXi::Zero(T0);
  Eigen::VectorXd betabar = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd dbar = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd bd = Eigen::VectorXd::Zero(p);
  Eigen::VectorXi B0 = Eigen::VectorXi::Zero(p);
  Eigen::MatrixXd h(n, n);
  Eigen::VectorXd cum_theta(n);
  Eigen::VectorXd cum_theta2(n);
  Eigen::VectorXd cum_theta3(n);
  vector<int> A(T0);
  vector<int> gr(T0);
  int mark=0;
  Eigen::VectorXd theta = X*beta0;
  for(int i=0;i<=n-1;i++) {
    if(theta(i)<-25.0) theta(i) = -25.0;
    if(theta(i)>25.0) theta(i) = 25.0;
  }
  theta = weights.array()*theta.array().exp();
  cum_theta(n-1) = theta(n-1);
  for(int k=n-2;k>=0;k--) {
    cum_theta(k) = cum_theta(k+1)+theta(k);
  }
  cum_theta2(0) = (status(0)*weights(0))/cum_theta(0);
  for(int k=1;k<=n-1;k++) {
    cum_theta2(k) = (status(k)*weights(k))/cum_theta(k)+cum_theta2(k-1);
  }
  cum_theta3(0) = (status(0)*weights(0))/pow(cum_theta(0),2);
  for(int k=1;k<=n-1;k++) {
    cum_theta3(k) = (status(k)*weights(k))/pow(cum_theta(k),2)+cum_theta3(k-1);
  }
  h = -cum_theta3.replicate(1, n);
  h = h.cwiseProduct(theta.replicate(1, n));
  h = h.cwiseProduct(theta.replicate(1, n).adjoint());
  for(int i=0;i<n;i++)
    for(int j=i+1;j<n;j++)
      h(j, i) = h(i, j);
  Eigen::VectorXd g = cum_theta2.cwiseProduct(theta) - weights.cwiseProduct(status);
  h.diagonal() = cum_theta2.cwiseProduct(theta) + h.diagonal();
  Eigen::VectorXd d0 = X.adjoint()*g;
  if(B00.size()<p){
    for(int i=0;i<B00.size();i++){
      d0(B00(i)-1)=0.0;
    }
  }
  for(int i=0;i<N-1;i++){
    Eigen::MatrixXd XG =X.middleCols(index(i)-1, index(i+1)-index(i));
    Eigen::MatrixXd XGbar = XG.adjoint()*h*XG;
    //cout<<XGbar<<" ";
    Eigen::MatrixXd phiG = -EigenR(XGbar);
    Eigen::MatrixXd invphiG = phiG.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(index(i+1)-index(i),index(i+1)-index(i)));
    betabar.segment(index(i)-1,index(i+1)-index(i)) = phiG*beta0.segment(index(i)-1,index(i+1)-index(i));
    dbar.segment(index(i)-1,index(i+1)-index(i)) = invphiG*d0.segment(index(i)-1,index(i+1)-index(i));
    //cout<<XGbar<<endl;
  }
  Eigen::MatrixXd XG =X.middleCols(index(N-1)-1,p-index(N-1)+1);
  Eigen::MatrixXd XGbar = XG.adjoint()*h*XG;
  Eigen::MatrixXd phiG = -EigenR(XGbar);
  Eigen::MatrixXd invphiG = phiG.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(p-index(N-1)+1,p-index(N-1)+1));
  betabar.segment(index(N-1)-1,p-index(N-1)+1) = phiG*beta0.segment(index(N-1)-1,p-index(N-1)+1);
  dbar.segment(index(N-1)-1,p-index(N-1)+1) = invphiG*d0.segment(index(N-1)-1,p-index(N-1)+1);
  bd = (betabar+dbar).array().abs();
  for(int k=0;k<T0;k++) {
    max_T=bd.maxCoeff(&A[k]);
    //cout<<bd(A[k])<<" ";
    bd(A[k])=0.0;
  }
  sort(A.begin(),A.end());
  for(int i=0;i<T0;i++){
    gr[i] = G(A[i]);
  }
  vector<int>gr_unique = uniqueR(gr);
  mark = 0;
  int gr_size = int(gr_unique.size());
  for(int i=0;i<gr_size-1;i++){
    B0.segment(mark, index(gr_unique[i])-index(gr_unique[i]-1)) = Eigen::VectorXi::LinSpaced(index(gr_unique[i])-index(gr_unique[i]-1), index(gr_unique[i]-1)-1, index(gr_unique[i])-2);
    mark = mark+index(gr_unique[i])-index(gr_unique[i]-1);
  }
  if(G(p-1)==gr_unique[gr_size-1]){
    B0.segment(mark, p-index(gr_unique[gr_size-1]-1)+1) = Eigen::VectorXi::LinSpaced(p-index(gr_unique[gr_size-1]-1)+1, index(gr_unique[gr_size-1]-1)-1, p-1);
    mark = mark+p-index(gr_unique[gr_size-1]-1)+1;
  }else{
    B0.segment(mark, index(gr_unique[gr_size-1])-index(gr_unique[gr_size-1]-1)) = Eigen::VectorXi::LinSpaced(index(gr_unique[gr_size-1])-index(gr_unique[gr_size-1]-1), index(gr_unique[gr_size-1]-1)-1, index(gr_unique[gr_size-1])-2);
    mark = mark+index(gr_unique[gr_size-1])-index(gr_unique[gr_size-1]-1);
  }
  Eigen::VectorXi B = B0.head(mark);
  for(int i=0;i<T0;i++)
    A_out(i) = A[i];
  return List::create(Named("A")=A_out,Named("B")=B,Named("max_T")=max_T,Named("gr_size")=gr_size);
}
