// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

arma::vec scad_group(arma::vec beta, double lam1){
  
  double scada = 3.7;
  arma::vec abs_beta=abs(beta);
  arma::uvec pos = find(abs_beta<lam1);
  int dim = abs_beta.size();
  arma::vec s(dim);
  s.fill(0);
  s.elem(pos) = lam1*abs_beta.elem(pos);
  
  arma::uvec pos2 = find(abs_beta>lam1);
  arma::uvec pos3 = find(abs_beta<=lam1*3.7);
  arma::uvec pos4 = intersect(pos2,pos3);
  
  arma::vec beta_pos4 = abs_beta(pos4);
  arma::vec beta_pos4_2 = pow(beta_pos4,2);
  arma::vec beta_pos4_3 = 2*scada*lam1*beta_pos4;
  arma::vec s_2 = beta_pos4_2-beta_pos4_3;
  double lam12 = pow(lam1,2);
  double dominator = 2*(scada-1);
  s.elem(pos4) = -(s_2+lam12)/dominator;
  
  double alam = scada*lam1;
  arma::uvec pos5 = find(abs_beta>=alam);
  double tmp = (scada+1)*pow(lam1,2)/2;
  s.elem(pos5).fill(tmp);
  
  
  return s;
}

double scad(double a, double lambda){
  double scada = 3.7;
  double abs_a = abs(a);
  double out = 0;
  if(abs_a<lambda){out = lambda;}
  else if((abs_a>lambda)&(abs_a<=scada*lambda)){
    out = (scada*lambda-abs_a)/(scada-1);
  }
  return out;
}

double L2norm(arma::rowvec a){
  arma::rowvec a2 = pow(a,2);
  double a2sum = sum(a2);
  return sqrt(a2sum);
}

arma::vec rowsums(arma::mat x) {
  
  int nRows = x.n_rows;
  arma::vec out(nRows);
  
  for( int i=0; i < nRows; i++ ) {
    arma::rowvec tmp = x.row(i);
    out(i) = L2norm(tmp);
  }
  
  return out;
  
}

// [[Rcpp::export]]
arma::mat D_Gradient(Rcpp::List xGroup,Rcpp::List zGroup,Rcpp::List LGroup,Rcpp::List yGroup,arma::vec b,int N,int q) {
  
  
  arma::mat out(q,q,fill::zeros); 
  
  for(int i=0;i<N;i++){
    arma::mat z = zGroup[i];
    arma::mat V = LGroup[i];
    arma::vec y = yGroup[i];
    arma::mat x = xGroup[i];
    
    arma::mat ztvztmp = z.t()*V;
    arma::mat ztvz = ztvztmp*z;
    
    arma::vec ri = y-x*b;
    
    arma::mat zv = z.t()*V;
    arma::mat zvri = zv*ri;
    
    arma::mat vv = zvri*zvri.t();
    
    out = out-ztvz+vv;
  }
  
  return out;
  
}

// [[Rcpp::export]]
arma::mat D_HessianMatrix(Rcpp::List xGroup,Rcpp::List zGroup,Rcpp::List LGroup,Rcpp::List yGroup,arma::vec b,int N,int q){
  
  int q2 = pow(q,2);
  arma::mat out(q2,q2,fill::zeros); 
  for(int i=0;i<N;i++){
    arma::mat z = zGroup[i];
    arma::mat V = LGroup[i];
    arma::vec y = yGroup[i];
    arma::mat x = xGroup[i];
    
    arma::mat ztvztmp = z.t()*V;
    arma::mat ztvz = ztvztmp*z;
    
    arma::vec ri = y-x*b;
    
    arma::mat zv = z.t()*V;
    arma::mat zvri = zv*ri;
    
    arma::mat ztvz2 = kron(ztvz,ztvz);
    
    arma::mat vv = zvri*zvri.t();
    arma::mat ztvzvv = kron(ztvz,vv);
    
    
    out = out-ztvz2+ztvzvv+ztvzvv.t();
  }
  
  return out;
  
}

// [[Rcpp::export]]
double ObjFunction(Rcpp::List xGroup, Rcpp::List yGroup, Rcpp::List LGroup, arma::vec b_nonpen, arma::mat L_nonpen, double lambda1,double lambda2,
                   Rcpp::String penalty_b, Rcpp::String penalty_L, double ll1){
  
  arma::uvec activeSet = find(b_nonpen!=0);
  arma::vec b_active = b_nonpen(activeSet);
  
  double rtvr=0;
  double ll2=0;
  
  for (int i = 0; i < xGroup.length(); ++i) {
    arma::mat x = xGroup[i];
    arma::vec y = yGroup[i];
    arma::mat V = LGroup[i];
    
    arma::mat x_active = x.cols(activeSet);
    arma::vec ri = y-x_active*b_active;
    
    arma::mat vitemp = ri.t()*V;
    double vi = as_scalar(vitemp*ri)*1/2;
    
    rtvr = rtvr+vi;
    
    double det_v = det(V);
    if(det_v==0){
      det_v =det_v+0.00000000000000000000000000000001;
    }
    double abs_det_v = std::abs(det_v);
    double logdetv = log(abs_det_v);
    ll2 = ll2-logdetv*1/2;
    
  }
  
  
  double pen_b=0;
  
  if(penalty_b=="lasso"){arma::vec abs_b = abs(b_nonpen); double sum_b = sum(abs_b); pen_b = lambda1*sum_b;} 
  else if(penalty_b=="scad") {arma::vec pen_b_presum = scad_group(b_nonpen,lambda1); pen_b = sum(pen_b_presum);}
  
  
  
  double pen_L=0;
  
  arma::vec l2norm = rowsums(L_nonpen);
  
  
  if(penalty_L=="lasso"){pen_L = lambda2*sum(l2norm);}
  else if(penalty_L=="scad"){arma::vec pen_L_presum = scad_group(l2norm,lambda2); pen_L = sum(pen_L_presum);}
  
  
  return rtvr+ll1+ll2+pen_b+pen_L;
}

// [[Rcpp::export]]
double ObjFunction_L(Rcpp::List xGroup, Rcpp::List yGroup, Rcpp::List zGroup, Rcpp::List zIdGrp,
                     arma::vec b_nonpen, arma::mat L, arma::uvec nonpen, double sigma, double lambda, Rcpp::String penalty, double ll1){
  
  arma::uvec activeSet = find(b_nonpen!=0);
  arma::vec b_active = b_nonpen(activeSet);
  
  double rtvr = 0;
  double ll2 = 0;
  
  for (int i = 0; i < xGroup.length(); ++i) {
    arma::mat x = xGroup[i];
    arma::vec y = yGroup[i];
    arma::mat z = zGroup[i];
    arma::mat ZID = zIdGrp[i];
    
    
    arma::mat D = L*L.t();
    arma::mat Vtmp = z*D;
    arma::mat Vtmp2 = Vtmp*z.t();
    arma::mat Sig = pow(sigma,2)*ZID;
    arma::mat V = Sig+Vtmp2;
    
    double V_det = det(V);
    
    arma::mat Vinv;
    if(V_det!=0){
      Vinv = inv(V);
    }else{
      Vinv = pinv(V);
    }
    
    
    arma::mat x_active = x.cols(activeSet);
    arma::vec ri = y-x_active*b_active;
    
    arma::mat vitemp = ri.t()*Vinv;
    
    
    double vi = as_scalar(vitemp*ri)*1/2;
    rtvr = rtvr+vi;
    
    double det_v = det(Vinv);
    if(det_v==0){
      det_v =det_v+0.00000000000000000000000000000001;
    }
    double abs_det_v = std::abs(det_v);
    double logdetv = log(abs_det_v);
    ll2 = ll2-logdetv*1/2;
    
  }
  
  
  
  double pen_L=0;
  
  arma::mat L_nonpen = L;
  L_nonpen.shed_rows(nonpen);
  
  arma::vec l2norm = rowsums(L_nonpen);
  
  
  if(penalty=="lasso"){pen_L = lambda*sum(l2norm);}
  else if(penalty=="scad"){arma::vec pen_L_presum = scad_group(l2norm,lambda); pen_L = sum(pen_L_presum);}
  
  
  return rtvr+ll1+ll2+pen_L;
  
}


// [[Rcpp::export]]
Rcpp::List ArmijoRule_L(Rcpp::List xGroup, Rcpp::List yGroup, Rcpp::List zGroup, arma::mat L, int l, int k, double grad, double hessian, 
                        arma::vec b, double sigma, Rcpp::List zIdGrp, bool linNonpen, arma::uvec nonpen, double lambda, Rcpp::String penalty, 
                        double ll1, double gamma, int maxArmijo, double a_init, double delta, double rho, int converged){
  arma::mat L_new = L;
  double Llk = L(l,k);
  double dk;
  
  arma::rowvec tmp = L.row(l);
  double l2norm = L2norm(tmp);
  if(linNonpen){
    dk = -grad/hessian;
  }else{
    
    if(penalty=="lasso"){
      dk = (-grad-lambda/l2norm*Llk)/(hessian+lambda/l2norm);
    }
    else if(penalty=="scad"){
      double group_scad = scad(l2norm,lambda);
      dk = (-grad-lambda/l2norm*Llk)/(hessian+lambda/group_scad);
    }
  }
  
  double fct = ObjFunction_L(xGroup, yGroup, zGroup, zIdGrp,
                             b, L, nonpen, sigma, lambda, penalty, ll1);
  
  
  
  double deltak =0;
  
  if(dk!=0){
    if (linNonpen) {
      deltak = dk*grad + gamma*pow(dk,2)*hessian;
    }else{
      arma::mat L_tmp = L;
      L_tmp(l,k) = L(l,k)+dk;
      arma::rowvec tmp = L_tmp.row(l);
      double l2norm_tmp = L2norm(tmp);
      deltak = dk*grad+gamma*pow(dk,2)*hessian+lambda*(l2norm_tmp-l2norm);
    }
    
    
    double fctOld = ObjFunction_L(xGroup, yGroup, zGroup, zIdGrp,
                                  b, L, nonpen, sigma, lambda, penalty, ll1);
    
    for(int j=0; j<=maxArmijo; j++){
      
      L_new(l,k) = Llk + a_init*pow(delta,j)*dk;
      
      double fctNew = ObjFunction_L(xGroup, yGroup, zGroup, zIdGrp,
                                    b, L_new, nonpen, sigma, lambda, penalty, ll1);
      
      double addDelta = a_init*pow(delta,j)*rho*deltak;
      
      if (fctNew <= fctOld + addDelta)
      {
        L(l,k) = Llk+a_init*pow(delta,j)*dk;
        fct = fctNew;
        break;
      }
      if (j==maxArmijo)
      {
        converged = converged + 2;
        fct = fctOld;
      }
      
    }
  }
  List armijo = List::create(Named("L")=L, Named("fct")=fct, Named("converged")=converged, Named("dk")=dk, Named("deltak")=deltak);
  
  return armijo;
}
