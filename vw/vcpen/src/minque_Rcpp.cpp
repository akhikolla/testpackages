#include <RcppArmadillo.h>

using namespace Rcpp; using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
List minque_Rcpp(arma::vec& y, arma::mat& X, const List& Kerns, arma::vec vc){
  // vc vec of init values of VCs
  
  int n_vc = Kerns.size();
  int n_subj = y.size();
  int n_beta  = X.n_cols;

  mat v(n_subj, n_subj, fill::zeros);
  mat vinv = v;
  mat Px = v;
  
  for(int i=0; i < n_vc; i++){
    v  = v + vc(i) * as<mat>(Kerns[i]);
  }
  vinv = inv(v);
  mat tx = trans(X) * vinv;

  Px  = inv(tx * X) * tx;
 
  mat beta = Px * y;
  mat resid = y - X * beta;
  
  mat Pv = X * Px;
  mat ident(n_subj, n_subj, fill::zeros);
  ident.eye();
  mat Qv = ident - Pv;
  mat e = Qv * y;

  vec uvec(n_vc, fill::zeros);
  mat e_vinv = trans(e) * vinv;
  mat mat_scalar(1,1);
  for(int i=0; i < n_vc; i++){
    mat_scalar = e_vinv * as<mat>(Kerns[i]) * trans(e_vinv);
    uvec(i) = mat_scalar(0,0);
  }
 
  mat smat(n_vc, n_vc, fill::zeros);
  mat Qvinv = trans(Qv) * vinv;
  
  for(int i=0; i < n_vc; i++){
    for(int j=i; j < n_vc; j++){
      smat(i,j) = trace(Qvinv * as<mat>(Kerns[i]) * trans(Qvinv) * as<mat>(Kerns[j]));
      smat(j,i) =  smat(i,j);
    }
  }
  

  mat minque_vc =  solve(smat, uvec);

  return Rcpp::List::create(Rcpp::Named("vc") = minque_vc,
			    Rcpp::Named("beta") = beta,
			    Rcpp::Named("residuals")  = resid);
}
