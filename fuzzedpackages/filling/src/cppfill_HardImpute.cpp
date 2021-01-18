#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::mat cpp_HardImpute_thrs(arma::vec s, const double lambda, const int n, const int p, const int rk){
  arma::mat output(n,p,fill::zeros);
  int m = 0;
  if (n < p){
    m = n;
  } else {
    m = p;
  }
  for (int i=0;i<rk;i++){
    output(i,i) = s(i);
  }
  return(output);
}

// [[Rcpp::export]]
arma::cube cpp_HardImpute(arma::mat& X, arma::mat& idmat, arma::vec& lambdas,
                          const double tol, const int maxiter, arma::cube& SoftRes, const int rk){
  // 1. settings
  const int N = X.n_rows;
  const int P = X.n_cols;
  const int K = lambdas.n_elem;

  // 2. set empty ones
  arma::cube result(N,P,K,fill::zeros); // where everything will be recorded
  //arma::mat  Zold = Zinit;
  arma::mat  Zold(N,P,fill::zeros);
  arma::mat  Ztmp(N,P,fill::zeros);
  arma::mat  Znew(N,P,fill::zeros);

  // 3. indexes
  arma::uvec idx_miss = find(idmat > 0.5);
  arma::uvec idx_obss = find(idmat <=0.5);

  // 4. main iteration
  double lambdak = 0.0;
  double incstop = 0.0;
  int iterk = 0;
  arma::mat U;
  arma::vec s;
  arma::mat V;

  for (int k=0;k<K;k++){
    // 4-1. current lambda value and incremental threshold
    lambdak = lambdas(k);
    incstop = 10000.0;
    iterk   = 0;
    Zold    = SoftRes.slice(k);

    // 4-2. REPEAT !
    while (incstop > tol){
      // 4-2-1. compute Ztmp before singular-value thresholding
      Ztmp(idx_obss) = X(idx_obss);
      Ztmp(idx_miss) = Zold(idx_miss);

      // 4-2-3. singular-value decomposition and thresholding
      svd(U, s, V, Ztmp);
      arma::mat thrS = cpp_HardImpute_thrs(s, lambdak, N, P, rk);
      Znew = U*thrS*V.t();

      // 4-2-4. compute incstop
      incstop = std::pow(arma::norm(Zold-Znew,"fro")/arma::norm(Znew,"fro"),2);

      // 4-2-5. assign
      Zold = Znew;

      // 4-2-6. iteration update
      iterk = iterk+1;
      if (iterk >= maxiter){
        break;
      }
    }

    // 4-3. assign
    result.slice(k) = Zold;
  }

  // 5. return resutls
  return(result);
}
