#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::mat cpp_SoftImpute_thrs(arma::vec s, const double lambda, const int n, const int p){
  arma::mat output(n,p,fill::zeros);
  int m = 0;
  if (n < p){
    m = n;
  } else {
    m = p;
  }
  for (int i=0;i<m;i++){
    if (s(i)>lambda){
      output(i,i) = s(i)-lambda;
    } else {
      output(i,i) = 0.0;
    }
  }
  return(output);
}

// [[Rcpp::export]]
arma::cube cpp_SoftImpute(arma::mat& X, arma::mat& idmat, arma::vec& lambdas,
                          const double tol, const int maxiter, arma::mat& Zinit){
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

    // 4-2. REPEAT !
    while (incstop > tol){
      // 4-2-1. compute Ztmp before singular-value thresholding
      Ztmp(idx_obss) = X(idx_obss);
      Ztmp(idx_miss) = Zold(idx_miss);

      // 4-2-3. singular-value decomposition and thresholding
      svd(U, s, V, Ztmp);
      arma::mat thrS = cpp_SoftImpute_thrs(s, lambdak, N, P);
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
