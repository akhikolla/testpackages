#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


void getMatrixBeta(const arma::mat & Z, 
                   const arma::vec & f, 
                   const int & k, 
                   const int & sel, 
                   arma::mat & betaOut, 
                   arma::mat & resOut, 
                   double & mseOut, 
                   double & critOut) {
  // This function finds the optimal beta, the residuals, mse and criterion corresponding to Z, f and k.
  // alpha is the last column of beta.
  // INPUT
  // Z: data matrix each ROW is a different time series
  // f: principal component
  // k: number of leads used
  // sel: LOO (1), AIC (2), BIC (3), BNG (4)
  // OUTPUT
  // betaOut: matrix of loadings and intercept corresponding to f. Last column are the
  // intercepts (alpha).
  // resOut: matrix of residuals
  // mseOut:  mean squared error (in N and m)
  // critOut: criterion used to evaluate the fit
  int m = Z.n_rows;
  int N = Z.n_cols;
  arma::mat Fmat = zeros(k + 1, N);
  arma::mat FF = zeros(k + 2, k + 2);
  arma::mat invFF = zeros(k + 2, k + 2);
  for ( int i = 0; i < N; i++){
    Fmat.col(i) = f.subvec(i, i + k);
  }
  Fmat.insert_rows(k + 1, ones(1, N));
  FF = Fmat * Fmat.t();
  double condition_FF = rcond(FF);
  if(condition_FF > 1e-10){
    invFF = inv_sympd(FF);
  } else{
    invFF = pinv(FF);
  }
  arma::mat Proj = Fmat.t() * invFF * Fmat;
  betaOut = Z * Fmat.t() * invFF.t();
  resOut = Z - betaOut * Fmat;
  mseOut = accu(pow(resOut, 2)) / N;
  if (sel == 1) {
    arma::vec weights = 1 / (1 - diagvec(Proj));
    critOut = accu( pow(resOut * diagmat(weights) , 2) ) / (N * m);
  } else if (sel == 2) {
    critOut = N * log(mseOut) + m * (k + 2) * 2;  
  } else if (sel == 3) {
    critOut = N * log(mseOut) + m * (k + 2) * log((double) N);  
  } else if (sel == 4) {
    arma::vec aux_min = zeros(2);
    aux_min(0) = N;
    aux_min(1) = m;
    double min_Nm = aux_min.min();
    critOut = min_Nm * log(mseOut) + (k + 1) * log(min_Nm);  
  }
  mseOut = mseOut / m;
}

arma::mat getMatrixC(const arma::subview_row<double> & rowZ,
                     const double & alpha,
                     const int & k) {
  // Constructs the matrix C correspoding to rowZ, alpha and k 
  int N = rowZ.n_elem;
  arma::mat C = zeros(N + k, k + 1);
  arma::vec liml = vec(2);
  liml[0] = 0;
  arma::vec limu = vec(2);
  limu[0] = k;
  for ( int t = 1; t <= N + k ; t++){
    liml[1] = t - N;
    limu[1] = t - 1;
    for (int q = 1; q <= k + 1; q++){
      if( (q >= max(liml) + 1) & (q <= min(limu) + 1)){
        C.at(t - 1, q - 1) = rowZ[t - q] - alpha;
      }
    }
  }
  return (C);
}

arma::mat getMatrixD(const arma::subview_row<double> & rowbeta,
                     const int & N,
                     const int & k) {
  // Constructs the matrix D corresponding to rowbeta, N and k 
  arma::mat beta_mat = zeros(N + k, N + k);
  arma::vec liml = vec(2);
  liml[1] = 1;
  arma::vec limu = vec(2);
  limu[1] = N;
  for ( int t = 1; t <= N + k ; t++){
    liml[0] = t - k;
    limu[0] = t;
    for (int r = max(liml); r <= min(limu); r++){
      for (int q = r; q <= r + k; q++){
        beta_mat.at(t - 1, q - 1) = beta_mat.at(t - 1, q - 1) + rowbeta[q - r] * rowbeta[t - r];
      }
    }
  }
  return (beta_mat);
}

arma::vec getF(const arma::mat & Z,
               const arma::mat & beta,
               const int & k) {
  // Get optimal f corresponding to Z, beta, alpha and k.
  int m = Z.n_rows;
  int N = Z.n_cols;
  arma::mat D = zeros(N + k, N + k);
  arma::vec fOut = zeros(N + k, 1);
  for (int j = 0 ; j < m ; j++){
    D = D + getMatrixD(beta(j, span(0, k)), N, k);
    fOut = fOut + getMatrixC(Z.row(j), beta(j, k + 1), k) * (beta(j, span(0, k))).t();
  }
  
  double condition_D = rcond(D);
  if (condition_D > 1e-10) {
    fOut = solve(D, fOut);
  } else {
    fOut = pinv(D) * fOut;
  }
  fOut = (fOut - mean(fOut)) / stddev(fOut);
  return(fOut);
}

// [[Rcpp::export]]
arma::mat getFitted(arma::vec & f_fin,
                    const arma::vec & f_ini,
                    const arma::mat & beta,
                    const arma::vec & alpha,
                    const int & k) {
  // Get fitted values associated with f and beta and alpha
  int N = f_fin.n_elem;
  if (k > 0){
    f_fin.insert_rows(0, f_ini);
  }
  arma::mat Fmat = zeros(k + 1, N);
  for ( int i = 0; i < N; i++){
    Fmat.col(i) = f_fin.subvec(i, i + k);
  }
  Fmat.insert_rows(k + 1, ones(1, N));
  arma::mat betalpha = fliplr(beta);
  betalpha.insert_cols(k + 1, alpha);
  arma::mat fit = Fmat.t() * betalpha.t();
  return(fit);
}

// [[Rcpp::export]]
arma::vec getFini(const arma::mat & Z,
                  const int & k) {
  // Get initial estimator: ordinary principal component with k leads
  int N = Z.n_cols;
  arma::vec f_iniOut = zeros(N + k, 1);
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::mat Z_trans = Z.t();
  arma::rowvec mean_Zt = mean(Z_trans);
  Z_trans.each_row() -= mean_Zt; 
  svd_econ(U, s, V, Z_trans, "right");
  f_iniOut.rows(0, N - 1) = Z_trans * V.col(0);
  if (k != 0) {
    f_iniOut.rows(N, N + k - 1) = zeros(k, 1) + f_iniOut(N - 1);
  }
  f_iniOut = (f_iniOut - mean(f_iniOut)) / stddev(f_iniOut);
  return(f_iniOut);
}


// [[Rcpp::export]]
List gdpc_priv(const arma::mat & Z,
               const int & k,
               const arma::vec & f_ini,
               const bool & passf_ini,
               const double & tol,
               const int & niter_max,
               const int & sel) {
// This function computes a single GDPC with a given number of leads.
// INPUT
// Z: data matrix each ROW is a different time series
// k: number of leads used
// f_ini: starting point for the iterations, optional.
// passf_ini: logical, indicates whether f_ini is being passed or not
// tol: relative precision, stopping criterion
// niter_max: maximum number of iterations
// sel: LOO (1), AIC (2), BIC (3), BNG (4)
// OUTPUT
// f: principal component
// beta: matrix of loadings and intercept corresponding to f. Last column are the
// intercepts (alpha).
// mse:  mean squared error (in N and m)
// crit: criterion used to evaluate the fit
// res: matrix of residuals
// conv: logical. Did the iterations converge?
  int m = Z.n_rows;
  int N = Z.n_cols;
  double criter = tol + 1;  //Stopping criterion for the iterations
  int niter = 0;
  arma::mat beta = zeros(m, k + 2);
  arma::mat res = zeros(m, N);
  arma::vec f = zeros(N + k, 1);
  double mse = 0;
  double crit = 0;
  bool conv = false;
  if (!passf_ini) {
    f = getFini(Z, k);  
  } else {
    f = f_ini;
  }
  getMatrixBeta(Z, f, k, sel, beta, res, mse, crit);
  double mse_ini = mse;
  while (niter < niter_max and criter > tol) {
    niter = niter + 1;
    f = getF(Z, beta, k);
    getMatrixBeta(Z, f, k, sel, beta, res, mse, crit);
    criter = 1 - mse / mse_ini;
    mse_ini = mse;
    if (niter % 10 == 0){
      Rcpp::checkUserInterrupt();
    }
  }
  if (niter < niter_max) { 
    conv = true;
  }
  List ret;
  ret["f"] = f;
  ret["k"] = k;
  ret["beta"] = beta;
  ret["mse"] = mse;
  ret["crit"] = crit;
  ret["res"] = res;
  ret["conv"] = conv;
  ret["niter"] = niter;
  return(ret);
}
