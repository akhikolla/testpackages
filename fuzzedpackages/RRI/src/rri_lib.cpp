#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(cpp11)]]

//' Fast least squares
//'
//' This functions fits the regression y ~ X using Armadillo \code{solve} function.
//' @param y Vector of outcomes.
//' @param X Matrix of covariates (first column should be 1's)
//' @return \code{List} of regression output with elements \code{coef}, \code{stderr}.
// [[Rcpp::export]]
Rcpp::List fastLm(const arma::vec & y, const arma::mat & X) {

  int n = X.n_rows, k = X.n_cols;

  arma::colvec coef = arma::solve(X, y);
  arma::colvec resid = y - X*coef;

  double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
  arma::colvec stderrest =
    arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );

  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("stderr")       = stderrest);
}

//' Fast least squares
//'
//' Fast OLS as in \link{fastLm} but returns only the fitted coefficients.
//'
//' @param y Vector of outcomes.
//' @param X Matrix of covariates (first column should be 1's)
//' @return Vector of coefficients.
// [[Rcpp::export]]
arma::vec OLS_c(arma::vec y,
                arma::mat X) {

  arma::vec coef = arma::solve(X, y);
  return coef;
}

double OLS_c_1d(arma::vec y,
                arma::vec x) {

  arma::vec coef = arma::solve(x, y);
  return coef[0];
}


//' Fast least squares with linear constraint
//'
//' This functions fits the regression y ~ X under a linear constraint on the
//' model parameters. The constraint is  \code{Q} * beta = \code{c} where beta
//' are the regression model parameters, and \code{Q, c} are inputs.
//'
//' @param y Vector of outcomes.
//' @param X Matrix of covariates (first column should be 1's)
//' @param bhat Unconstrained OLS-fitted coefficients.
//' @param Q Matrix of linear constraints (k x p).
//' @param c Vector of constraint values (k x 1).
//'
//' @return Vector of fitted OLS coefficients under linear constraint.
//' @seealso Advanced Econometrics (Section 1.4, Takeshi Amemiya, 1985)
// [[Rcpp::export]]
arma::mat restricted_OLS_c(arma::vec y,
                           arma::mat X,
                           arma::vec bhat,
                           arma::mat Q,
                           double c) {
  arma::mat Sx = inv(trans(X) * X);

  return bhat - Sx * Q * inv(trans(Q) * Sx * Q) * (trans(Q) * bhat - c);
}

double Tn_c(arma::vec eps,
            arma::mat X,
            arma::vec lam) {

  arma::vec coef = fastLm(eps, X+0)[0]; // 0th index: coefficients, 1st index: stderr
  return dot(lam, coef);
}

arma::vec g_c(Rcpp::List cluster_eps,
              bool use_perm,
              bool use_sign) {
  // inputted partition will be a list of list, e.g. {{1, 2}, {3, 4, 5}, {6, 7}}
  // this is partitioned "er"
  int n_cluster = cluster_eps.size();
  // initilizing list with size n_cluster
  Rcpp::List permuted_cluster(n_cluster);

  // initializing vector that will be returned
  arma::vec out; //

  // IF there is only 1 cluster, just shuffle
  if (n_cluster == 1) {
    // Rprintf("entering here");
    arma::vec temp = cluster_eps[0];
    out = temp;
    if (use_perm) {
      out = arma::shuffle(temp);
    }
  }

  // IF more than 1 cluster, permute and randomly change sign on cluster level
  else {
    // 1. samples random signs (as many clusters)
    arma::vec rand_int = arma::randi<arma::vec>(n_cluster, arma::distr_param(0, 1)); // samples 0 and 1 randomly
    rand_int.replace(0, -1); // replace 0 with -1
    // rand_int.print();

    // 2. permutes the elements within the specified clusters
    // and multiply by randomly generated sign on cluster level
    for (int i=0; i<n_cluster; ++i) {
      // 1. grab i^th cluster
      arma::vec temp = cluster_eps[i];

      // 2. randomly shuffle i^th cluster, multiply random sign, and append to out ("er")
      if (use_perm && use_sign) {
        out = arma::join_cols(out, arma::shuffle(temp) * rand_int[i]);
      }
      else if (use_perm && !use_sign) {
        out = arma::join_cols(out, arma::shuffle(temp));
      }
      else if (!use_perm && use_sign) {
        out = arma::join_cols(out, temp * rand_int[i]);
      }
    }
    // print(cluster);
    // out.print();
  }
  return out;
}

//' Residual randomization test
//'
//' Implements the residual randomization test. The hypothesis tested is
//'
//'     H0: lam' beta = lam[1] * beta[1] + ... + lam[p] * beta[p] = lam0.
//'
//' @param y Vector of outcomes (n x 1).
//' @param X Matrix of covariates (n x p). First column should be 1's.
//' @param lam Vector of coefficients in linear H0 (p x 1).
//' @param lam0 Scalar value for linear H0.
//' @param cluster_eps_r A \code{List} with restricted residuals. See \link{get_clustered_eps}.
//' @param use_perm \code{Boolean} flag whether to use permutations within clusters.
//' @param use_sign \code{Boolean} flag whether to use sign flips across clusters.
//' @param num_R Integer of how many randomization values to calculate.
//' @return A \code{List} with the observed test statistic value (\code{tobs}), and the randomization values (\code{tvals})
// [[Rcpp::export]]
Rcpp::List r_test_c(arma::vec y,
                    arma::mat X,
                    arma::vec lam,
                    double lam0,
                    Rcpp::List cluster_eps_r,
                    bool use_perm,
                    bool use_sign,
                    int num_R) {
  if (y.size() != X.n_rows) {
    Rcpp::Rcout << "Error: length of y and nrow of X not matching, returning Inf" << std::endl;
    return std::numeric_limits<double>::infinity();
  }

  if (lam.size() != X.n_cols) {
    Rcpp::Rcout << "Error: length of lam and ncol of X not matching, returning Inf" << std::endl;
    return std::numeric_limits<double>::infinity();
  }

  arma::vec bhat = fastLm(y, X+0)[0]; // 0th index: coefficients, 1st index: stderr
  double tobs = dot(lam, bhat) - lam0;


  arma::vec tvals = arma::zeros(num_R + 1);
  for (int i=0; i<num_R; ++i) {
    arma::vec er_new = g_c(cluster_eps_r, use_perm, use_sign);
    // Rcpp::print(er_new_partitioned);
    tvals(i) = Tn_c(er_new, X, lam);
  }
  tvals(num_R) = tobs;
  // Return test statistic values.
  return Rcpp::List::create(Rcpp::Named("tobs") = tobs,
                            Rcpp::Named("tvals") = tvals);
}

