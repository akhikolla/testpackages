// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::export]]
Rcpp::List Rcpp_D_4cell(const arma::mat & Xmat, const arma::mat & Q_res, const arma::mat & q_res) {
        int n = Xmat.n_rows, m = Q_res.n_cols;
        arma::field<arma::mat> D_00(m, m);
        arma::field<arma::mat> D_01(m, m);
        arma::field<arma::mat> D_10(m, m);
        arma::field<arma::mat> D_11(m, m);

        for (int i = 0; i < m; ++i) {
                for (int j = 0; j < i + 1; ++j) {
                        D_00(i, j) = Xmat.t() * arma::diagmat(Q_res.col(i) % Q_res.col(j)) * Xmat/n;
                        D_00(j, i) = D_00(i, j);

                        D_11(i, j) = Xmat.t() * arma::diagmat(q_res.col(i) % q_res.col(j)) * Xmat/n;
                        D_11(j, i) = D_11(i, j);
                }
        }

        for (int i = 0; i < m; ++i) {
                for (int j = 0; j < m; ++j) {
                        D_10(i, j) = Xmat.t() * arma::diagmat(Q_res.col(i) % q_res.col(j)) * Xmat/n;
                        D_01(j, i) =  D_10(i, j);
                }
        }

        return Rcpp::List::create(Rcpp::Named("D_00") = D_00,
                                  Rcpp::Named("D_01") = D_01,
                                  Rcpp::Named("D_10") = D_10,
                                  Rcpp::Named("D_11") = D_11);
}


