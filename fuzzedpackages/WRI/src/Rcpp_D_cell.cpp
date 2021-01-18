// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::export]]
arma::field<arma::mat> Rcpp_D_cell(const arma::mat & Xmat, const arma::mat & Q_res) {
        int n = Xmat.n_rows, m = Q_res.n_cols;
        arma::field<arma::mat> res(m, m);

        for (int i = 0; i < m; ++i) {
                for (int j = 0; j < i + 1; ++j) {
                        res(i, j) = Xmat.t() * arma::diagmat(Q_res.col(i) % Q_res.col(j)) * Xmat/n;
                        res(j, i) = res(i, j);

                }
        }

        return res;
}
