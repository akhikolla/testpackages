// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::export]]

arma::mat kernel_partialF(const arma::mat & Xc, const arma::mat & Q_res, const arma::mat & left_mat, int r) {
        int n = Xc.n_rows, m = Q_res.n_cols, rxm = r*m;
        arma::vec rowindex(r);
        arma::vec colindex(r);
        arma::mat temp(r, r);
        arma::mat kernel(rxm, rxm);
                for (int i = 0; i < m; ++i) {
                        for (int j = 0; j < m; ++j) {
                                rowindex = arma::regspace(i, m, rxm);
                                colindex = arma::regspace(j, m, rxm);
                                temp = left_mat * Xc.t() * arma::diagmat(Q_res.col(i) % Q_res.col(j)) * Xc * left_mat.t()/n;
                                for (int p = 0; p < r; ++p) {
                                        for (int q = 0; q < r; ++q){
                                                kernel(rowindex(p), colindex(q)) = temp(p, q);
                                        }
                                }
                        }
                }
        return kernel;
}
