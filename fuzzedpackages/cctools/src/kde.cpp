#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec kern_epan_1d(const arma::vec& x,
                       const double& b)
{
    vec out(x);
    unsigned int n = x.n_rows;

    for (unsigned int i = 0; i < n; ++i) {
        if (std::fabs(x[i] / b) >= 1.0) {
            out[i] = 0;
        } else {
            out[i] = 3.0/4.0 * (1 - pow(x[i] / b, 2)) / b;
        }
    }
    return out;
}


arma::vec kern_epan(const arma::mat& x,
                    const arma::vec& b)
{
    vec out(x.n_rows);
    unsigned int d = x.n_cols;

    out = kern_epan_1d(x.col(0), b[0]);
    for (unsigned int j = 1; j < d; ++j) {
        out = out % kern_epan_1d(x.col(j), b[j]);
    }
    return out;
}

// [[Rcpp::export]]
arma::vec eval_mvkde(const arma::mat& xev,
                     const arma::mat& x,
                     const arma::vec& b)
{
    vec out(xev.n_rows);

    for (unsigned int i = 0; i < xev.n_rows; ++i) {
        out[i] = mean(kern_epan(x - repmat(xev.row(i), x.n_rows, 1), b));
    }

    return out;
}

// for (initially) discrete data
// [[Rcpp::export]]
arma::vec lcv_mvkde_disc(const arma::mat& x,
                         const arma::mat& x_cc,
                         const arma::vec& b)
{
    vec KK(1), tmpll(1);
    rowvec tmpx_cc(x.n_cols);
    mat xx = x_cc;
    double n = x.n_rows;

    KK[0] = 0.0;
    for (unsigned int i = 0; i < n; ++i) {
        tmpx_cc = xx.row(i);
        xx.shed_row(i);
        tmpll = log(eval_mvkde(x.row(i), xx, b));
        tmpll[0] = fmax(tmpll[0], -15);  // ensure that values remain finite
        KK += tmpll / n;
        xx.insert_rows(i, tmpx_cc);
    }

    return KK;
}

// other bandwidth selection methods below (unused) --------------------

// arma::vec convkern_epan_1d(const arma::vec& x,
//                            const double& b)
// {
//     vec out(x);
//     unsigned int n = x.n_rows;
//
//     for (unsigned int i = 0; i < n; ++i) {
//         if (std::fabs(x[i] / b) >= 2.0) {
//             out[i] = 0;
//         } else {
//             out[i] = 3.0/160.0 * pow(2.0 - x[i] / b, 3) * (pow(x[i] / b, 2) + 6 * x[i] + 4) / b;
//         }
//     }
//     return out;
// }
//
// arma::vec convkern_epan(const arma::mat& x,
//                         const arma::vec& b)
// {
//     vec out(x.n_rows);
//     unsigned int d = x.n_cols;
//
//     out = kern_epan_1d(x.col(0), b[0]);
//     for (unsigned int j = 1; j < d; ++j) {
//         out = out % convkern_epan_1d(x.col(j), b[j]);
//     }
//     return out;
// }

// //' @export
// // see Bruce Hansen's lecture notes (Section 2.16)
// // [[Rcpp::export]]
// arma::vec lscv_mvkde(const arma::mat& x,
//                      const arma::mat& b)
// {
//     vec tmp, crit(b.n_rows);
//     mat xx = x;
//     double n = x.n_rows;
//
//     for (unsigned int j = 0; j < b.n_rows; ++j) {
//         double KK = 0.0;
//         for (unsigned int i = 0; i < x.n_rows; ++i) {
//             tmp = xx.row(i);
//             xx.shed_row(i);
//             KK += sum(convkern_epan(xx - repmat(tmp, xx.n_rows, 1), b.row(j)));
//             KK -= sum(2 * kern_epan(xx - repmat(tmp, xx.n_rows, 1), b.row(j)));
//             xx.insert_rows(i, tmp);
//         }
//         crit[j] = KK / pow(n, 2);
//         crit[j] += pow(3.0/5.0, b.n_cols) / (n * prod(b.row(j)));
//     }
//
//     return crit;
// }
//
// //' @export
// // for (initially) discrete data
// // [[Rcpp::export]]
// arma::vec lscv_mvkde_disc(const arma::mat& x,
//                           const arma::mat& x_cc,
//                           const arma::mat& b,
//                           const arma::mat& grid)
// {
//     vec tmpx_cc, est, sq, KK(1), crit(b.n_rows);
//     mat xx = x_cc;
//     double n = x.n_rows;
//
//     for (unsigned int j = 0; j < b.n_rows; ++j) {
//         KK[0] = 0.0;
//         est = eval_mvkde(grid, x_cc, b.row(j));
//         sq = sum(pow(est, 2));
//         for (unsigned int i = 0; i < n; ++i) {
//             tmpx_cc = xx.row(i);
//             xx.shed_row(i);
//             KK += eval_mvkde(x.row(i), xx, b.row(j)) / n;
//             xx.insert_rows(i, tmpx_cc);
//         }
//         crit[j] = sq[0] - 2 * KK[0];
//     }
//
//     return crit;
// }
