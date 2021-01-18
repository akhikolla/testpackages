#include "functions.h"

//function to draw from a multivariate normal
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat cholSigma)
{
    int ncols = cholSigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * cholSigma;
}

