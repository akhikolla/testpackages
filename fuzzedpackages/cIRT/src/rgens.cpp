#include <RcppArmadillo.h>
#include "rgens.h"

//' Generate Random Multivariate Normal Distribution
//'
//' Creates a random Multivariate Normal when given number of obs, mean,
//' and sigma. 
//' 
//' @param n  An `integer`, which gives the number of observations.  (> 0) 
//' @param mu A `vector` length m that represents the means of the normals. 
//' @param S  A `matrix` with dimensions m x m that provides Sigma,
//'           the covariance matrix. 
//' @return 
//' A `matrix` that is a Multivariate Normal distribution.
//' 
//' @author 
//' James Joseph Balamuta 
//' 
//' @seealso
//' [TwoPLChoicemcmc()] and [probitHLM()]
//' 
//' @export
//' @examples
//' # Call with the following data: 
//' rmvnorm(2, c(0,0), diag(2))
// [[Rcpp::export]]
arma::mat rmvnorm(unsigned int n, const arma::vec &mu, const arma::mat &S)
{
    unsigned int ncols = S.n_cols;
    arma::mat Y(n, ncols);
    Y.imbue(norm_rand);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(S);
}

//' Generate Random Wishart Distribution
//'
//' Creates a random wishart distribution when given degrees of freedom
//' and a sigma matrix. 
//' 
//' @param df An `integer`, which gives the degrees of freedom of the Wishart.
//'           (> 0) 
//' @param S  A `matrix` with dimensions m x m that provides Sigma, the
//'           covariance matrix. 
//'           
//' @return 
//' A `matrix` that is a Wishart distribution, aka the sample covariance
//' matrix of a Multivariate Normal Distribution 
//' 
//' @author 
//' James Joseph Balamuta
//'  
//' @seealso 
//' [riwishart()] and [probitHLM()]
//' 
//' @export 
//' @examples 
//' # Call with the following data:
//' rwishart(3, diag(2))
//'
//' # Validation
//' set.seed(1337)
//' S = toeplitz((10:1)/10)
//' n = 10000
//' o = array(dim = c(10,10,n))
//' for(i in 1:n){
//' o[,,i] = rwishart(20, S)
//' }
//' mR = apply(o, 1:2, mean)
//' Va = 20*(S^2 + tcrossprod(diag(S)))
//' vR = apply(o, 1:2, var)
//' stopifnot(all.equal(vR, Va, tolerance = 1/16))
// [[Rcpp::export]]
arma::mat rwishart(unsigned int df, const arma::mat &S)
{
    // Dimension of returned wishart
    unsigned int m = S.n_rows;

    // Z composition:
    // sqrt chisqs on diagonal
    // random normals below diagonal
    // misc above diagonal
    arma::mat Z(m, m);

    // Fill the diagonal
    for (unsigned int i = 0; i < m; i++) {
        Z(i, i) = sqrt(R::rchisq(df - i));
    }

    // Fill the lower matrix with random guesses
    for (unsigned int j = 0; j < m; j++) {
        for (unsigned int i = j + 1; i < m; i++) {
            Z(i, j) = R::rnorm(0, 1);
        }
    }

    // Lower triangle * chol decomp
    arma::mat C = arma::trimatl(Z).t() * arma::chol(S);

    // Return random wishart
    return C.t() * C;
}

//' Generate Random Inverse Wishart Distribution
//'
//' Creates a random inverse wishart distribution when given
//' degrees of freedom and a sigma matrix. 
//' 
//' @param df An `integer` that represents the degrees of freedom. (> 0) 
//' @param S  A `matrix` with dimensions m x m that provides Sigma, 
//'           the covariance matrix. 
//' 
//' @return 
//' A `matrix` that is an inverse wishart distribution. 
//' 
//' @author
//' James Joseph Balamuta 
//'  
//' @seealso 
//' [rwishart()] and [TwoPLChoicemcmc()]
//' 
//' @export 
//' @examples
//' #Call with the following data:
//' riwishart(3, diag(2))
// [[Rcpp::export]]
arma::mat riwishart(unsigned int df, const arma::mat &S)
{
    return rwishart(df, S.i()).i();
}
