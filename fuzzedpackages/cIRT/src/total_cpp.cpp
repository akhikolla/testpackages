#include <RcppArmadillo.h>
#include "total_cpp.h"

//' Calculate Tabulated Total Scores
//'
//' Internal function to -2LL
//' 
//' @param N An `integer`, which gives the number of observations.  (> 0)
//' @param J An `integer`, which gives the number of items.  (> 0)
//' @param Y A N by J `matrix` of item responses.
//' 
//' @return 
//' A vector of tabulated total scores.
//' 
//' @author 
//' Steven Andrew Culpepper
//' 
//' @export
// [[Rcpp::export]]
arma::uvec Total_Tabulate(unsigned int N, unsigned int J, const arma::mat Y)
{
    arma::vec T(N);
    T = sum(Y, 1);
    arma::uvec H = arma::hist(T, arma::linspace<arma::vec>(0, J, J + 1));
    return H;
}

//' Generate Observed Data from choice model
//'
//' Generates observed cognitive and choice data from the IRT-Thurstone model.
//' 
//' @param N                  An `integer` for the number of observations.
//' @param J                  An `integer` for the number of items.
//' @param K                  An `integer` for the number of paired comparisons.
//' @param theta              A `vector` of latent cognitive variables.
//' @param as                 A `vector` of length J with item discriminations.
//' @param bs                 A `vector` of length J with item locations.
//' @param zeta               A `matrix` with dimensions N x V containing
//'                           random parameter estimates. 
//' @param gamma              A `vector` with dimensions P x 1 containing
//'                           fixed parameter estimates, where \eqn{P = P_1 + P_2} 
//' @param X                  A `matrix` with dimensions N*K x P_1 containing 
//'                           fixed effect design matrix without theta. 
//' @param W                  A `matrix` with dimensions N*K x V containing
//'                           random effect variables. 
//' @param subject_ids        A `vector` with length NK x 1 containing
//'                           subject-choice IDs. 
//' @param unique_subject_ids A `vector` with length N x 1 containing
//'                           unique subject IDs. 
//' 
//' @return 
//' A \code{list} that contains: 
//' 
//' \describe{
//'    \item{\code{Y}}{A `matrix` of dimension N by J}
//'    \item{\code{C}}{A `vector` of length NK}
//' } 
//' 
//' @author 
//' Steven Andrew Culpepper and James Joseph Balamuta 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List Generate_Choice(unsigned int N, unsigned int J, unsigned int K,
                           const arma::vec &theta, const arma::vec &as,
                           const arma::vec &bs, const arma::mat &zeta,
                           const arma::vec &gamma, const arma::mat &X,
                           const arma::mat &W, const arma::vec &subject_ids,
                           const arma::vec &unique_subject_ids)
{

    // Generate Y matrix
    arma::mat Y = arma::zeros<arma::mat>(N, J);
    arma::vec oneN = arma::ones<arma::vec>(N);
    arma::vec oneK = arma::ones<arma::vec>(K);
    arma::vec oneJ = arma::ones<arma::vec>(J);
    arma::mat etaY = theta * as.t() - oneN * bs.t();
    arma::mat ZY = arma::randn<arma::mat>(N, J);
    // generating Y
    Y.elem(find(etaY > ZY)).ones();

    // Generate choice vector,Cs
    //# of cases and choices
    unsigned int nk = X.n_rows;
    //  arma::vec Cs= arma::zeros<arma::vec>(nk);
    arma::mat Cs = arma::zeros<arma::mat>(N, K);
    // random z's for c
    //  arma::vec ZC = arma::randn<arma::vec>(nk);
    arma::mat ZC = arma::randn<arma::mat>(N, K);
    arma::vec Wzeta_0(nk);
    arma::mat etaC_mat(N, K);

    for (unsigned int i = 0; i < N; i++) {
        arma::uvec current_subject = find(subject_ids == unique_subject_ids(i));
        arma::mat Wi = W.rows(current_subject);
        arma::mat Xi = X.rows(current_subject);
        Wzeta_0.elem(current_subject) = Wi * trans(zeta.row(i)); // nk x 1
        etaC_mat.row(i) = gamma.t() * Xi.t() + zeta.row(i) * Wi.t();
    }
    
    arma::vec etaC = X * gamma + Wzeta_0;
    Cs.elem(find(etaC_mat > ZC)).ones();

    //  arma::vec nC = Cs*oneK;
    //  arma::vec nY = Y*oneJ;
    arma::uvec nC = Total_Tabulate(N, K, Cs);
    arma::uvec nY = Total_Tabulate(N, J, Y);

    return Rcpp::List::create(Rcpp::Named("nC", nC), Rcpp::Named("nY", nY),
                              Rcpp::Named("Cs", Cs));
}
