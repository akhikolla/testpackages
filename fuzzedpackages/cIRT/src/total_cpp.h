#ifndef TOTAL_CPP_H
#define TOTAL_CPP_H
Rcpp::List Generate_Choice(unsigned int N, unsigned int J, unsigned int K,
                           const arma::vec &theta, const arma::vec &as,
                           const arma::vec &bs, const arma::mat &zeta,
                           const arma::vec &gamma, const arma::mat &X,
                           const arma::mat &W, const arma::vec &subject_ids,
                           const arma::vec &unique_subject_ids);

arma::uvec Total_Tabulate(unsigned int N, unsigned int J, const arma::mat Y);

#endif
