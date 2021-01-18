#ifndef test_gamma_pois_h
#define test_gamma_pois_h

double find_pnorm(arma::vec x);

arma::field<arma::vec> test_gamma_pois(arma::vec gamma_try, arma::mat Theta_last_0, arma::imat votes_mat, arma::vec cutoff_seq);

#endif
