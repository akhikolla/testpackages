#ifndef grad_est_h
#define grad_est_h

arma::vec grad_est_comb(arma::vec x, arma::vec params, double del, arma::mat Theta_last_0, arma::imat votes_mat, arma::vec cutoff_seq);

arma::vec grad_est(arma::vec x, arma::vec params, double del, arma::mat Theta_last_0, arma::imat votes_mat, arma::vec cutoff_seq);

double grad_est_1(arma::vec x, arma::vec params, double del, arma::mat Theta_last_0, arma::imat votes_mat, arma::vec cutoff_seq);

double grad_est_2(arma::vec x, arma::vec params, double del, arma::mat Theta_last_0, arma::imat votes_mat, arma::vec cutoff_seq);

double grad_est_12(arma::vec x, arma::vec params, double del, arma::mat Theta_last_0, arma::imat votes_mat, arma::vec cutoff_seq);

#endif
