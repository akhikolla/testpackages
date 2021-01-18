// [[Rcpp::depends(RcppArmadillo)]]
                                     /* this is called an Rcpp attribute, it is
                                        a comment from the perspective of C++,
                                        but sourceCpp() knows to look at them
                                        and act in a certain way when it sees
                                        them

                                        the intervention here is making sure the
                                        RcppArmadillo.h file can be found
                                      */

#include <RcppArmadillo.h>

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "test_gamma_pois.h"

arma::vec grad_est_comb(arma::vec x,
			arma::vec params,
			double del,
			arma::mat Theta_last_0,
			arma::imat votes_mat,
			arma::vec cutoff_seq) {
  using namespace arma;
  double dev1_p, dev1_n, dev2_p, dev2_n, dev0;

  vec tmp_x(2);

  tmp_x << x(0) + del << params(1) << endr;

  // Rcpp::Rcout << "tmp_x: " << tmp_x << endl;
  
  field<vec> tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);

  dev1_p = (tmp_out(0))(0);

  tmp_x(0) -= 2 * del;

  // Rcpp::Rcout << "test_gamma_pois 1 output: " << dev1_p << endl;

  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);

  dev1_n = (tmp_out(0))(0);

  // Rcpp::Rcout << "test_gamma_pois 2 output" << dev1_n << endl;

  tmp_x(0) = params(0);
  tmp_x(1) = x(1) + del;

  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);

  dev2_p = (tmp_out(0))(0);

  // Rcpp::Rcout << "test_gamma_pois 3 output: " << dev2_p << endl;

  tmp_x(0) -= 2 * del;

  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);

  dev2_n = (tmp_out(0))(0);

  // Rcpp::Rcout << "test_gamma_pois 4 output: " << dev2_n << endl;

  tmp_out = test_gamma_pois(x, Theta_last_0, votes_mat, cutoff_seq);

  dev0 = (tmp_out(0))(0);

  // Rcpp::Rcout << "test_gamma_pois5 output: " << dev0 << endl;
  
  vec out(4);
  out(0) = (dev1_p - dev1_n) / (4 * del);
  out(1) = (dev2_p - dev2_n) / (4 * del);
  out(2) = (dev1_p - 2 * dev0 + dev1_n) / (2 * del);
  out(3) = (dev2_p - 2 * dev0 + dev2_n) / (2 * del);

  return out;

}


arma::vec grad_est(arma::vec x,
		   arma::vec params,
		   double del,
		   arma::mat Theta_last_0,
		   arma::imat votes_mat,
		   arma::vec cutoff_seq) {
  
  using namespace arma;
  double grad1;
  double grad2;


  vec tmp_x(2);
  tmp_x << x(0) + del << params(1) << endr;
  field<vec> tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad1 = (tmp_out(0))(0);
  tmp_x(0) -= 2 * del;
  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad1 -= (tmp_out(0))(0);
  grad1 /= (4 * del);

  tmp_x(0) = params(0);
  tmp_x(1) = x(1) + del;
  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad2 = (tmp_out(0))(0);
  tmp_x(1) -= 2 * del;
  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad2 -= (tmp_out(0))(0);
  grad2 /= (4 * del);

  vec out;
  out << grad1 << grad2 << endr;

  return out;

}

double grad_est_1(arma::vec x,
		  arma::vec params,
		  double del,
		  arma::mat Theta_last_0,
		  arma::imat votes_mat,
		  arma::vec cutoff_seq) {
  
  using namespace arma;
  double grad1;

  vec tmp_x(2);
  tmp_x << x(0) + del << params(1) << endr;
  field<vec> tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad1 = (tmp_out(0))(0);
  tmp_x(0) -= 2 * del;
  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad1 += (tmp_out(0))(0);
  tmp_out = test_gamma_pois(x, Theta_last_0, votes_mat, cutoff_seq);
  grad1 -= 2 * (tmp_out(0))(0);
  grad1 /= (2 * del);

  

  return grad1;

}

double grad_est_2(arma::vec x,
		  arma::vec params,
		  double del,
		  arma::mat Theta_last_0,
		  arma::imat votes_mat,
		  arma::vec cutoff_seq) {
  
  using namespace arma;
  double grad2;

  vec tmp_x(2);
  tmp_x << params(0) << x(1) + del << endr;
  field<vec> tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad2 = (tmp_out(0))(0);
  tmp_x(1) -= 2 * del;
  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad2 += (tmp_out(0))(0);
  tmp_out = test_gamma_pois(x, Theta_last_0, votes_mat, cutoff_seq);
  grad2 -= 2 * (tmp_out(0))(0);
  grad2 /= (2 * del);

  return grad2;

}

double grad_est_12(arma::vec x,
		   arma::vec params,
		   double del,
		   arma::mat Theta_last_0,
		   arma::imat votes_mat,
		   arma::vec cutoff_seq) {
  
  using namespace arma;
  double grad12;

  vec tmp_x(2);
  tmp_x = x + del;
  field<vec> tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad12 = (tmp_out(0))(0);
  tmp_x(1) -= 2 * del;
  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad12 -= (tmp_out(0))(0);
  tmp_x(0) -= 2 * del;
  tmp_x(1) += 2 * del;
  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad12 -= (tmp_out(0))(0);
  tmp_x = x - del;
  tmp_out = test_gamma_pois(tmp_x, Theta_last_0, votes_mat, cutoff_seq);
  grad12 += (tmp_out(0))(0);
  grad12 /= (4 * del);

  return grad12;

}

