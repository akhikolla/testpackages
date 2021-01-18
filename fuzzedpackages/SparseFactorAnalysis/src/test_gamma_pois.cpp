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

double find_pnorm(arma::vec x) {
  using namespace arma;
  
  double tmp_a = x(0);
  double tmp_b = x(1);

  // Rcpp::Rcout << "inside find_pnorm?" << endl;
  
  vec coefs(8);
  coefs << -1.82517672 << 0.51283415 << -0.81377290 << -0.02699400 << -0.49642787 << -0.33379312 << -0.24176661 << 0.03776971;

  // Rcpp::Rcout << "inside find_pnorm? 2" << endl;
  
  vec tmp_x(8);
  tmp_x(0) = 1; tmp_x(1) = tmp_a; tmp_x(2) = tmp_b; tmp_x(3) = pow(tmp_a, 2);
  tmp_x(4) = pow(tmp_b, 2); tmp_x(5) = log(abs(tmp_a - tmp_b));
  tmp_x(6) = pow(log(abs(tmp_a - tmp_b)), 2); tmp_x(7) = tmp_a * tmp_b;

  return as_scalar(tmp_x.t() * coefs);

}


arma::field<arma::vec> test_gamma_pois(arma::vec gamma_try,
				       arma::mat Theta_last_0,
				       arma::imat votes_mat,
				       arma::vec cutoff_seq
				       ) {
  
  using namespace arma;

  // Rcpp::Rcout << "gamma_try: " << gamma_try << endl;
  // Rcpp::Rcout << "cutoff_seq: " << cutoff_seq << endl;

  // 
  // inputs

  uvec tmp_ind;
  vec dev_out(1);
  uvec::iterator a;
  uvec::iterator b;

  if (any(gamma_try > 2)) {
    tmp_ind = find(gamma_try > 2);
    a = tmp_ind.begin();
    b = tmp_ind.end();
    for (uvec::iterator i=a; i!=b; ++i) {
      gamma_try(*i) = 2;
    }
  }

  if (any(gamma_try < -2)) {
    tmp_ind = find(gamma_try < -2);
    a = tmp_ind.begin();
    b = tmp_ind.end();
    for (uvec::iterator i=a; i!=b; ++i) {
      gamma_try(*i) = -2;
    }
  }

  gamma_try = exp(gamma_try);

  vec count_seq;
  if (cutoff_seq.n_elem == 0) {
    count_seq = linspace<vec>(-1, votes_mat.max() + 2, 1);
  } else {
    count_seq = cutoff_seq;
  }
  
  vec taus_try = count_seq * 0;
  tmp_ind = uvec(find(vectorise(votes_mat) == 0));
  double a_int = static_cast<double>(tmp_ind.n_elem) / static_cast<double>(votes_mat.n_elem);

  tmp_ind = uvec(find(count_seq >= 0));
  taus_try(tmp_ind) = R::qnorm(a_int, 0, 1, 1, 0) + gamma_try(0) * vectorise(pow(count_seq(tmp_ind), gamma_try(1)));

  // Rcpp::Rcout << "count_seq: " << count_seq.t();

  tmp_ind = uvec(find(count_seq < 0));

  a = tmp_ind.begin();
  b = tmp_ind.end();
  for (uvec::iterator i=a; i!=b; ++i) {
    taus_try(*i) = - datum::inf;
  }

  vec tmp_vec = taus_try(uvec(find_finite(taus_try)));
  tmp_vec = sort(tmp_vec);
  vec tmp_vec2 = taus_try(uvec(find_nonfinite(taus_try)));

  taus_try.tail(tmp_vec.n_elem) = tmp_vec;
  taus_try.head(tmp_vec2.n_elem) = tmp_vec2;
  taus_try(0) = - datum::inf; 
  // taus_try(0) = -100000;

  vec lik(votes_mat.n_elem);
  vec log_lik(votes_mat.n_elem);
  vec theta_vec = vectorise(Theta_last_0);
  int i_tmp = 0;
  uvec votes_vec = conv_to<uvec>::from(vectorise(votes_mat));
  a = votes_vec.begin();
  b = votes_vec.end();
  for (uvec::iterator i=a; i!=b; ++i) {
    lik(i_tmp) = R::pnorm(taus_try(*i + 2 - 1) - theta_vec(i_tmp), 0, 1, 1, 0) - R::pnorm(taus_try(*i + 1 - 1) - theta_vec(i_tmp), 0, 1, 1, 0);
    i_tmp++;
  }

  log_lik = log(lik);

  tmp_ind = uvec(find_nonfinite(log_lik));

  // Rcpp::Rcout << "tmp_ind.n_elem: " << tmp_ind.n_elem << endl;
  // if (tmp_ind.n_elem > 0) {
    mat tmp_mat(tmp_ind.n_elem, 2);
    tmp_mat.col(0) = taus_try(votes_vec(tmp_ind) + 2 - 1) - theta_vec(tmp_ind);
    tmp_mat.col(1) = taus_try(votes_vec(tmp_ind) + 1 - 1) - theta_vec(tmp_ind);

  //   Rcpp::Rcout << votes_vec(tmp_ind(0)) + 2 << endl;
  //   Rcpp::Rcout << votes_vec(tmp_ind(0)) + 1 << endl;
  // }

  // Rcpp::Rcout << "where?0" << endl;

  int mat_b = tmp_mat.n_rows;
  // Rcpp::Rcout << "mat_b = " << mat_b << endl;
  i_tmp = 0;
  for (int j=0; j<mat_b; ++j) {
    // Rcpp::Rcout << "j = " << j << endl;
    // Rcpp::Rcout << "tmp_ind(i_tmp) = " << tmp_ind(i_tmp) << endl;
    // if (j == 0) {
    //   // Rcpp::Rcout << "tmp_mat.row(j).t(): " << endl;
    //   // Rcpp::Rcout << tmp_mat.row(j).t();
    //   // Rcpp::Rcout << "tmp_mat.row(j).t() end" << endl;
    // }
    log_lik(tmp_ind(i_tmp)) = find_pnorm(tmp_mat.row(j).t());
    i_tmp++;
  }

  tmp_ind = uvec(find_nonfinite(log_lik));

  // Rcpp::Rcout << "tmp_ind 2: " << endl;
  // Rcpp::Rcout << log_lik(tmp_ind(0));
  // Rcpp::Rcout << "tmp_ind 2 end" << endl;


  // Rcpp::Rcout << "where?" << endl;
  
  double threash = -1e30;
  if (threash < min(log_lik)) threash = min(log_lik);
  // if (any(log_lik < threash)) {
  log_lik.elem( uvec(find(log_lik < threash)) ).fill(threash);
  // }

  // Rcpp::Rcout << "where? 2" << endl;

  dev_out(0) = -2 * sum(log_lik) + sum(pow(log(gamma_try), 2));


  // output
  field<vec> out(2);
  out(0) = dev_out;
  out(1) = taus_try;


  return out;

}
