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
#include "grad_est.h"
#include "test_gamma_pois.h"


// [[Rcpp::export]]
Rcpp::List hmc4sfa(int num_HMC,
		   arma::vec q,
		   arma::vec p,
		   arma::vec params,
		   double del,
		   double step_size,
		   arma::mat Theta_last_0,
		   arma::imat votes_mat,
		   arma::vec cutoff_seq,
		   double current_U,
		   arma::vec q0
		   ) {

  using namespace arma;
  using namespace Rcpp;
  using std::endl;

  mat epsilon(2, 2);
  mat hess_mat(2, 2);
  uvec tmp_ind;
  tmp_ind << 1 << 2 << endr;
  vec tmp_out(4);
  field<vec> tmp_out2;
  vec tmp_vec(2);
  bool revert_q;

  Rcout << "Update hessian here" << endl;

  for (int i_HMC = 0; i_HMC < num_HMC; ++i_HMC) {
    tmp_out = grad_est_comb(q, params, del, Theta_last_0, votes_mat, cutoff_seq);
    // Rcout << "tmp_out 1" << endl;
    // Rcout << tmp_out << endl;
    hess_mat(0, 0) = tmp_out(2);
    hess_mat.elem(tmp_ind).fill(grad_est_12(q, params, del, Theta_last_0, votes_mat, cutoff_seq));
    hess_mat(1, 1) = tmp_out(3);
    hess_mat.elem(uvec(find_nonfinite(hess_mat))).fill(0);
    hess_mat = .5 * (hess_mat + hess_mat.t());
    epsilon = pinv(hess_mat) * step_size / 10;
    // Rcout << "pinv * step_size done" << endl;
    epsilon.elem(uvec(find_nonfinite(epsilon))).fill(0);

    q = q + epsilon * p;
    tmp_vec = grad_est(q, params, del, Theta_last_0, votes_mat, cutoff_seq);
    tmp_vec(uvec(find_nonfinite(tmp_vec))).fill(0);
    p = p - epsilon * tmp_vec;
    
    tmp_out2 = test_gamma_pois(q, Theta_last_0, votes_mat, cutoff_seq);

    if (abs( (tmp_out2(0))(0) / 2 - current_U) > 10000) {
      revert_q = true;
      q = q0;
      break;
    }
    revert_q = false;

  }
  q = q + epsilon * p;

  tmp_vec = grad_est(q, params, del, Theta_last_0, votes_mat, cutoff_seq);
  tmp_vec(uvec(find_nonfinite(tmp_vec))).fill(0);

  p = p - ((epsilon * tmp_vec) / 2);
  p = -p;

  if (revert_q) q = q0;

  return List::create(_["p"] = p,
		      _["q"] = q);
}
