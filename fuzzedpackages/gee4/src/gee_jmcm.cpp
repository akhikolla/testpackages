#include "gee_jmcm.h"
#include "utils.h"

namespace gee {

  //// NEW FUNCTION 
  void gee_jmcm::UpdateGEES(const arma::vec &x) {
    int debug = 0;

    bool update = true;         // update or not?

    switch (free_param_) {
    case 0:
      if (arma::min(x == tht_) == 1) update = false;
      break;
    case 1:
      if (arma::min(x == bta_) == 1) update = false;
      break;
    case 2:
      if (arma::min(x == lmd_) == 1) update = false;
      break;
    case 3:
      if (arma::min(x == gma_) == 1) update = false;
      break;
    default:
      Rcpp::Rcout << "Wrong value for free_param_" << std::endl;
    }

    if (update) {
      UpdateParam(x);
      UpdateModel();
    } else {
      if (debug) Rcpp::Rcout << "Hey, I did save some time!:)" << std::endl;
    }
  }

  //// NEW FUNCTION
  void gee_jmcm::UpdateParam(const arma::vec &x) {
    int n_bta = X_.n_cols;
    int n_lmd = Z_.n_cols;
    int n_gma = W_.n_cols;

    switch (free_param_) {
    case 0:
      tht_ = x;
      bta_ = x.rows(0, n_bta - 1);
      lmd_ = x.rows(n_bta, n_bta + n_lmd - 1);
      gma_ = x.rows(n_bta + n_lmd, n_bta + n_lmd + n_gma - 1);
      break;

    case 1:
      tht_.rows(0, n_bta - 1) = x;
      bta_ = x;
      break;

    case 2:
      tht_.rows(n_bta, n_bta + n_lmd - 1) = x;
      lmd_ = x;
      break;

    case 3:
      tht_.rows(n_bta + n_lmd, n_bta + n_lmd + n_gma - 1) = x;
      gma_ = x;
      break;

    default:
      Rcpp::Rcout << "Wrong value for free_param_" << std::endl;
    }
  }

  //// NEW FUNCTION
  void gee_jmcm::UpdateModel() {
    int debug = 0;

    switch (free_param_) {
    case 0:
      // if (cov_only_) Xbta_ = mean_;
      // else Xbta_ = X_ * beta_;
      Xbta_ = X_ * bta_;
      Zlmd_ = Z_ * lmd_;
      Wgma_ = W_ * gma_;
      Resid_ = Y_ - Xbta_;
      if (debug) Rcpp::Rcout << "update Xbta Zlmd Wgam r" << std::endl;
      
      /*
        if (debug) Rcpp::Rcout << "UpdateG(x)" << std::endl;
        UpdateG();
        if (debug) Rcpp::Rcout << "UpdateTResid(x)" << std::endl;
        UpdateTResid();
        if (debug) Rcpp::Rcout << "Update Finished.." << std::endl;
      */
      break;

    case 1:
      // if (cov_only_) Xbta_ = mean_;
      // else Xbta_ = X_ * beta_;
      Xbta_ = X_ * bta_;
      Resid_ = Y_ - Xbta_;
      if (debug) Rcpp::Rcout << "gees_jmcm::UpdateModel(): Resid_ updated" << std::endl;

      /*
        UpdateG();
        UpdateTResid();
      */
      break;

    case 2:
      Zlmd_ = Z_ * lmd_;
      if (debug) Rcpp::Rcout << "gees_jmcm::UpdateModel(): Zlmd_ updated" << std::endl;
      
      break;

    case 3:
      Wgma_ = W_ * gma_;
      if (debug) Rcpp::Rcout << "gees_jmcm::UpdateModel(): Wgma_ updated" << std::endl;
      
      /*
        UpdateTResid();
      */
      break;

    default:
      Rcpp::Rcout << "Wrong value for free_param_" << std::endl;
    }
  }
  
  //// NEW FUNCTION 
  void gee_jmcm::UpdateBeta() {
    arma::uword n_sub = m_.n_elem;
    arma::uword n_bta = X_.n_cols;
    
    arma::mat bta_lhs = arma::zeros<arma::mat>(n_bta, n_bta);
    arma::vec bta_rhs = arma::zeros<arma::vec>(n_bta);
    for (arma::uword i = 0; i != n_sub; ++i) {
      arma::vec ri = get_Resid(i);
      arma::mat Sigmai_inv = get_Sigma_inv(i);
      if (use_ipw_) {
        arma::mat H_sqrt = get_weights_sqrt(i);
        Sigmai_inv = H_sqrt * Sigmai_inv * H_sqrt;
      }
      
      arma::mat deriv1;
      arma::vec yi_tilde;
      if (link_mode_ == identity_link) {
        deriv1 = get_X(i).t();
        yi_tilde = get_Y(i);
      }
      bta_lhs += deriv1 * Sigmai_inv * deriv1.t();
      bta_rhs += deriv1 * Sigmai_inv * yi_tilde;
    }

    arma::vec beta = bta_lhs.i() * bta_rhs;

    set_beta(beta);
  }

  /////// TO DO!!!!
  void gee_jmcm::UpdateLambda() {
    arma::uword n_sub = m_.n_elem;
    arma::uword n_lmd = Z_.n_cols;
    
    arma::mat lmd_lhs = arma::zeros<arma::mat>(n_lmd, n_lmd);
    arma::vec lmd_rhs = arma::zeros<arma::vec>(n_lmd);
    for (arma::uword i = 0; i != n_sub; ++i) {
      arma::uword mi = m_(i);
      arma::vec ri = get_Resid(i);
      arma::mat Ti = get_T(i);
      arma::mat Di = get_D(i);
      arma::vec di = Di.diag();
      arma::vec di2 = arma::pow(Di.diag(), 1);
      arma::vec logdi2 = arma::log(arma::pow(Di.diag(), 1));
      arma::mat Di_inv = arma::diagmat(arma::pow(Di.diag(), -1));
    
      arma::mat Ai_sqrt_inv = 1 / arma::datum::sqrt2 * Di_inv;
      arma::mat Ri_inv;
    
      if (corr_mode_ == Identity_corr) Ri_inv = arma::eye(mi, mi);
      else if (corr_mode_ == CompSymm_corr) Ri_inv = dragonwell::corr_cs(rho_, mi).i();
      else if (corr_mode_ == AR1_corr) Ri_inv = dragonwell::corr_ar1(rho_, mi).i();
    
      arma::mat cov_inv = Ai_sqrt_inv * Ri_inv * Ai_sqrt_inv;
      if (use_ipw_) {
        arma::mat H_sqrt = get_weights_sqrt(i);
        //H_sqrt = H_sqrt * H_sqrt * H_sqrt;
        cov_inv = H_sqrt * cov_inv * H_sqrt;
      }

      arma::mat deriv2 = get_Z(i).t();
      for (arma::uword t = 1; t <= mi; ++t) {
        deriv2.col(t - 1) *= arma::as_scalar(di(t - 1));
      }
    
      arma::vec epsi2 = arma::pow(Ti * ri, 2);
      arma::vec epsi2_tilde = epsi2 - di2 + Di * logdi2;
    
      lmd_lhs += deriv2 * cov_inv * deriv2.t();
      lmd_rhs += deriv2 * cov_inv * epsi2_tilde;
    }

    arma::vec lambda = lmd_lhs.i() * lmd_rhs;

    set_lambda(lambda);    
  }

  //////// TO DO!!!!
  void gee_jmcm::UpdateGamma() {
    arma::uword n_sub = m_.n_elem;
    arma::uword lgma = W_.n_cols;

    arma::mat gma_lhs = arma::zeros<arma::mat>(lgma, lgma);
    arma::vec gma_rhs = arma::zeros<arma::vec>(lgma);
    for (arma::uword i = 0; i != n_sub; ++i) {
      arma::uword mi = m_(i);
      arma::mat Wi = get_W(i);
      arma::vec ri = get_Resid(i);
      arma::mat Di = get_D(i);
      arma::mat Di_inv = arma::diagmat(arma::pow(Di.diag(), -1));
      if (use_ipw_) {
        arma::mat H_sqrt = get_weights_sqrt(i);
        //H_sqrt = H_sqrt * H_sqrt * H_sqrt;
        Di_inv = H_sqrt * Di_inv * H_sqrt;
      }

      arma::uword rindex = 0;
      
      arma::mat deriv3 = arma::zeros<arma::mat>(lgma, mi);
      for (arma::uword j = 2; j <= mi; ++j) {
        for (arma::uword k = 1; k <= (j - 1); ++k) {
          deriv3.col(j - 1) += ri(k - 1) * Wi.row(rindex).t();
          ++rindex;
        }
      }
      
      gma_lhs += deriv3 * Di_inv * deriv3.t();
      gma_rhs += deriv3 * Di_inv * ri;
    }

    arma::vec gamma = gma_lhs.i() * gma_rhs;

    set_gamma(gamma);
  }

  arma::vec gee_jmcm::operator()(const arma::vec &x) {
    int debug = 0;

    if (debug) Rcpp::Rcout << "gee_jmcm::operator(): initializing..." << std::endl;
    
    set_params(x);

    arma::uword n_sub = m_.n_elem;
    arma::uword lbta = X_.n_cols;
    arma::uword llmd = Z_.n_cols;
    arma::uword lgma = W_.n_cols;

    arma::vec gee_bta = arma::zeros<arma::vec>(lbta);
    arma::vec gee_lmd = arma::zeros<arma::vec>(llmd);
    arma::vec gee_gma = arma::zeros<arma::vec>(lgma);

    if (debug) Rcpp::Rcout << "gee_jmcm::operator(): calculating three GEEs..." << std::endl; 
    for (arma::uword i = 0; i < n_sub; ++i) {
      // if (debug) Rcpp::Rcout << "iter " << i << std::endl;

      arma::uword mi = m_(i);
      arma::vec ri = get_Resid(i);
      arma::mat Ti = get_T(i);
      arma::mat Di = get_D(i);
      arma::vec di = Di.diag();
      arma::vec di2 = arma::pow(Di.diag(), 1);
      arma::vec logdi2 = arma::log(arma::pow(Di.diag(), 1));
      arma::mat Di_inv = arma::diagmat(arma::pow(Di.diag(), -1));
      arma::mat Sigmai_inv = get_Sigma_inv(i);
 
      arma::mat deriv1;
      if (link_mode_ == identity_link) {
        deriv1 = get_X(i).t();
      }
      gee_bta += deriv1 * Sigmai_inv * ri;

      // if (debug) Rcpp::Rcout << "update gee_gma" << std::endl;

      arma::mat Ai_sqrt_inv = 1 / arma::datum::sqrt2 * Di_inv;
      arma::mat Ri_inv;
      //if (debug) Ai_sqrt_inv.print("Ai_inv = ");
      if (corr_mode_ == Identity_corr) Ri_inv = arma::eye(mi, mi);
      if (corr_mode_ == CompSymm_corr) Ri_inv = dragonwell::corr_cs(rho_, mi).i();
      //  Ri_inv = dragonwell::corr_cs(rho_, mi).i();
      if (corr_mode_ == AR1_corr) Ri_inv = dragonwell::corr_ar1(rho_, mi).i();
      arma::mat cov_inv = Ai_sqrt_inv * Ri_inv * Ai_sqrt_inv;
      arma::mat deriv2 = get_Z(i).t();
      for (arma::uword t = 1; t <= mi; ++t) {
        deriv2.col(t - 1) *= arma::as_scalar(di(t - 1));
      }
      
      arma::vec epsi2 = arma::pow(Ti * ri, 2);
      arma::vec epsi2_tilde = epsi2 - di2 + Di * logdi2;
      

      gee_lmd += deriv2 * cov_inv * (epsi2 - di2);
      // if (debug) Rcpp::Rcout << "update gee_lmd" << std::endl;

      arma::mat Wi = get_W(i);
      arma::uword rindex = 0;
      arma::mat deriv3 = arma::zeros<arma::mat>(lgma, mi);
      for (arma::uword j = 2; j <= mi; ++j) {
        for (arma::uword k = 1; k <= (j - 1); ++k) {
          deriv3.col(j - 1) += ri(k - 1) * Wi.row(rindex).t();
          ++rindex;
        }
      }

      gee_gma += deriv3 * Di_inv * (Ti * ri);
      // if (debug) Rcpp::Rcout << "update gee_gma" << std::endl;
      // if (debug) Rcpp::Rcout << "m length = " << m_.n_elem << std::endl;

    }
    arma::vec result = dragonwell::join_vecs({gee_bta,gee_lmd,gee_gma});
    return result;
  }

  bool gee_jmcm::learn(const arma::uvec& m, const arma::mat& Y,
                       const arma::mat& X, const arma::mat& Z, const arma::mat& W,
                       const gee_link_mode& link_mode,
                       const gee_corr_mode& corr_mode, const double rho,
                       const arma::vec& start, const arma::uword fs_iter,
                       const bool print_mode) {
    int debug = 1;

    m_ = m;
    Y_ = Y;
    X_ = X;
    Z_ = Z;
    W_ = W;

    if (debug) Rcpp::Rcout << "initialization" << std::endl;

    const bool link_mode_ok = (link_mode == identity_link);
    const bool corr_mode_ok =
      (corr_mode == Identity_corr) || (corr_mode == CompSymm_corr) || (corr_mode == AR1_corr);

    if (!link_mode_ok)
      Rcpp::Rcerr << "gee_jmcm::learn(): unknown link_mode" << std::endl;
    if (!corr_mode_ok)
      Rcpp::Rcerr << "gee_jmcm::learn(): unknown corr_mode" << std::endl;

    if (fs_iter > 0) {
      bool status = false;
      status = fs_iterate(link_mode, corr_mode, rho, start, fs_iter, print_mode);
      if (status == false) {
        Rcpp::Rcerr << "gee_jmcm::learn(): quasi-Fisher scoring algorithm failed"
                    << std::endl;
      }
    }

    return true;
  }

  bool gee_jmcm::fs_iterate(const gee_link_mode& link_mode,
                            const gee_corr_mode& corr_mode, const double rho,
                            const arma::vec& start, arma::uword max_iter,
                            const bool verbose) {
    int debug = 0;
    arma::vec x = start;
    if (debug) Rcpp::Rcout << "entering set_params()" << std::endl;
    set_params(start);

    if (verbose) {
      Rcpp::Rcout << "0: " << std::endl;
      x.t().print();
    }

    const double kEpsilon =
      std::numeric_limits<double>::epsilon();  // Machine precision
    const double kTolX = 4 * kEpsilon;  // Convergence criterion on x values
    const int n_pars = x.n_rows;        // number of parameters
    arma::vec p;
    for (arma::uword iter = 1; iter <= max_iter; ++iter) {
      arma::vec x2 = x;
      if (debug) Rcpp::Rcout << "before fs_update" << std::endl;
      fs_update_params(link_mode, corr_mode, rho);
      if (debug) Rcpp::Rcout << "after fs_update" << std::endl;

      x = dragonwell::join_vecs({bta_, lmd_, gma_});
      if (verbose) {
        Rcpp::Rcout << iter << ": " << std::endl;
        x.t().print();
      }
      p = x - x2;
      // Test for convergence on Delta x
      double test = 0.0;
      for (int i = 0; i != n_pars; ++i) {
        double temp = std::abs(p(i)) / std::max(std::abs(x(i)), 1.0);
        if (temp > test) test = temp;
      }
      if (test < kTolX) {
        break;
      }
    }

    return true;
  }

  void gee_jmcm::fs_update_params(const gee_link_mode& link_mode,
                                  const gee_corr_mode& corr_mode,
                                  const double rho) {
    int debug = 0;
    arma::uword n_sub = m_.n_elem;
    arma::uword lbta = X_.n_cols;
    arma::uword llmd = Z_.n_cols;
    arma::uword lgma = W_.n_cols;

    bta_ = arma::zeros<arma::vec>(lbta);
    lmd_ = arma::zeros<arma::vec>(llmd);
    gma_ = arma::zeros<arma::vec>(lgma);

    arma::mat bta_lhs = arma::zeros<arma::mat>(lbta, lbta);
    arma::vec bta_rhs = arma::zeros<arma::vec>(lbta);
    arma::mat lmd_lhs = arma::zeros<arma::mat>(llmd, llmd);
    arma::vec lmd_rhs = arma::zeros<arma::vec>(llmd);
    arma::mat gma_lhs = arma::zeros<arma::mat>(lgma, lgma);
    arma::vec gma_rhs = arma::zeros<arma::vec>(lgma);
    if (debug) Rcpp::Rcout << "fs_update before for" << std::endl;
    for (arma::uword i = 0; i < n_sub; ++i) {
      if (debug) Rcpp::Rcout << "fs_update for update beta" << std::endl;

      arma::uword mi = m_(i);
      if (debug) Rcpp::Rcout << "i= " << i << " mi= " << mi << std::endl;
      arma::mat Sigmai_inv = get_Sigma_inv(i);
      arma::mat deriv1;
      arma::vec yi_tilde;
      if (link_mode == identity_link) {
        deriv1 = get_X(i).t();
        yi_tilde = get_Y(i);
      }
      bta_lhs += deriv1 * Sigmai_inv * deriv1.t();
      bta_rhs += deriv1 * Sigmai_inv * yi_tilde;

      if (debug) Rcpp::Rcout << "fs_update for update lmd" << std::endl;

      arma::vec ri = get_Resid(i);
      arma::mat Ti = get_T(i);
      arma::mat Di = get_D(i);
      arma::vec di = Di.diag();
      arma::vec logdi = arma::log(di);
      // arma::vec di2 = arma::pow(Di.diag(), 1);
      // arma::vec logdi2 = arma::log(arma::pow(Di.diag(), 1));
      arma::mat Di_inv = arma::diagmat(arma::pow(Di.diag(), -1));

      if (debug) Rcpp::Rcout << "fs_update for update gma" << std::endl;

      arma::mat Ai_sqrt_inv = 1 / arma::datum::sqrt2 * Di_inv;
      arma::mat Ri_inv;
      //if (debug) Ai_sqrt_inv.print("Ai_inv = ");
      if (corr_mode_ == Identity_corr) Ri_inv = arma::eye(mi, mi);
      if (corr_mode == CompSymm_corr) Ri_inv = dragonwell::corr_cs(rho_, mi).i();
      if (corr_mode == AR1_corr) Ri_inv = dragonwell::corr_ar1(rho_, mi).i();
      arma::mat cov_inv = Ai_sqrt_inv * Ri_inv * Ai_sqrt_inv;
      arma::mat deriv2 = get_Z(i).t();
      for (arma::uword t = 1; t <= mi; ++t) {
        deriv2.col(t - 1) *= arma::as_scalar(di(t - 1));
      }
      if (debug) Rcpp::Rcout << "fs_update for update gma part1" << std::endl;
      arma::vec epsi2 = arma::pow(Ti * ri, 2);
      arma::vec epsi2_tilde = (epsi2 - di) + Di * logdi;
      if (debug) Rcpp::Rcout << "fs_update for update gma part2" << std::endl;

      lmd_lhs += deriv2 * cov_inv * deriv2.t();
      lmd_rhs += deriv2 * cov_inv * epsi2_tilde;

      arma::mat Wi = get_W(i);
      arma::uword rindex = 0;
      arma::mat deriv3 = arma::zeros<arma::mat>(lgma, mi);
      for (arma::uword j = 2; j <= mi; ++j) {
        for (arma::uword k = 1; k <= (j - 1); ++k) {
          deriv3.col(j - 1) += ri(k - 1) * Wi.row(rindex).t();
          ++rindex;
        }
      }

      gma_lhs += deriv3 * Di_inv * deriv3.t();
      gma_rhs += deriv3 * Di_inv * ri;
    }

    bta_ = bta_lhs.i() * bta_rhs;
    lmd_ = lmd_lhs.i() * lmd_rhs;
    gma_ = gma_lhs.i() * gma_rhs;

    set_params(dragonwell::join_vecs({bta_, lmd_, gma_}));
  }

}  // namespace gee
