#ifndef GEE_JMCM_H_
#define GEE_JMCM_H_

#include "utils.h"
#include <RcppArmadillo.h>

struct gee_link_mode {
  const arma::uword id_;
  inline explicit gee_link_mode(const arma::uword id) : id_(id) {}
};

inline bool operator==(const gee_link_mode &a, const gee_link_mode &b) {
  return (a.id_ == b.id_);
}

inline bool operator!=(const gee_link_mode &a, const gee_link_mode &b) {
  return (a.id_ != b.id_);
}

struct gee_link_identity : public gee_link_mode {
  inline gee_link_identity() : gee_link_mode(1) {}
};

static const gee_link_identity identity_link;

struct gee_corr_mode {
  arma::uword id_;
  void setid(const arma::uword id) {id_ = id;}
  inline explicit gee_corr_mode(const arma::uword id) : id_(id) {}
};

inline bool operator==(const gee_corr_mode &a, const gee_corr_mode &b) {
  return (a.id_ == b.id_);
}

inline bool operator!=(const gee_corr_mode &a, const gee_corr_mode &b) {
  return (a.id_ != b.id_);
}

struct gee_corr_Identity : public gee_corr_mode {
  inline gee_corr_Identity() : gee_corr_mode(1) {}
};

struct gee_corr_CompSymm : public gee_corr_mode {
  inline gee_corr_CompSymm() : gee_corr_mode(2) {}
};

struct gee_corr_AR1 : public gee_corr_mode {
  inline gee_corr_AR1() : gee_corr_mode(3) {}
};

static const gee_corr_Identity Identity_corr;
static const gee_corr_CompSymm CompSymm_corr;
static const gee_corr_AR1 AR1_corr;

namespace gee {

  class gee_jmcm {
  public:
  gee_jmcm(const arma::uvec &m,
           const arma::vec &Y,
           const arma::mat &X,
           const arma::mat &Z,
           const arma::mat &W,
           const double rho,
           const gee_link_mode &link_mode,
           const gee_corr_mode &corr_mode):
    m_(m),Y_(Y),X_(X),Z_(Z),W_(W), rho_(rho),
      link_mode_(link_mode), corr_mode_(corr_mode){

      int debug = 0;

      int N = Y_.n_rows;
      int n_bta = X_.n_cols;
      int n_lmd = Z_.n_cols;
      int n_gma = W_.n_cols;

      H_ = arma::ones<arma::vec>(N);
      use_ipw_ = false;

      tht_ = arma::zeros<arma::vec>(n_bta + n_lmd + n_gma);
      bta_ = arma::zeros<arma::vec>(n_bta);
      lmd_ = arma::zeros<arma::vec>(n_lmd);
      gma_ = arma::zeros<arma::vec>(n_gma);

      Xbta_ = arma::zeros<arma::vec>(N);
      Zlmd_ = arma::zeros<arma::vec>(N);
      Wgma_ = arma::zeros<arma::vec>(W_.n_rows);
      Resid_ = arma::zeros<arma::vec>(N);

      free_param_ = 0;

      if (debug) Rcpp::Rcout << "gee_jmcm object created" << std::endl;
    }

    ~gee_jmcm(){};

    void set_free_param(const int n);

    void set_params(const arma::vec &x);
    arma::vec get_theta() const;
    void set_beta(const arma::vec &beta);
    arma::vec get_beta() const;
    void set_lambda(const arma::vec &lambda);
    arma::vec get_lambda() const;
    void set_gamma(const arma::vec &gamma);
    arma::vec get_gamma() const;
    void set_weights(const arma::vec &H);

    arma::vec operator()(const arma::vec &x);

    void UpdateGEES(const arma::vec &x);
    void UpdateParam(const arma::vec &x);
    void UpdateModel();

    void UpdateBeta();
    void UpdateLambda();
    void UpdateGamma();

    arma::uword get_m(const arma::uword i) const;
    arma::vec get_Y(const arma::uword i) const;
    arma::mat get_X(const arma::uword i) const;
    arma::mat get_Z(const arma::uword i) const;
    arma::mat get_W(const arma::uword i) const;

    arma::vec get_Resid(const arma::uword i) const;
    arma::mat get_D(const arma::uword i) const;
    arma::mat get_T(const arma::uword i) const;

    arma::vec get_mu(const arma::uword i) const;
    arma::mat get_Sigma(const arma::uword i) const;

    arma::mat get_Sigma_inv(const arma::uword i) const;
    arma::mat get_weights_sqrt(const arma::uword i) const;

    arma::mat get_fim() const; // get the fisher information matrix
    arma::vec get_sd() const;  // get the standand deviation of params

    bool learn(const arma::uvec &m, const arma::mat &Y, const arma::mat &X,
               const arma::mat &Z, const arma::mat &W,
               const gee_link_mode &link_mode, const gee_corr_mode &corr_mode,
               const double rho, const arma::vec &start,
               const arma::uword fs_iter, const bool print_mode = false);

    inline double get_quasi_likelihood(const arma::vec &x) {
      set_params(x);
      arma::uword n_sub = m_.n_elem;

      double result = 0.0;
      for(arma::uword i = 0; i < n_sub; ++i) {
        arma::mat Sigmai_inv = get_Sigma_inv(i);
        arma::vec ri = get_Resid(i);
        result += arma::as_scalar(ri.t() * Sigmai_inv * ri);
        //result += arma::as_scalar(ri.t() * ri);
      }
      result *= -0.5;

      return result;
    }

  private:
    arma::uvec m_;
    arma::vec Y_;
    arma::mat X_;
    arma::mat Z_;
    arma::mat W_;

    bool use_ipw_;
    arma::vec H_;

    arma::vec tht_;
    arma::vec bta_;
    arma::vec lmd_;
    arma::vec gma_;

    arma::vec Xbta_;
    arma::vec Zlmd_;
    arma::vec Wgma_;
    arma::vec Resid_;

    int free_param_;

    double rho_;
    gee_link_mode link_mode_;
    gee_corr_mode corr_mode_;

    bool fs_iterate(const gee_link_mode &link_mode,
                    const gee_corr_mode &corr_mode, const double rho,
                    const arma::vec &start, const arma::uword max_iter,
                    const bool verbose);
    void fs_update_params(const gee_link_mode &link_mode,
                          const gee_corr_mode &corr_mode, const double rho);
  };

  inline void gee_jmcm::set_free_param(const int n) { free_param_ = n; }

  inline void gee_jmcm::set_params(const arma::vec &x) {

    int fp2 = free_param_;
    free_param_ = 0;
    UpdateGEES(x);
    free_param_ = fp2;
  }

  inline arma::vec gee_jmcm::get_theta() const { return tht_; }

  inline void gee_jmcm::set_beta(const arma::vec& x) {
    int fp2 = free_param_;
    free_param_ = 1;
    UpdateGEES(x);
    free_param_ = fp2;
  }

  inline arma::vec gee_jmcm::get_beta() const { return bta_; }

  inline void gee_jmcm::set_lambda(const arma::vec& x) {
    int fp2 = free_param_;
    free_param_ = 2;
    UpdateGEES(x);
    free_param_ = fp2;
  }

  inline arma::vec gee_jmcm::get_lambda() const { return lmd_; }

  inline void gee_jmcm::set_gamma(const arma::vec& x) {
    int fp2 = free_param_;
    free_param_ = 3;
    UpdateGEES(x);
    free_param_ = fp2;
  }

  inline arma::vec gee_jmcm::get_gamma() const { return gma_; }

  inline void gee_jmcm::set_weights(const arma::vec &H) {
    use_ipw_ = true;
    H_ = H;
  }

  inline arma::uword gee_jmcm::get_m(const arma::uword i) const { return m_(i); }

  inline arma::vec gee_jmcm::get_Y(const arma::uword i) const {
    arma::vec Yi;
    if (i == 0)
      Yi = Y_.subvec(0, m_(0) - 1);
    else {
      arma::uword vindex = arma::sum(m_.subvec(0, i - 1));
      Yi = Y_.subvec(vindex, vindex + m_(i) - 1);
    }
    return Yi;
  }

  inline arma::mat gee_jmcm::get_X(const arma::uword i) const {
    arma::mat Xi;
    if (i == 0)
      Xi = X_.rows(0, m_(0) - 1);
    else {
      arma::uword rindex = arma::sum(m_.subvec(0, i - 1));
      Xi = X_.rows(rindex, rindex + m_(i) - 1);
    }
    return Xi;
  }

  inline arma::mat gee_jmcm::get_Z(const arma::uword i) const {
    arma::mat Zi;
    if (i == 0)
      Zi = Z_.rows(0, m_(0) - 1);
    else {
      arma::uword rindex = arma::sum(m_.subvec(0, i - 1));
      Zi = Z_.rows(rindex, rindex + m_(i) - 1);
    }
    return Zi;
  }

  inline arma::mat gee_jmcm::get_W(const arma::uword i) const {
    arma::mat Wi;
    if (m_(i) != 1) {
      if (i == 0) {
        arma::uword first_index = 0;
        arma::uword last_index = m_(0) * (m_(0) - 1) / 2 - 1;
        Wi = W_.rows(first_index, last_index);
      } else {
        arma::uword first_index = 0;
        for (arma::uword idx = 0; idx != i; ++idx) {
          first_index += m_(idx) * (m_(idx) - 1) / 2;
        }
        arma::uword last_index = first_index + m_(i) * (m_(i) - 1) / 2 - 1;
        Wi = W_.rows(first_index, last_index);
      }
    }
    return Wi;
  }

  inline arma::vec gee_jmcm::get_Resid(const arma::uword i) const {
    arma::vec ri;
    if (i == 0)
      ri = Resid_.subvec(0, m_(0) - 1);
    else {
      arma::uword vindex = arma::sum(m_.subvec(0, i - 1));
      ri = Resid_.subvec(vindex, vindex + m_(i) - 1);
    }
    return ri;
  }

  inline arma::mat gee_jmcm::get_D(const arma::uword i) const {
    arma::mat Di = arma::eye(m_(i), m_(i));
    if (i == 0)
      Di = arma::diagmat(arma::exp(Zlmd_.subvec(0, m_(0) - 1)));
    else {
      arma::uword rindex = arma::sum(m_.subvec(0, i - 1));
      Di = arma::diagmat(arma::exp(Zlmd_.subvec(rindex, rindex + m_(i) - 1)));
    }
    return Di;
  }

  inline arma::mat gee_jmcm::get_T(const arma::uword i) const {
    arma::mat Ti = arma::eye(m_(i), m_(i));
    if (m_(i) != 1) {
      if (i == 0) {
        arma::uword first_index = 0;
        arma::uword last_index = m_(0) * (m_(0) - 1) / 2 - 1;
        Ti = dragonwell::ltrimat(m_(0), -Wgma_.subvec(first_index, last_index));
      } else {
        arma::uword first_index = 0;
        for (arma::uword idx = 0; idx != i; ++idx) {
          first_index += m_(idx) * (m_(idx) - 1) / 2;
        }
        arma::uword last_index = first_index + m_(i) * (m_(i) - 1) / 2 - 1;
        Ti = dragonwell::ltrimat(m_(i), -Wgma_.subvec(first_index, last_index));
      }
    }
    return Ti;
  }

  inline arma::vec gee_jmcm::get_mu(const arma::uword i) const {
    arma::vec mui;
    if (i == 0)
      mui = Xbta_.subvec(0, m_(0) - 1);
    else {
      int index = arma::sum(m_.subvec(0, i - 1));
      mui = Xbta_.subvec(index, index + m_(i) - 1);
    }
    return mui;
  }

  inline arma::mat gee_jmcm::get_Sigma(const arma::uword i) const {
    int debug = 0;

    arma::mat Ti = get_T(i);
    arma::mat Ti_inv = arma::pinv(Ti);
    arma::mat Di = get_D(i);

    if (debug) {
      Ti.print("Ti = ");
      Di.print("Di = ");
    }
    return Ti_inv * Di * Ti_inv.t();
  }

  inline arma::mat gee_jmcm::get_Sigma_inv(const arma::uword i) const {
    arma::mat Ti = get_T(i);
    arma::mat Di = get_D(i);
    arma::mat Di_inv = arma::diagmat(arma::pow(Di.diag(), -1));
    return Ti.t() * Di_inv * Ti;
  }

  inline arma::mat gee_jmcm::get_weights_sqrt(const arma::uword i) const {

    int debug = 0;

    arma::mat result = arma::eye(m_(i), m_(i));
    if (i == 0) {
      if (debug) Rcpp::Rcout << "i = " << i << std::endl;
      arma::vec Hi = H_.subvec(0, m_(0) - 1);
      result = arma::diagmat(arma::sqrt(Hi));
      //result = arma::diagmat(Hi);
      if (debug) result.print("result = ");
    } else {
      if (debug) Rcpp::Rcout << "i = " << i << std::endl;
      arma::uword vindex = arma::sum(m_.subvec(0, i - 1));
      arma::vec Hi = H_.subvec(vindex, vindex + m_(i) - 1);
      result = arma::diagmat(arma::sqrt(Hi));
      //result = arma::diagmat(Hi);
      if (debug) result.print("result = ");
    }
    return result;
  }

  inline arma::mat gee_jmcm::get_fim() const {
    arma::uword n_sub = m_.n_elem;
    arma::uword n_bta = X_.n_cols;
    arma::uword n_lmd = Z_.n_cols;
    arma::uword n_gma = W_.n_cols;
    arma::uword n_tht = n_bta + n_lmd + n_gma;

    arma::mat fim     = arma::zeros<arma::mat>(n_tht, n_tht);
    arma::mat fim_bta = arma::zeros<arma::mat>(n_bta, n_bta);
    arma::mat fim_lmd = arma::zeros<arma::mat>(n_lmd, n_lmd);
    arma::mat fim_gma = arma::zeros<arma::mat>(n_gma, n_gma);
    for (arma::uword i = 0; i != n_sub; ++i) {

      // Compute fisher information matrix(beta)
      // arma::vec ri = get_Resid(i);
      arma::mat Sigmai_inv = get_Sigma_inv(i);
      if (use_ipw_) {
        arma::mat H_sqrt = get_weights_sqrt(i);
        Sigmai_inv = H_sqrt * Sigmai_inv * H_sqrt;
      }

      arma::mat deriv1;
      // arma::vec yi_tilde;
      if (link_mode_ == identity_link) {
        deriv1 = get_X(i).t();
        // yi_tilde = get_Y(i);
      }
      fim_bta += deriv1 * Sigmai_inv * deriv1.t();

      // Compute fisher information matrix(lambda)
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
        cov_inv = H_sqrt * cov_inv * H_sqrt;
      }

      arma::mat deriv2 = get_Z(i).t();
      for (arma::uword t = 1; t <= mi; ++t) {
        deriv2.col(t - 1) *= arma::as_scalar(di(t - 1));
      }

      arma::vec epsi2 = arma::pow(Ti * ri, 2);
      arma::vec epsi2_tilde = epsi2 - di2 + Di * logdi2;

      fim_lmd += deriv2 * cov_inv * deriv2.t();

      // Compute fisher information matrix(gamma)
      //arma::uword mi = m_(i);
      arma::mat Wi = get_W(i);
      //arma::vec ri = get_Resid(i);
      //arma::mat Di = get_D(i);
      //arma::mat Di_inv = arma::diagmat(arma::pow(Di.diag(), -1));
      if (use_ipw_) {
        arma::mat H_sqrt = get_weights_sqrt(i);
        Di_inv = H_sqrt * Di_inv * H_sqrt;
      }

      arma::uword rindex = 0;

      arma::mat deriv3 = arma::zeros<arma::mat>(n_gma, mi);
      for (arma::uword j = 2; j <= mi; ++j) {
        for (arma::uword k = 1; k <= (j - 1); ++k) {
          deriv3.col(j - 1) += ri(k - 1) * Wi.row(rindex).t();
          ++rindex;
        }
      }

      fim_gma += deriv3 * Di_inv * deriv3.t();
    }

    fim.submat(0, 0, n_bta-1, n_bta-1) = fim_bta;
    fim.submat(n_bta, n_bta, n_bta+n_lmd-1, n_bta+n_lmd-1) = fim_lmd;
    fim.submat(n_bta+n_lmd, n_bta+n_lmd, n_bta+n_lmd+n_gma-1, n_bta+n_lmd+n_gma-1) = fim_gma;

    return fim;
  }

  inline arma::vec gee_jmcm::get_sd() const {
    arma::uword n_bta = X_.n_cols;
    arma::uword n_lmd = Z_.n_cols;
    arma::uword n_gma = W_.n_cols;

    arma::mat fim = get_fim();
    arma::mat fim_bta = fim.submat(0, 0, n_bta-1, n_bta-1);
    arma::mat fim_lmd = fim.submat(n_bta, n_bta, n_bta+n_lmd-1, n_bta+n_lmd-1);
    arma::mat fim_gma = fim.submat(n_bta+n_lmd, n_bta+n_lmd,
                                   n_bta+n_lmd+n_gma-1, n_bta+n_lmd+n_gma-1);

    arma::vec bta_sd = arma::sqrt(arma::diagvec(fim_bta.i()));
    arma::vec lmd_sd = arma::sqrt(arma::diagvec(fim_lmd.i()));
    arma::vec gma_sd = arma::sqrt(arma::diagvec(fim_gma.i()));

    arma::vec sd = dragonwell::join_vecs({bta_sd, lmd_sd, gma_sd});
    return sd;
  }
}  // namespace gee4

#endif  // GEE_JMCM_H_
