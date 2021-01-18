#ifndef IPW_H_
#define IPW_H_

#include <RcppArmadillo.h>

namespace gee {

  class ipw {
  public:
  ipw(const arma::uvec &m, const arma::vec &Y, arma::uword order)
    : m_(m), Y_(Y), order_(order) {
      arma::uword ntobs = Y_.n_elem;
      p_ = arma::zeros<arma::vec>(ntobs);
      Pi_ = arma::zeros<arma::vec>(ntobs);
    }

    arma::vec get_Y(const arma::uword i) const {
      arma::vec Yi;
      if (i == 0)
        Yi = Y_.subvec(0, m_(0) - 1);
      else {
        arma::uword vindex = arma::sum(m_.subvec(0, i - 1));
        Yi = Y_.subvec(vindex, vindex + m_(i) - 1);
      }
      
      return Yi;
    }

    arma::vec get_Z(const arma::uword i, const arma::uword j) const {
      arma::vec Yi = get_Y(i);
      arma::vec Z_ij = arma::zeros<arma::vec>(order_+1);
      if ((j+1) <= order_) {
        Z_ij(0) = 1;
      }
      else {
        Z_ij(0) = 1;
        for (arma::uword k = 1; k <= order_; ++k) Z_ij(k) = Yi(j - k); 
      }

      return Z_ij;
    }
    
    arma::vec operator()(const arma::vec &alpha) {
      arma::uword max_obs = m_.max();
      arma::uword nsub = m_.n_elem;
      arma::uword npar = alpha.n_elem;

      arma::vec result = arma::zeros<arma::vec>(npar);
      for (arma::uword i = 0; i != nsub; ++i) {
        for (arma::uword j = 1; j != m_(i); ++j) {
          arma::vec Z_ij = get_Z(i, j);          
          double tmp = arma::as_scalar(arma::exp(Z_ij.t() * alpha));
          double p_ij =  tmp / (1 + tmp);
          result += 1 * (1 - p_ij) * Z_ij;
        }
        if (m_(i) != max_obs) {
          arma::vec Z_ij = get_Z(i, m_(i));
          double tmp = arma::as_scalar(arma::exp(Z_ij.t() * alpha));
          double p_ij =  tmp / (1 + tmp);
          result += 1 * (0 - p_ij) * Z_ij;
        }
      }

      return result;
    }

    arma::vec CalWeights(const arma::vec &alpha) {
      int debug = 0;
      
      arma::uword nsub  = m_.n_elem;
      arma::uword index = 0, pindex = 0, Piindex = 0;
      arma::vec result = arma::zeros<arma::vec>(Y_.n_elem);
      for (arma::uword i = 0; i != nsub; ++i) {
        arma::vec Pi_i = arma::zeros<arma::vec>(m_(i));
        for (arma::uword j = 0; j != m_(i); ++j, ++index) {
          if (j == 0) {
            double p_ij = 1;
            Pi_i(0) = p_ij;
            result(index) = 1 / Pi_i(0);

            p_(pindex++) = p_ij;
          } else {
            arma::vec Z_ij = get_Z(i, j);            
            double tmp = arma::as_scalar(arma::exp(Z_ij.t() * alpha));
            double p_ij =  tmp / (1 + tmp);
            Pi_i(j) = Pi_i(j-1) * p_ij;

            if (debug && i == (nsub - 1)) {
              Rcpp::Rcout << "j = " << j+1 << std::endl;
              Rcpp::Rcout << "Pi_ij = " << Pi_i(j) << std::endl;
              Rcpp::Rcout << std::endl;
            }
	    
            /* if (Pi_i(j) < 1e-7) Pi_i(j) = 1e-7;   */
	    
            result(index) = 1 / Pi_i(j);

            p_(pindex++) = p_ij;
          }
        } // for loop j

        Pi_.subvec(Piindex, Piindex + m_(i) - 1) = Pi_i;
        Piindex += m_(i);
	
      } // for loop i
      return result;
    } 

    arma::vec get_p() const { return p_; }
    arma::vec get_Pi() const { return Pi_; }
    
  private:
    arma::uvec m_;
    arma::vec Y_;
    arma::uword order_;

    arma::vec p_, Pi_;
  };
}

#endif
