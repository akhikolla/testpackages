#ifndef NEWTON_H_
#define NEWTON_H_

#include "linesearch.h"

#include <RcppArmadillo.h>
#include <limits>

namespace dragonwell {
  template <typename T>
    struct NRfmin {
      arma::vec fvec;
      T &func;
    NRfmin(T &funcc) : func(funcc) {}
      double operator()(const arma::vec &x) {
        /* arma::uword n = x.n_elem; */
        fvec = func(x);
        double sum = arma::as_scalar(fvec.t() * fvec);
        return 0.5 * sum;
      }
    };

  template <typename T>
    struct NRfdjac {
      const double kEpsilon;
      T &func;
    NRfdjac(T &funcc) : kEpsilon(1.0e-8), func(funcc) {}
      arma::mat operator()(const arma::vec &x, arma::vec &fvec) {
        arma::uword n = x.n_elem;
        arma::mat df = arma::zeros<arma::mat>(n, n);
        arma::vec xh = x;
        for (arma::uword j = 0; j < n; ++j) {
          double temp = xh(j);
          double h = kEpsilon * std::abs(x(j));
          if (h == 0.0) h = kEpsilon;
          xh(j) = temp + h;
          h = xh(j) - temp;
          arma::vec f = func(xh);
          xh(j) = temp;
          df.col(j) = (f - fvec) / h;
        }
        return df;
      }
    };

  template <typename T>
    class Newton : public LineSearch<NRfmin<T>> {
  public:
  Newton(T &func) : func_(func) {n_iter_ = 0; f_min_ = 0.0;}
    ~Newton(){};

    bool Optimize(arma::vec &x, const double grad_tol = 1e-6,
                  const bool print_mode = false) {
      bool check;
      const arma::uword kMaxIters = 2000;
      const double kTolF = 1.0e-8;
      const double kTolMin = 1.0e-12;
      const double kStpMax = 100.0;
      //const double kTolX = std::numeric_limits<double>::epsilon();
      const double kTolX = 1.0e-10;

      const arma::uword n = x.n_elem;
      NRfmin<T> fmin(func_);
      NRfdjac<T> fdjac(func_);
      arma::vec &fvec = fmin.fvec;
      double f = fmin(x);

      double test = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        if (std::abs(fvec(i)) > test) test = std::abs(fvec(i));
      }
      if (test < 0.01 * kTolF) {
        check = false;
        return check;
      }
      double sum = arma::as_scalar(x.t() * x);
      double stpmax = kStpMax * std::max(std::sqrt(sum), (double)n);

      // double fold = 0.0;
      arma::mat fjac = arma::zeros<arma::mat>(n, n);
      for (arma::uword iter = 0; iter < kMaxIters; ++iter) {
        fjac = fdjac(x, fvec);
        arma::vec g = fjac.t() * fvec;
        arma::vec xold = x;
        // fold = f;
        arma::vec p = -fvec;

        //p = arma::solve(fjac, p);
        arma::vec ptmp;
        if (arma::solve(ptmp, fjac, p)) p = ptmp;
        else p = arma::pinv(fjac.t() * fjac) * fjac.t() * p;
        check = this->GetStep(fmin, &f, &x, g, p, stpmax);

        if (print_mode) {
          Rcpp::Rcout << iter << ": " << f << ": " << x.t() << std::endl;
        }

        test = 0.0;
        for (arma::uword i = 0; i < n; ++i) {
          if (std::abs(fvec(i)) > test) test = std::abs(fvec(i));
        }
        if (test < 0.01 * kTolF) {
          check = false;
          return check;
        }
        if (check) {
          test = 0.0;
          double den = std::max(f, 0.5*n);
          for(arma::uword i = 0; i < n; ++i) {
            double temp = std::abs(g(i)) * std::max(std::abs(x(i)),1.0)/den;
            if(temp > test) test = temp;
          }
          check = (test < kTolMin);
          return check;
        }

        test = 0.0;
        for (arma::uword i = 0; i < n; ++i) {
          double temp = std::abs(x(i) - xold(i)) / std::max(std::abs(x(i)), 1.0);
          if (temp > test) test = temp;
        }
        if (test < kTolX) return check;
      }
      return check;
    }

    arma::uword n_iters() {return n_iter_;}
    double f_min() {return f_min_;}
  private:
    T func_;
    arma::uword n_iter_;
    double f_min_;
  };

}  // namespace dragon

#endif  // ROOTS_MULTIDIM_H_
