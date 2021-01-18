#ifndef PSQN_REPORTER_H
#define PSQN_REPORTER_H

#include <Rcpp.h>
#include <ios>
#include <iostream>
#include <iomanip>

namespace PSQN {
/** class used to print to the console during estimation */
struct R_reporter {
  static void cg_it(int const trace, size_t const iteration,
                    size_t const maxit, double const r_norm,
                    double const threshold) {
    if(trace > 3L and iteration % (maxit / 5L) == 0L)
      Rcpp::Rcout << "      Conjugate gradient iteration " << iteration
                  << ". Residual norm is " << r_norm
                  << " (threshold is " << threshold << ")\n";
  }

  static void cg(int const trace, size_t const iteration,
                 size_t const n_cg, bool const successful) {
    if(trace > 0)
      Rcpp::Rcout << "Conjugate gradient "
                  << (successful ? "succeeded" : "failed")
                  << " in itteration " << iteration + 1L
                  << '\n';

    if(trace > 2L)
      Rcpp::Rcout << "    " << n_cg
                  << " conjugate itterations have been used\n";
  }

  static void line_search
  (int const trace, size_t const iteration, size_t const n_eval,
   size_t const n_grad, double const fval_old,
   double const fval, bool const successful, double const step_size,
   double const *new_x, size_t const n_global) {
    if(trace > 0)
      Rcpp::Rcout << "Line search "
                  << (successful ? "succeeded" : "failed")
                  << '\n';

    if(trace > 1L){
      std::streamsize n_digits(9 - std::log10(fval_old)),
                      old_size = Rcpp::Rcout.precision();
      Rcpp::Rcout << "  New (old) function value is "
                  << std::fixed
                  << std::setprecision(n_digits)
                  << fval << " (" << fval_old << ")\n";

      // see https://stackoverflow.com/a/35173760/5861244
      Rcpp::Rcout.unsetf(std::ios_base::floatfield);
      Rcpp::Rcout << std::setprecision(old_size);
    }

    if(trace > 2L){
      Rcpp::Rcout << "    step size is " << step_size
                  << " and new global parameters are\n      ";
      for(size_t i = 0; i < n_global; ++i)
        Rcpp::Rcout << new_x[i] << " ";

      Rcpp::Rcout << "\n    " << n_eval
                  << " function evaluations and "
                  << n_grad << " gradient evaluations have been used\n";
    }

    if(trace > 0L)
      Rcpp::Rcout << '\n';
  }

  static void line_search_inner
  (int const trace, double const a_old, double const a_new,
   double const f_new, bool const is_zoom, double const d_new,
   double const a_high) {
    if(trace < 4L)
      return;

    if(is_zoom)
      Rcpp::Rcout << "      (zoom) lower: " << a_old
                  << " upper: " << a_high
                  << " new value: " << a_new
                  << " f new: " << f_new
                  << " d new: " << d_new << '\n';
    else
      Rcpp::Rcout << "      a_prev: " << a_old
                  << " new value: " << a_new
                  << " f new: " << f_new
                  << " d new: " << d_new << '\n';
  }
};

class R_interrupter {
public:
  static bool check_interrupt() {
    try {
      Rcpp::checkUserInterrupt();
    }
    catch(Rcpp::internal::InterruptedException e)
    {
      return true;
    }
    return false;
  }
};
} // namespace PSQN

#endif
