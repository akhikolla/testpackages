#include <Rcpp.h>

inline void eval_marker_inner
  (double * __restrict__ o, double const * const __restrict__ m_start,
   double const * __restrict__ b, size_t const lm,
   size_t const lo) noexcept {
  for(size_t i = 0; i < lo; ++i, ++o){
    double const * mj = m_start;
    for(size_t j = 0; j < lm; ++j, ++mj, ++b)
      *o += *mj * *b;
  }
}

// [[Rcpp::export(name = "eval_marker_cpp", rng = false)]]
void eval_marker(SEXP B, SEXP m, SEXP Sout){
  bool const out_is_mat = Rf_isMatrix(Sout);
  if(Rf_isMatrix(B) and out_is_mat){
    if(Rf_isMatrix(m)){
      size_t const nr = Rf_nrows(B),
                   nc = Rf_ncols(B),
            n_col_out = Rf_nrows(m),
                   nm = Rf_ncols(m);

      bool const B_m_ok = nr == nm,
                 out_ok = static_cast<size_t>(Rf_ncols(Sout)) == n_col_out;
      if(B_m_ok and out_ok){
        double * o = REAL(Sout);
        double const * const m_start = REAL(m),
                     * const b_start = REAL(B);

        for(size_t i = 0; i < n_col_out; ++i){
          double const * const o_end = o + nc,
                       *           b = b_start;
          for(; o != o_end; ++o){
            double const *           mi = m_start + i,
                         * const mi_end = mi + n_col_out * nm;

            for(; mi != mi_end; mi += n_col_out, ++b)
              *o += *mi * *b;
          }
        }

        return;

      } else
        throw std::invalid_argument("eval_marker: dims do not match");

    } else if(Rf_isVector(m)){
      size_t const nr = Rf_nrows(B),
                   nc = Rf_ncols(B),
                   nm = XLENGTH(m);

      bool const B_m_ok = nr == nm,
                 out_ok = Rf_ncols(Sout) == 1L;
      if(B_m_ok and out_ok){
        double const *b = REAL(B),
               *m_start = REAL(m);

        eval_marker_inner(REAL(Sout), m_start, b, nm, nc);

        return;

      } else
        throw std::invalid_argument("eval_marker: dims do not match");

    }
  }

  throw std::invalid_argument("eval_marker: B and Sout must be a matrix. m must be a vector or a Matrix");
}
