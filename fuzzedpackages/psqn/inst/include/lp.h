#ifndef LPSQN_P_H
#define LPSQN_P_H
#include <cstddef>
#include <algorithm>
#include "constant.h"

namespace lp {

inline void copy(double * PSQN_RESTRICT x, double const * PSQN_RESTRICT y,
                 size_t const dim) noexcept {
  for(size_t i = 0; i < dim; ++i, ++x, ++y)
    *x = *y;
}

inline void vec_diff
(double const * PSQN_RESTRICT x, double const * PSQN_RESTRICT y,
 double * PSQN_RESTRICT res, size_t const n) noexcept {
  for(size_t i = 0; i < n; ++i, ++x, ++y, ++res)
    *res = *x - *y;
}

inline double vec_dot(double const *x, size_t const n) noexcept {
  double out(0.);
  for(size_t i = 0; i < n; ++i, ++x)
    out += *x * *x;
  return out;
}

inline double vec_dot
(double const * PSQN_RESTRICT x, double const * PSQN_RESTRICT y,
 size_t const n) noexcept {
  double out(0.);
  for(size_t i = 0; i < n; ++i, ++x, ++y)
    out += *x * *y;
  return out;
}

/**
 computes b <- b + Xx where is is a n x n symmetric matrix containing only
 the upper triangular and x is a n-dimensional vector.
 */
inline void mat_vec_dot
(double const * PSQN_RESTRICT  X, double const * const PSQN_RESTRICT x,
 double * const PSQN_RESTRICT res, size_t const n) noexcept {
  double const * xj = x;
  double * rj = res;
  for(size_t j = 0; j < n; ++j, ++xj, ++rj){
    double const *xi = x;
    double * ri = res;

    for(size_t i = 0L; i < j; ++i, ++X, ++ri, ++xi){
      *ri += *X * *xj;
      *rj += *X * *xi;
    }
    *rj += *X++ * *xi;
  }
}

/**
 computes b <- b + Xx where b and x are seperated into an nb1 and bn2
 dimensional vector.
 */
inline void mat_vec_dot
(double const * PSQN_RESTRICT X, double const * PSQN_RESTRICT x1,
 double const * PSQN_RESTRICT x2, double * const PSQN_RESTRICT r1,
 double * const PSQN_RESTRICT r2,
 size_t const n1, size_t const n2) noexcept {
  size_t const n = n1 + n2;
  auto loop_body =
    [&](double const xj, double &rj, size_t const j) -> void {
      size_t i = 0L;
      {
        double       * ri = r1;
        double const * xi = x1;
        size_t const iend = std::min(j, n1);
        for(; i < iend; ++i, ++X, ++ri, ++xi){
          *ri += *X *  xj;
           rj += *X * *xi;
        }
        if(i < n1)
          rj += *X++ * *xi;
      }

      if(i == n1){ // still work to do
        double       * ri = r2;
        double const * xi = x2;
        for(; i < j; ++i, ++X, ++ri, ++xi){
          *ri += *X *  xj;
           rj += *X * *xi;
        }
        rj += *X++ * *xi;

      }
  };

  {
    double const *xj = x1;
    double       *rj = r1;
    for(size_t j = 0; j < n1; ++j, ++xj, ++rj)
      loop_body(*xj, *rj, j);
  }

  {
    double const *xj = x2;
    double       *rj = r2;
    for(size_t j = n1; j < n; ++j, ++xj, ++rj)
      loop_body(*xj, *rj, j);
  }
}

/**
 computes b <- b + Xx where b and x are seperated into an nb1 and bn2
 dimensional vector but excluding the first n1 x n1 block of X.
 */
inline void mat_vec_dot_excl_first
(double const * PSQN_RESTRICT X, double const * PSQN_RESTRICT x1,
 double const * PSQN_RESTRICT x2, double * const PSQN_RESTRICT r1,
 double * const PSQN_RESTRICT r2,
 size_t const n1, size_t const n2) noexcept {
  size_t const n = n1 + n2;
  auto loop_body =
    [&](double const xj, double &rj, size_t const j,
        bool const excl) -> void {
      size_t const end_first = std::min(j, n1);
      size_t i = excl ? end_first : 0L;
      if(excl)
        X += end_first + (j < n1);
      else {
        double       * ri = r1;
        double const * xi = x1;
        size_t const iend = end_first;
        for(; i < iend; ++i, ++X, ++ri, ++xi){
          *ri += *X *  xj;
           rj += *X * *xi;
        }
        if(i < n1)
          rj += *X++ * *xi;
      }

      if(i == n1){ // still work to do
        double       * ri = r2;
        double const * xi = x2;
        for(; i < j; ++i, ++X, ++ri, ++xi){
          *ri += *X *  xj;
           rj += *X * *xi;
        }
        rj += *X++ * *xi;

      }
  };

  {
    double const *xj = x1;
    double       *rj = r1;
    for(size_t j = 0; j < n1; ++j, ++xj, ++rj)
      loop_body(*xj, *rj, j, true);
  }

  {
    double const *xj = x2;
    double       *rj = r2;
    for(size_t j = n1; j < n; ++j, ++xj, ++rj)
      loop_body(*xj, *rj, j, false);
  }
}

/**
 performs a rank one update X <- X + scal * x.x^T where X is a symmetric
 matrix contaning only the upper triangular.
 */
inline void rank_one_update
(double * PSQN_RESTRICT X, double const * PSQN_RESTRICT x,
 double const scal, size_t const n)
noexcept {
  double const * xj = x;
  for(size_t j = 0; j < n; ++j, ++xj){
    double const * xi = x;
    for(size_t i = 0; i <= j; ++i, ++xi, ++X)
      *X += scal * *xj * *xi;
  }
}

/***
 performs the update
   H <- (I - cx.y^T).H.(I - cy.x^T) + cx.x^T
      = H - cx.y^T.H - cH.y.x^T + c^2x.(y^T.H.y).x^T + cx.x^T
 where X is a symmetric matrix contaning only the upper triangular.
 */
inline void bfgs_update
  (double * PSQN_RESTRICT X, double const * PSQN_RESTRICT x,
   double const * PSQN_RESTRICT H_y, double const y_H_y,
   double const scal, size_t const n)
  noexcept {
  double const * xj = x,
             * H_yj = H_y;
  double const y_H_y_p_scal = scal * (scal * y_H_y + 1);
  for(size_t j = 0; j < n; ++j, ++xj, ++H_yj){
    double const * xi = x,
               * H_yi = H_y;
    for(size_t i = 0; i <= j; ++i, ++xi, ++H_yi, ++X){
      *X += y_H_y_p_scal * *xj * *xi;
      *X -= scal * (*H_yj * *xi + *H_yi * *xj);
    }
  }
}
} // namespace lp

#endif
