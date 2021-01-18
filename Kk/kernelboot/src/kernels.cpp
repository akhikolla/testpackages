
#include "kernels.h"
#include "shared.h"

// Constants

const double SQRT_3 = 1.732050807568877193177;
const double SQRT_5 = 2.236067977499789805051;
const double SQRT_6 = 2.449489742783177881336;
const double SQRT_7 = 2.645751311064590716171;
const double PI_SQ = M_PI * M_PI;
const double COS_C = 1.0 / std::sqrt(1.0/3.0 - 2.0/PI_SQ);
const double OPTCOS_C = 1.0 / std::sqrt(1.0 - 8.0/PI_SQ);

// kernel RNGs

double rng_epan() {
  double u1, u2, u3, au1, au2, au3;
  u1 = (rng_unif() * 2.0 - 1.0);
  u2 = (rng_unif() * 2.0 - 1.0);
  u3 = (rng_unif() * 2.0 - 1.0);
  au1 = std::abs(u1);
  au2 = std::abs(u2);
  au3 = std::abs(u3);

  if (au3 >= au2 && au3 >= au1) {
    return u2 * SQRT_5;
  } else {
    return u3 * SQRT_5;
  }
}

double rng_cosine() {
  return (R::rbeta(3.3575, 3.3575) * 2.0 - 1.0) * COS_C;
}

double rng_optcos() {
  double u;
  u = (rng_unif() * 2.0 - 1.0);
  return (2.0 * std::acos(u)/M_PI - 1.0) * OPTCOS_C;
}

double rng_triang() {
  double u, v;
  u = rng_unif();
  v = rng_unif();
  return (u-v) * SQRT_6;
}

double rng_rect() {
  return (rng_unif() * 2.0 - 1.0) * SQRT_3;
}

double rng_biweight() {
  return (R::rbeta(3.0, 3.0) * 2.0 - 1.0) * SQRT_7;
}

