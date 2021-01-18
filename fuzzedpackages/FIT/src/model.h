// model.h
#ifndef MODEL_H_
#define MODEL_H_

#define _USE_MATH_DEFINES
#include <cmath>
double constexpr pi_ = M_PI; // 3.14159265358979323846;
#undef _USE_MATH_DEFINES

namespace model {

////////////////////////////////////////////////////////////////
// Functions defining the model

int constexpr periodC = 1440;
double constexpr _2pi_over_periodC = 2 * pi_ / periodC;

inline Rcpp::NumericVector clockCos(Rcpp::NumericVector t){
  return cos(_2pi_over_periodC * t) * std::sqrt(2);
}
inline Rcpp::NumericVector clockSin(Rcpp::NumericVector t){
  return sin(_2pi_over_periodC * t) * std::sqrt(2);
}

//inline double clockC(int const t, double const phase) {
//  return 0.5 * std::cos(_2pi_over_periodC * (t-phase));
//}

inline double envF(double const w, double const envAmp, double const envTh) {
  double x = std::exp(envAmp) * (w - envTh);
  return (x > 0) ? std::tanh(x)*std::sqrt(std::exp(-2*envAmp)+1) : 0;
}

// - exploit the fact that t and gatePhase appears in the combination (t-gatePh)
// - must be 1440 periodic in t nad (t-gatePh)
// - val between 0 and 1 (must be non-negative)
int constexpr periodG = 1440;
double constexpr _2pi_over_periodG = 2 * pi_ / periodG;

inline double gateG(int const t, double const gatePh,
                    double const gateAmp, double const gateTh) {

  double gateClock = std::cos(_2pi_over_periodG * (t-gatePh)) - gateTh;
  if (gateClock == 0) return 0.5;
  double expGateAmp = std::exp(gateAmp);
  if (expGateAmp == 0) return 0.0;

  return (std::tanh(expGateAmp * gateClock) - std::tanh(expGateAmp*(-1-gateTh)))
                / (std::tanh(expGateAmp*(1-gateTh)) - std::tanh(expGateAmp*(-1-gateTh)));
}

////////////////////////////////////////////////////////////////
} // namespace model
#endif // MODEL_H_
