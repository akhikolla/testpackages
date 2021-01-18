#ifndef INTRAPOLATE_H
#define INTRAPOLATE_H
#include <cmath>
#include <limits>
#include <algorithm>

namespace PSQN {
/***
 performs intrapolation.
 */
class intrapolate {
  double const f0, d0;
  double xold = std::numeric_limits<double>::quiet_NaN(),
         fold = std::numeric_limits<double>::quiet_NaN(),
         xnew, fnew;
  bool has_two_values = false;

public:
  intrapolate(double const f0, double const d0, double const x,
              double const f): f0(f0), d0(d0), xnew(x), fnew(f) { }

  inline double get_value(double const v1, double const v2)
  const noexcept {
    double const a = std::min(v1, v2),
                 b = std::max(v1, v2),
             small = .01 * (b - a);

    double const val = ([&]() -> double {
      if(!has_two_values)
        return - d0 * xnew * xnew / 2 / (fnew - f0 - d0 * xnew);

      double fac = xnew * xnew * xold * xold * (xnew - xold);
      if(fac == 0)
        fac = 1;
      double const f1 = fnew - f0 - d0 * xnew,
                   f2 = fold - f0 - d0 * xold,
                   vb = (-xold * xold * xold  * f1 +
                     xnew * xnew * xnew * f2) / fac,
                   va = (xold * xold * f1 - xnew * xnew * f2) / fac,
                deter = vb * vb - 3. * va * d0;
      if(deter < 0)
        return a + (b - a) / 2.;

      return (-vb + std::sqrt(deter)) / 3. / va;
    })();

    if(val < a + small or val > b - small)
      return a + (b - a) / 2.;

    return val;
  }

  inline void update(double const x, double const f) noexcept {
    xold = xnew;
    fold = fnew;
    xnew = x;
    fnew = f;
    has_two_values = true;
  }
};
} // namespace PSQN

#endif
