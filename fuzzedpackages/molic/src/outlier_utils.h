#ifndef OUTLIER_UTILS_H
#define OUTLIER_UTILS_H

#include "molic_types.h"

VS     matpr(Rcpp::CharacterMatrix A);
RIV    n_a(RCM & A);
RIV    n_b(RIV & na, RIV & b);
VD     subtract_one(VD x);
VD     Gx_(VD x);
VD     Hx_(VD x);
RCM    subM( RCM & A, RCV & x );
RL     a_marginals( RCM A, RL & am );
double TY(RCV y, RL & C_marginals, RL & S_marginals);

#endif
