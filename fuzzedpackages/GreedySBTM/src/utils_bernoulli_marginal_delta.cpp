#include "utils_bernoulli_marginal_delta.h"

double BernoulliMarginalDelta(double N2, double eta2, double N1, double eta1, double a, double b)
{
    double res = 0;
    if (N1 > 0 && N2 > 0) res += lgamma(a+eta2) - lgamma(a+eta1) + lgamma(b+N2-eta2) - lgamma(b+N1-eta1) - lgamma(a+b+N2) + lgamma(a+b+N1);
    else if (N1 == 0) res += lgamma(a+b) - lgamma(a) - lgamma(b) + lgamma(a+eta2) + lgamma(b+N2-eta2) - lgamma(a+b+N2);
    else if (N2 == 0) res -= lgamma(a+b) - lgamma(a) - lgamma(b) + lgamma(a+eta1) + lgamma(b+N1-eta1) - lgamma(a+b+N1);
    return(res);
}


