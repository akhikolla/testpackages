#ifndef GAUSSLOGLIKH
#define GAUSSLOGLIKH

#include<vector>
#include <RcppArmadillo.h>

const double LN2PI = log(2*PI);

void MDataGaussLogLik(const int n, const int p, const int Cf, const arma::mat& X, const arma::vec& u, 
                      arma::mat* Sigmap, arma::mat* SigmaInvp, double* lndetSigp, std::vector<double>& res, bool& validsol, 
//                      const double maxlnk2, const double MaxSctlnEgvlRt, const bool chksing=true);
                      const double maxlnk2, const bool chksing=true);

#endif
