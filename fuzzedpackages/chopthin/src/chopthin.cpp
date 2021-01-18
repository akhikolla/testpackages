#include <Rcpp.h>

using namespace Rcpp;

int intrand(const int n) { return floor(unif_rand()*n); }
inline double myrunif(){return unif_rand();}
#define chopthin_error(x) throw Rcpp::exception((x))
#include "chopthin_internal.h"

/////
//' The Chopthin Resampler
//'
//' A fast implementation of the Chopthin resampler. Can be used as the
//' resampling step in particle filters and in sequential Monte Carlo.
//'
//' @param w a vector of weights
//' @param N target number of particles
//' @param eta upper bound on the ratio between the weights. Must be >=4.
//'        If eta=Inf then only thinnig is performed, requiring the number of weights to be at least N.
//' @param normalise Flag for controlling if the returned weights should
//'        be normalised. If TRUE (default) then the sum of the returned
//'        weights will sum to N. If FALSE then the returned weights have
//'        the same sum as the original weights (within the numerical precision).
//' @param checks Flag controlling if checks on the input and the result should be performed. Default TRUE.
//' @return A list with two elements: new weights and indices
//' of the ancestors of the new particles. The weights are normalised to add up to N.
//'
//' @references
//' A Gandy and F. D-H Lau. The chopthin algorithm for
//' resampling.  IEEE Transactions on Signal Processing, 64(16):4273--4281, 2016
//'
//' @examples
//' chopthin(runif(10),N=10)
//' chopthin(runif(10),N=20,4)
//' chopthin(runif(10),N=5)
//' chopthin(runif(10),N=1)
//' @export
// [[Rcpp::export]]
List chopthin(std::vector<double>& w, int N, double eta=5.828427, bool normalise=true, bool checks=true){
  std::vector<double> wres(N);
  std::vector<int> ires(N);
  if (checks)
    chopthin_internal<true>(w,N,wres,ires,eta,normalise);
  else
    chopthin_internal<false>(w,N,wres,ires,eta,normalise);
  List res;
  res["weights"]=NumericVector(wres.begin(),wres.end());
  res["indices"]=IntegerVector(ires.begin(),ires.end());
  return res;
}

