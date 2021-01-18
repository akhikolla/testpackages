#include <Rcpp.h>
#include <random>
using namespace Rcpp;

// Function to seed and generate samples from a normal distribution

// [[Rcpp::export]]
NumericVector rng(int n, unsigned seed)
{
  // Vector to return
  NumericVector v(n);

  // Set seed and generator
  std::default_random_engine generator;
  generator.seed(seed);
  std::normal_distribution<double> distribution(0.0,1.0);

  // Generate values
  for (int i=0; i<n; ++i) {
    v[i] = distribution(generator);
  }

  return v;
}
