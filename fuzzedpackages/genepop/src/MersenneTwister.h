#ifndef MERSENNETWISTER_H
#define MERSENNETWISTER_H

#include <random>
class MTRand {
  std::mt19937 MTRAND;
  std::uniform_real_distribution<double> runif;
public:
  MTRand() {
    MTRAND=std::mt19937();
    runif=std::uniform_real_distribution<double>(0.0,1.0);
  }  
  void seed(unsigned long alea_seed) {MTRAND.seed(alea_seed);}
  double randExc(const double& r) { return runif(MTRAND)*r;}
  unsigned long int randInt(const unsigned long int& n ) {
    std::uniform_int_distribution<unsigned long int> sample(0, n);
    return(sample(MTRAND));
  }
  double operator()() { 
    return runif(MTRAND); 
    //return(double(MTRAND()) * (1.0/4294967296.0)); // similar to, but does not replicate the results with RJWagner's implementation.
  }  
};

#endif  // MERSENNETWISTER_H
