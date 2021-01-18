#include "bbd.h"

// [[Rcpp::export]]
std::vector<double> SIR_Cpp(const double t, const double alpha, const double beta,
    const long int S0, const long int I0, const int Ap1, const int Bp1, const int direction, 
    const int nblocks, const double tol, const int computeMode, const int nThreads) {
      
      const int matsize = Ap1*Bp1;
      std::vector<double> lambda1(matsize), lambda2(matsize);
      
      for (int a=0; a<Ap1; ++a) {
        for (int b=0; b<Bp1; ++b) {
          double Spop = std::max(0.0, S0-a+0.0);
          double Ipop = std::max(0.0, a+I0-b+0.0);
          lambda1[a + b*Ap1] = beta*Spop*Ipop; // Infection rate is beta*S*I
          lambda2[a + b*Ap1] = alpha*Ipop;
        }
      }
      
      return(bb_lt_invert_Cpp(t, lambda1, lambda2, Ap1, Bp1, direction,
            nblocks, tol, computeMode, nThreads));
}
