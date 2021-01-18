#include "bbd.h"

void bb_lt_Cpp(const mytype::ComplexNumber s, const std::vector<double>& lambda1, 
    const std::vector<double>& lambda2, const int Ap1, const int Bp1, const int direction,
    const std::vector<double>& yvec, mytype::ComplexVector& f) {  

/////////////////////////////////////////////////////////////////////////////////
// The following code computes Laplace transform of the transition probabilities      
///////////////////////////////////////////////////////////////////////////////// 

  if (direction == 0) { // Forward direction
    f[0] = 1/(s + yvec[0]);
    
    for (int i=0; i<(Ap1-1); ++i) {
      f[(i+1)*Bp1] = lambda1[i]*f[i*Bp1]/(s + yvec[i+1]);
    }
    
    for (int j=0; j<(Bp1-1); ++j) {
      f[j+1] = lambda2[j*Ap1]*f[j]/(s + yvec[(j+1)*Ap1]);
    }
    
    for (int i=0; i<(Ap1-1); ++i) {
      for (int j=0; j<(Bp1-1); ++j) {
        f[(i+1)*Bp1 + j+1] = (lambda1[i + (j+1)*Ap1]*f[i*Bp1 + j+1] + lambda2[i+1 + j*Ap1]*f[(i+1)*Bp1 + j])/(s + yvec[i+1 + (j+1)*Ap1]);
      }
    }
  } else { // Backward direction
    f[Ap1*Bp1-1] = 1/(s + yvec[Ap1*Bp1-1]);
    
    for (int i=(Ap1-2); i>=0; --i) {
      f[i*Bp1 + Bp1-1] = lambda1[i + (Bp1-1)*Ap1]*f[(i+1)*Bp1 + Bp1-1]/(s + yvec[i+ (Bp1-1)*Ap1]);
    }
    
    for (int j=(Bp1-2); j>=0; --j) {
      f[(Ap1-1)*Bp1 + j] = lambda2[Ap1-1 + j*Ap1]*f[(Ap1-1)*Bp1 + j+1]/(s + yvec[Ap1-1 + j*Ap1]);
    }
    
    for (int i=(Ap1-2); i>=0; --i) {
      for (int j=(Bp1-2); j>=0; --j) {
        f[i*Bp1 + j] = (lambda1[i + j*Ap1]*f[(i+1)*Bp1 + j] + lambda2[i + j*Ap1]*f[i*Bp1 + j+1])/(s + yvec[i + j*Ap1]);
      }
    }
  } 
  
}  
