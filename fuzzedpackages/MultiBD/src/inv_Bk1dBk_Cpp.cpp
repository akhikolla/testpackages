#include "bbd.h"

void inv_Bk1dBk_Cpp(const int Bp1, const std::vector<double>& xvec, const mytype::ComplexVector& yvec, 
    mytype::ComplexVector& inv_Bk1dBk) {

////////////////////////
// Compute Y_{k}/Y_{k-1}
////////////////////////
  
  mytype::ComplexNumber Dj = zero, Dj1 = zero;
    for (int j=0; j<Bp1; ++j) {
        Dj = yvec[j] + xvec[j]*Dj1;
        if (Dj == zero) Dj = tiny;
        Dj1 = one/Dj;
        inv_Bk1dBk[j] = Dj;
    }
}
