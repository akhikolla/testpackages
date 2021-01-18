#include "bbd.h"

void BidBj_Cpp(const int Bp1, const std::vector<double>& xvec, const mytype::ComplexVector& yvec, 
    const mytype::ComplexVector& inv_Bk1dBk, mytype::ComplexVector& BidBj) {
      
///////////////////////
// Compute Y_{i}/Y_{j}
///////////////////////
    
  std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1-1),
//  unroll::for_each_4(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1-1),
    [&](int i) {
      BidBj[Trimat(i,i)] = one;
      BidBj[Trimat(i,i+1)] = one/inv_Bk1dBk[i];
      for (int j=(i+2); j<Bp1; ++j) {
        mytype::ComplexNumber tmp = yvec[j-1]/BidBj[Trimat(i,j-1)] + xvec[j-1]/BidBj[Trimat(i,j-2)];
        BidBj[Trimat(i,j)] = one/tmp;
        if (BidBj[Trimat(i,j)]==zero) {std::fill_n(&BidBj[Trimat(i,j)],Bp1-j,zero);break;}  
      }
    });
      
  BidBj[Trimat(Bp1-1,Bp1-1)] = one;  
}
