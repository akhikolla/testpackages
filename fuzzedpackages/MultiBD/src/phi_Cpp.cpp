#include "bbd.h"

void phi_Cpp (const mytype::ComplexNumber s, const int a0, const int b0, const std::vector<double>& lambda2, 
    const std::vector<double>& mu2, const int A, const int Bp1, const int maxdepth, 
    mytype::ComplexVector& phi, const std::deque<std::vector<double>>& prod_mu2, 
    const std::deque<std::vector<double>>& prod_lambda2, const std::deque<std::vector<double>>& xvec, 
    const std::deque<std::vector<double>>& yvec_minus_s, mytype::ComplexVector& yvec, 
    mytype::ComplexVector& lentz_plus_invBk1dBk, mytype::ComplexVector& inv_Bk1dBk, 
    mytype::ComplexVector& BidBj) {  

/////////////////////////////////////////////////////////////////
// Compute continued fractions phi_{a,b}^m using Lentz algorithm
/////////////////////////////////////////////////////////////////
      
  for (int a=0; a<(A-a0+1); ++a) {
    for (int i=0; i<(Bp1 + maxdepth); ++i) yvec[i] = s + yvec_minus_s[a][i];
    
    inv_Bk1dBk_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk);
    lentz_plus_invBk1dBk_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk,lentz_plus_invBk1dBk);
    BidBj_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk,BidBj);
    
      
    std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1*(Bp1+1)/2),
//    unroll::for_each_4(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1*(Bp1+1)/2),
      [&](int k) {
        int j = (int)((-1+sqrt(8*k+1))/2);  
        int i = k - j*(j+1)/2;
        mytype::ComplexNumber tmp = BidBj[k]/lentz_plus_invBk1dBk[j];
        // When tmp == 0, set phi = 0 to avoid "0 x infinity" 
        // in the case prod_lambda2 or prod_mu2 is overflow
        if (tmp == zero) {
          phi[get_phi(a,i,j,Bp1)] = 0;            
          phi[get_phi(a,j,i,Bp1)] = 0;             
        } else {
          phi[get_phi(a,i,j,Bp1)] = prod_mu2[a][k] * tmp;            
          phi[get_phi(a,j,i,Bp1)] = prod_lambda2[a][k] * tmp;             
        }
      }
    );
    
//    unroll::vectorized_for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1*(Bp1+1)/2),
//      [&](int k) {
//        int j1 = (int)((-1+sqrt(8*k+1))/2);  
//        int i1 = k - j1*(j1+1)/2;
//        int j2 = (int)((-1+sqrt(8*(k+1)+1))/2);  
//        int i2 = k+1 - j2*(j2+1)/2;
//        
//        Complex4d BidBj_vec(BidBj[k],BidBj[k+1]);
//        Complex4d lentz_plus_invBk1dBk_vec(lentz_plus_invBk1dBk[j1],lentz_plus_invBk1dBk[j2]);
//        Complex4d tmp_vec = BidBj_vec/lentz_plus_invBk1dBk_vec;
//        phi[get_phi(a,i1,j1,Bp1)] = prod_mu2[a][k] * tmp_vec.get_low();            
//        phi[get_phi(a,j1,i1,Bp1)] = prod_lambda2[a][k] * tmp_vec.get_low();
//        phi[get_phi(a,i2,j2,Bp1)] = prod_mu2[a][k+1] * tmp_vec.get_high();            
//        phi[get_phi(a,j2,i2,Bp1)] = prod_lambda2[a][k+1] * tmp_vec.get_high();          
//      }
//    );
    
  }
}

