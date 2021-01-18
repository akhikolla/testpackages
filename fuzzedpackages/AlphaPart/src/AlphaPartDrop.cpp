#include "AlphaPartDrop.h"

SEXP AlphaPartDrop(SEXP c1_, SEXP c2_, SEXP nI_, SEXP nP_, SEXP nT_, SEXP y_, SEXP P_, SEXP Px_)
{

  using namespace Rcpp ;
  //' @export

  // --- Temp ---
      
  int i, j, k, t, p;
  
  // --- Inputs ---
      
  double c1 = Rcpp::as<double>(c1_);  
  double c2 = Rcpp::as<double>(c2_); 
  int nI = Rcpp::as<int>(nI_); 
  int nP = Rcpp::as<int>(nP_);
  int nT = Rcpp::as<int>(nT_);
  Rcpp::NumericMatrix ped(y_);
  Rcpp::IntegerVector P(P_);  
  Rcpp::IntegerVector Px(Px_);
  
  // --- Outputs ---
      
  Rcpp::NumericMatrix pa(nI+1, nT);    // parent average
  Rcpp::NumericMatrix  w(nI+1, nT);    // Mendelian sampling
  Rcpp::NumericMatrix xa(nI+1, nP*nT); // Parts
  
  // --- Compute ---
      
  for(i = 1; i < nI+1; i++) {
    for(t = 0; t < nT; t++) {
      // Parent average (PA)
      pa(i, t) = c1 * ped(ped(i, 1), 3+t) +
                 c2 * ped(ped(i, 2), 3+t);
    
      // Mendelian sampling (MS)
      w(i, t) = ped(i, 3+t) - pa(i, t);
    
      // Parts

      // ... for the MS part
      j = Px[t] + P[i];
      xa(i, j) = w(i, t);

      // ... for the PA parts
      for(p = 0; p < nP; p++) {
        j = Px[t] + p;
        xa(i, j) += c1 * xa(ped(i, 1), j) +
                    c2 * xa(ped(i, 2), j);
      }
    }
  }
  
  // --- Return ---

  return Rcpp::List::create(Rcpp::Named("pa", pa),
                            Rcpp::Named("w",  w),
                            Rcpp::Named("xa", xa));

}
