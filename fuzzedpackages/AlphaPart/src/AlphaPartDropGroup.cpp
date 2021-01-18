#include "AlphaPartDropGroup.h"

SEXP AlphaPartDropGroup(SEXP c1_, SEXP c2_, SEXP nI_, SEXP nP_, SEXP nT_, SEXP nG_, SEXP y_, SEXP P_, SEXP Px_, SEXP g_)
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
  int nG = Rcpp::as<int>(nG_);
  Rcpp::NumericMatrix ped(y_);
  Rcpp::IntegerVector P(P_);  
  Rcpp::IntegerVector Px(Px_);
  Rcpp::IntegerVector g(g_);  
  
  // --- Outputs ---
      
  Rcpp::NumericMatrix pa(nI+1, nT);    // parent average
  Rcpp::NumericMatrix  w(nI+1, nT);    // Mendelian sampling
  Rcpp::NumericMatrix xa(nI+1, nP*nT); // Parts
  Rcpp::NumericMatrix xg(nG+1, nP*nT); // Parts for groups

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

      // ... for the PA part
      for(p = 0; p < nP; p++) {
        j = Px[t] + p;
        xa(i, j) += c1 * xa(ped(i, 1), j) +
                    c2 * xa(ped(i, 2), j);
        // Rprintf("Animal: %i, Trait: %i, Path: %i (%i), Group: %i\\n", i, t, p, j, g(i));
        xg(g(i), j) += xa(i, j);
      }

    }
  }
  
  // --- Return ---

  return xg;

}
