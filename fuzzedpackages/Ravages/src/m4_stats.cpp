// ************* TAKEN FROM GASTON ***************
#include <Rcpp.h>
#include "gaston/matrix4.h"

using namespace Rcpp;

uint8_t N0[256] = {
4, 3, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
3, 2, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 
2, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};

uint8_t N1[256] = {
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
2, 3, 2, 2, 3, 4, 3, 3, 2, 3, 2, 2, 2, 3, 2, 2, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
1, 2, 1, 1, 2, 3, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 1, 0, 0, 1, 2, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0};


List geno_stats_snps(matrix4 & A, LogicalVector chr_xy, LogicalVector sexf) {

  // création vecteur de masques sex
  uint8_t * sex_ = new uint8_t[A.true_ncol];
  std::fill(sex_, sex_+A.true_ncol, 0);  
  int nbm = 0;
  for(size_t j = 0; j < A.true_ncol; j++) {
    if(!sexf(4*j)) { sex_[j] |= 3; nbm++; }
    if(4*j+1 < A.ncol && !sexf(4*j+1)) { sex_[j] |= 12;  nbm++; }
    if(4*j+2 < A.ncol && !sexf(4*j+2)) { sex_[j] |= 48;  nbm++; }
    if(4*j+3 < A.ncol && !sexf(4*j+3)) { sex_[j] |= 192; nbm++; }
  }

  // Stats SNP
  IntegerMatrix SN(4,A.nrow);
  IntegerMatrix SNf(4,A.nrow); // à remplir seulement pour chr X

  for(size_t i = 0; i < A.nrow; i++) {
    bool xy = chr_xy(i);
    for(size_t j = 0; j < A.true_ncol; j++) {
      uint8_t d = A.data[i][j];
      SN(0,i) += N0[d];
      SN(3,i) += N0[255-d];
      SN(1,i) += N1[d];
      SN(2,i) += N1[255-d];

      if(xy) { // Stats SNPs en ne prenant en compte que les femmes
        d |= sex_[j];
        SNf(0,i) += N0[d];
        SNf(3,i) += N0[255-d];
        SNf(1,i) += N1[d];
        SNf(2,i) += N1[255-d];
      }
    }
    // Des NAs en trop (la bordure)
    SN(3,i) -= (4*A.true_ncol - A.ncol);
    // Pour les femmes la bordure + le nbre d'hommes
    SNf(3,i) -= (4*A.true_ncol - A.ncol) + nbm;
  }
  List L;

  L["snps"] = DataFrame::create(Named("N0") = SN(0,_), Named("N1")  = SN(1,_), 
                                Named("N2") = SN(2,_), Named("NAs") = SN(3,_),
                                Named("N0.f") = SNf(0,_), Named("N1.f") = SNf(1,_),
                                Named("N2.f") = SNf(2,_), Named("NAs.f") = SNf(3,_));

  delete [] sex_;
  return L;
}

List geno_stats_snps(XPtr<matrix4> p_A,  LogicalVector chr_xy, LogicalVector sexf) {
  return geno_stats_snps(*p_A, chr_xy, sexf);
}

// geno_stats_snp
RcppExport SEXP gg_geno_stats_snps(SEXP p_ASEXP, SEXP chr_xySEXP, SEXP sexfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type chr_xy(chr_xySEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type sexf(sexfSEXP);
    __result = Rcpp::wrap(geno_stats_snps(p_A, chr_xy, sexf));
    return __result;
END_RCPP
}

