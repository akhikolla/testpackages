#include <Rcpp.h>
#include <iostream>
#include <map>
using namespace Rcpp;


// CHR, begin, end = positions de début / fin des genes
// chr, pos = position des SNPs à étiquetter
// dans le résultat, les labels commencent à 1 [0 = pas de label]
// [[Rcpp::export]]
IntegerVector label_genes(IntegerVector CHR, IntegerVector begin, IntegerVector end, IntegerVector chr, IntegerVector pos) {
  if(begin.size() != end.size() || begin.size() != CHR.size() || chr.size() != pos.size())
    stop("brx");

  // on fait une map de ( (chr, begin), (end, index) )
  std::map< std::pair<int, int>, std::pair<int, int> > intervals;
  for(int i = 0; i < begin.size(); i++)
    intervals.insert( std::pair< std::pair<int, int>, std::pair<int, int>>( std::pair<int, int>(CHR[i], begin[i]), std::pair<int, int>(end[i], i+1)) );

  IntegerVector R( pos.size() );
  for(int j = 0; j < pos.size(); j++) {
    // on cherche le premier gène dont la position est > que celle du SNP
    auto it = intervals.upper_bound( std::pair<int, int>(chr[j], pos[j]) );
    it--;  // on revient au gene d'avant
    auto x1 = it->first;  // (chr, begin)
    auto x2 = it->second; // (end, index)
    // si le SNP fait partie de ce gène : on est bons.
    if(x1.first == chr[j] && x2.first >= pos[j]) R[j] = x2.second;
  }
  return R;
}

RcppExport SEXP label_genes(SEXP CHRSEXP, SEXP beginSEXP, SEXP endSEXP, SEXP chrSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type CHR(CHRSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type begin(beginSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type end(endSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(label_genes(CHR, begin, end, chr, pos));
    return rcpp_result_gen;
END_RCPP
}

