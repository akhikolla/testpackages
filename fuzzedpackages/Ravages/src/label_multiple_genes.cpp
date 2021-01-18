#include <Rcpp.h>
#include <iostream>
#include <boost/icl/interval_map.hpp>
// [[Rcpp::depends(BH)]]
using namespace Rcpp;

namespace icl = boost::icl;

typedef std::set<int> gene;
typedef std::pair<int, int> gene_pos;

// [[Rcpp::export]]
List label_multiple_genes(IntegerVector GENES_CHR, IntegerVector GENES_begin, IntegerVector GENES_end, IntegerVector VARIANTS_chr, IntegerVector VARIANTS_pos) {
  if(GENES_begin.size() != GENES_end.size() || GENES_begin.size() != GENES_CHR.size() || VARIANTS_chr.size() != VARIANTS_pos.size() )
    stop("Wrong dimensions");

   icl::interval_map< gene_pos, gene > MAP;
   for(int i = 0; i < GENES_begin.size(); i++) {
     gene G;
     G.insert(i+1);
     MAP += std::make_pair( icl::interval<gene_pos>::closed( gene_pos(GENES_CHR[i], GENES_begin[i]), gene_pos(GENES_CHR[i], GENES_end[i]) ), G );
   }

   List L(VARIANTS_chr.size());
   for(int i = 0; i < VARIANTS_chr.size(); i++) {
     gene_pos p( VARIANTS_chr[i], VARIANTS_pos[i] );
     auto r = MAP(p);
     std::vector<int> V(r.begin(), r.end());
     L[i] = wrap(V);
   }
   return L;
}


RcppExport SEXP label_multiple_genes(SEXP GENES_CHRSEXP, SEXP GENES_beginSEXP, SEXP GENES_endSEXP, SEXP VARIANTS_chrSEXP, SEXP VARIANTS_posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type GENES_CHR(GENES_CHRSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type GENES_begin(GENES_beginSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type GENES_end(GENES_endSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type VARIANTS_chr(VARIANTS_chrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type VARIANTS_pos(VARIANTS_posSEXP);
    rcpp_result_gen = Rcpp::wrap(label_multiple_genes(GENES_CHR, GENES_begin, GENES_end, VARIANTS_chr, VARIANTS_pos));
    return rcpp_result_gen;
END_RCPP
}
