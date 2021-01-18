//[[Rcpp::depends(gRbase,RcppArmadillo,RcppEigen)]]

// tabMarg__ tabMult__ tabDiv0__

#include <RcppArmadillo.h>
#include "gRbase.h"

using namespace Rcpp;
using namespace gRbase;


// //[[Rcpp::export]]
// NumericVector check__(List cqpotList_){
//   NumericVector out = tabMult__(cqpotList_[0], cqpotList_[1]);
//   return out;
// }


//[[Rcpp::export]]
List propagateLS__(List cqpotList_, List rip){
  List cqpotList = clone(cqpotList_);

  List cliq  = rip["cliques"], seps = rip["separators"], childList = rip["childList"];
  IntegerVector pa = rip["parents"];
  
  NumericVector sp_pot, cq_pot, pa_pot, tmpd;
  CharacterVector cq, sp;
  IntegerVector ch, tmp;
  double normConst, zd;
  int i, j, ncliq = cliq.length(), nch, idx;

  for (i = ncliq - 1; i > 0; --i){
    //std::cout << " collect    i= " << i << std::endl;
    cq = cliq[i];            //Rprintf("cq\n");     Rf_PrintValue( cq );
    sp = seps[i];            //Rprintf("sp\n");     Rf_PrintValue( sp );
    cq_pot = cqpotList[ i ]; //Rprintf("cq_pot\n"); Rf_PrintValue( cq_pot );
    if (sp.size() >= 1){
      idx    = pa[ i ] - 1;
      pa_pot = cqpotList[ idx ];           //Rprintf("pa_pot\n"); Rf_PrintValue(pa_pot);
      //sp_pot = tabMarg__(cq_pot, sp);     //Rprintf("sp_pot\n"); Rf_PrintValue(sp_pot);
      sp_pot = tab_marg_(cq_pot, sp);     //Rprintf("sp_pot\n"); Rf_PrintValue(sp_pot);
      //cqpotList[ i ]   = tabDiv0__(cq_pot, sp_pot);
      cqpotList[ i ]   = tab_div0_(cq_pot, sp_pot);
      //cqpotList[ idx ] = tabMult__(pa_pot, sp_pot);
      cqpotList[ idx ] = tab_mult_(pa_pot, sp_pot);
    } else {
      zd = sum( cq_pot );
      tmpd = cqpotList[0];
      tmpd = zd * tmpd;
      cqpotList[0] = tmpd;
      tmpd = cqpotList[i];
      tmpd = tmpd / zd; 
      cqpotList[i] = tmpd;
    }
  }

  //Rcout << "JJJJJJJJJJJJJJJJJJJJ" << std::endl;
  tmpd = cqpotList[0];
  //Rf_PrintValue(tmpd);
  normConst = sum( tmpd );
  tmpd = tmpd / normConst;
  cqpotList[0] = tmpd;
  //Rf_PrintValue(tmpd);	
  
  for (i=0; i < ncliq; ++i){
    // std::cout << " distribute i= " << i << std::endl;
    ch  = childList[i];
    nch = ch.size();
    if (nch > 0){
      for (j=0; j < nch; ++j){
	idx = ch[ j ] - 1;
	sp  = seps[ idx ];
	if(sp.size() > 0){
	  // Rf_PrintValue(sp);
	  //sp_pot = tabMarg__(cqpotList[i], sp);
	  sp_pot = tab_marg_(cqpotList[i], sp);
	  //cqpotList[ idx ] = tabMult__( cqpotList[ idx ], sp_pot);
	  cqpotList[ idx ] = tab_mult_( cqpotList[ idx ], sp_pot);
	}
      }
    }
  }
  
  cqpotList.attr("pEvidence") = normConst;
  return cqpotList;
}


