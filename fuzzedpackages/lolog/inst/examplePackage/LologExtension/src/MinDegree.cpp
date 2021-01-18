// [[Rcpp::depends(lolog)]]
#include "MinDegree.h"
#include <lolog.h>


//[[Rcpp::export()]]
void registerMinDegree(){
  Rcpp::XPtr< lolog::AbstractStat<lolog::Undirected> > ps1(new lologext::UndirectedMinDegree());
  REGISTER_UNDIRECTED_STATISTIC(ps1);
}
