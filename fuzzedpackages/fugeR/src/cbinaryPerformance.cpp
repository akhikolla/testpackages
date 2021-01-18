#include "cbinaryPerformance.h"

SEXP cbinaryPerformance( SEXP lstPredicted,
				         SEXP lstActual,
                         SEXP threshold ) {
				  
  using namespace Rcpp ;
  NumericVector xlstPredicted ( lstPredicted );
  NumericVector xlstActual    ( lstActual );
  NumericVector xthreshold    ( threshold );
  
  int _nbCase    = xlstActual.size();
  double _thresh = xthreshold[0];
  
  int truePos  = 0;
  int trueNeg  = 0;
  int falsePos = 0;
  int falseNeg = 0;
  
  double accu  = 0.0;
  double sensi = 0.0;
  double speci = 0.0;

  int tmpPred = 0;
  int tmpAct  = 0;
  for( int i = 0; i < _nbCase; i++ )
  {
    tmpPred = ( xlstPredicted[i] < _thresh ) ? (0) : (1);
    tmpAct  = ( xlstActual[i]    < _thresh ) ? (0) : (1);
    
    if( tmpPred && tmpAct)
    {
      truePos++;
    }
    else if ( !tmpPred && !tmpAct )
    {
      trueNeg++;
    }
    else if ( !tmpPred && tmpAct )
    {
      falseNeg++;
    }
    else if ( tmpPred && !tmpAct )
    {
      falsePos++;
    }
  }

  accu  = (double)( truePos + trueNeg ) /
      (double)(truePos + trueNeg + falsePos + falseNeg);
  sensi = (double)( truePos ) / (double)( truePos + falseNeg );
  speci = (double)( trueNeg ) / (double)( trueNeg + falsePos );

  return List::create( Named("accu")  = accu,
                       Named("sensi") = sensi,
                       Named("speci") = speci);
}
