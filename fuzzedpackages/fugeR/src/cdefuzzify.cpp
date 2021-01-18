#include "cdefuzzify.h"

SEXP cdefuzzify( SEXP nbCase,
				 SEXP nbRule,
                 SEXP inRule,
				 SEXP lstMf,
				 SEXP lstMfId,
                 SEXP defaultMfId,
				 SEXP lstActivation ) {
  
  using namespace Rcpp ;
  NumericVector xnbCase        ( nbCase );
  NumericVector xnbRule        ( nbRule );
  LogicalVector xinRule        ( inRule );
  NumericVector xlstMf         ( lstMf );
  NumericVector xlstMfId       ( lstMfId );
  NumericVector xdefaultMfId   ( defaultMfId );
  List          xlstActivation ( lstActivation );
  
  const int _nbCase = xnbCase[0];
  const int _nbRule = xnbRule[0];
  
  //vector<vector<point> > a(10, vector<point>(10));
  std::vector<NumericVector> rulesActivation( _nbRule - 1 );
  //NumericVector rulesActivation[_nbRule - 1];
  for( int j = 0; j < _nbRule - 1; j++ )
    rulesActivation[j] = as<SEXP>(xlstActivation[j]);

  //Prediction to return
  NumericVector res   ( _nbCase, 0.0 );

  double act;
  double maxAct;
  double num;
  double denum;

  for( int i = 0; i < _nbCase; i++ )
  {    
    //Reset
    act = 0.0;
    maxAct = 0.0;
    num = 0.0;
    denum = 0.0;

    for( int j = 0; j < _nbRule - 1; j++ )
    {
      act = rulesActivation[j][i];
      if( act <= 1.0 && xinRule[j] )
      {
        num    = num   + ( act * xlstMf[xlstMfId[j]-1] );
        denum  = denum + ( act );
        if( act > maxAct )
        {
          maxAct = act;
        }
      }
    }

    //Apply default Rule 1-MaxActivation
    num   = num   + ( ( 1.0 - maxAct ) * xlstMf[xdefaultMfId[0]-1] );
    denum = denum + ( 1.0 - maxAct );

    res[i] = num / denum;
  }

  return res;
}
