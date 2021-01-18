#include "cfuzzifyVar.h"

SEXP cfuzzifyVar( SEXP mfId,
				  SEXP mfPoints,
                  SEXP values,
				  SEXP minValues ) {
				  
  using namespace Rcpp ;
  NumericVector mf( mfId );
  NumericVector mfList( mfPoints );
  NumericVector data( values );
  NumericVector minVal( minValues );
  int id = mf[0] - 1;
  int nbValues = data.size();
  NumericVector res(nbValues, 2.0);

  //Sort mfList ONLY 2 VALUES
  if ( mfList[0] > mfList[1] ){
    double tmp = mfList[0];
    mfList[0] = mfList[1];
    mfList[1] = tmp;
  }

  for( int i = 0; i < nbValues; i++ )
  {
    /*
    if( data[i] != data[i] )
    {
      continue;
    }
    */
    if(ISNAN(data[i]))
      continue;

    if( id == 0 )
    {
      if( data[i] <= mfList[id] )
      {
        res[i] = 1.0;
      }
      else if ( data[i] >= mfList[id+1] )
      {
        res[i] = 0.0;
      }
      else
      {
        res[i] = 1.0 - ( ( data[i] - mfList[id] )  /
                 ( mfList[id+1] - mfList[id] ) );
      }
    }
    else if( id == 1)
    {
      if( data[i] <= mfList[id-1] )
      {
        res[i] = 0.0;
      }
      else if ( data[i] >= mfList[id] )
      {
        res[i] = 1.0;
      }
      else
      {
        res[i] = ( ( data[i] - mfList[id-1] )  /
                 ( mfList[id] - mfList[id-1] ) );
      } 
    }

    //We keep only the min value because we apply the AND (min) operator
    if( minVal[i]  < res[i] )
    {
      res[i] = minVal[i];
    }
  }
  return res;
}
