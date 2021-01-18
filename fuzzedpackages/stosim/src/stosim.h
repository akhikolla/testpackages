// stosim.h file
 /*
 * Author: David Silkworth
 *         (c) 2011-2018 OpenReliability.org
 */

#ifndef _stosim_H
#define _stosim_H

#include <Rcpp.h>

RcppExport SEXP withStorage(SEXP a, SEXP b, SEXP c);

RcppExport  SEXP SimulationHistory(SEXP arg1, SEXP arg2, SEXP arg3,
      SEXP arg4,  SEXP arg5, SEXP arg6,
        SEXP arg7,  SEXP arg8, SEXP arg9,
          SEXP arg10,  SEXP arg11, SEXP arg12);

 RcppExport  SEXP addWpush(SEXP arg1, SEXP arg2, SEXP arg3,
      SEXP arg4,  SEXP arg5, SEXP arg6,
        SEXP arg7,  SEXP arg8);

  RcppExport  SEXP addOverlay(SEXP arg1, SEXP arg2, SEXP arg3,
      SEXP arg4,  SEXP arg5, SEXP arg6,
        SEXP arg7,  SEXP arg8);
        
   RcppExport  SEXP DetailOpLinesCPP(SEXP arg1, SEXP arg2, SEXP arg3,
      SEXP arg4); 
      
   RcppExport  SEXP MultiTrainSingleBU4a(SEXP arg1, SEXP arg2, SEXP arg3,
      SEXP arg4,  SEXP arg5, SEXP arg6,
        SEXP arg7,  SEXP arg8);

   RcppExport  SEXP MultiTrainWithInventoryCPP(SEXP arg1, SEXP arg2, SEXP arg3,
      SEXP arg4); 		

#endif
