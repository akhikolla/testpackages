#ifndef __HIPREAD_TYPES__
#define __HIPREAD_TYPES__

#include "datasource.h"
#include <Rcpp.h>


typedef Rcpp::XPtr<DataSource> XPtrDataSource;

// I don't really understand why, but this block is required
// to compile on macOS. I learned about it from
// https://github.com/RcppCore/Rcpp/pull/847
#ifdef FALSE
  #undef FALSE
#endif

#endif
