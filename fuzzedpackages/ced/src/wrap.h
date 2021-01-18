#include <RcppCommon.h>
#include "ced.h"

namespace Rcpp {
template <>
::Encoding as(SEXP);
template <>
::Language as(SEXP);
template <>
std::vector<::Encoding> as(SEXP);
template <>
std::vector<::Language> as(SEXP);

} // Rcpp

#include <Rcpp.h>
