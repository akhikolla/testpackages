#include "wrap.h"

namespace Rcpp {
template <>
::Encoding as(SEXP x) {
  if (Rf_isNull(x)) {
    return UNKNOWN_ENCODING;
  }
  ::Encoding res;
  EncodingFromName(Rcpp::as<std::string>(x).c_str(), &res);
  return res;
}

template <>
std::vector<::Encoding> as(SEXP x) {
  if (Rf_isNull(x)) {
    return std::vector<::Encoding>(1, UNKNOWN_ENCODING);
  }
  if (TYPEOF(x) != STRSXP) {
    Rcpp::stop("Input must be NULL or character vector. Provided %s", Rcpp::type2name(x));
  }
  std::vector<std::string> s = Rcpp::as<std::vector<std::string>>(x);
  std::vector<::Encoding> res(s.size());
  for(size_t i = 0; i < s.size(); ++i) {
    EncodingFromName(s[i].c_str(), &res[i]);
  }
  return res;
}

template <>
::Language as(SEXP x) {
  if (Rf_isNull(x)) {
    return UNKNOWN_LANGUAGE;
  }
  ::Language res;
  LanguageFromCode(Rcpp::as<std::string>(x).c_str(), &res);
  return res;
}

template <>
std::vector<::Language> as(SEXP x) {
  if (Rf_isNull(x)) {
    return std::vector<::Language>(1, UNKNOWN_LANGUAGE);
  }
  if (TYPEOF(x) != STRSXP) {
    Rcpp::stop("Input must be NULL or character vector. Provided %s", Rcpp::type2name(x));
  }
  std::vector<std::string> s = Rcpp::as<std::vector<std::string>>(x);
  std::vector<::Language> res(s.size());
  for(size_t i = 0; i < s.size(); ++i) {
    LanguageFromCode(s[i].c_str(), &res[i]);
  }
  return res;
}
} // Rcpp
