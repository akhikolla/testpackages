#include <Rcpp.h>
#include "wrap.h"

Rcpp::String detect(const char* x, size_t n, const Encoding& enc, const Language& lang) {
  if (n == 0) {
    return NA_STRING;
  }
  bool is_reliable;
  int bytes_consumed;
  Encoding res = DetectEncoding(
    x,                 // input string
    n,                 // input length
    nullptr,           // hint: web-page URL or domain
    nullptr,           // hint: web page HTTP header charset
    nullptr,           // hint: web page meta tag charset
    enc,               // hint: encoding
    lang,              // hint: language
    QUERY_CORPUS,      // type of corpus
    false,             // ignore 7=bit mail encodings
    &bytes_consumed,   // how many bytes used
    &is_reliable       // encoding is at least 2^10 time more probable then the second-best encoding
  );
  // if (!is_reliable) {
  //   return NA_STRING;
  // }
  return Rcpp::wrap(MimeEncodingName(res));
}

//' @title
//' Detect Encoding
//'
//' @description
//' Detect charset encoding of the character or raw vector.
//'
//' @param x Raw or character vector.
//' @param enc_hint Character vector with encoding hint.
//' @param lang_hint Character vector with langauge code hint.
//'
//' @return Character vector with suggested encodings.
//'
//' @export
//'
//' @example man-roxygen/ex-detect.R
//'
// [[Rcpp::export(rng = false)]]
SEXP ced_enc_detect(SEXP x, SEXP enc_hint = R_NilValue, SEXP lang_hint = R_NilValue) {
  size_t n = LENGTH(x);
  if (n == 0) {
    return Rcpp::StringVector();
  }
  if (TYPEOF(x) == RAWSXP) {
    const char* s = reinterpret_cast<const char*>(RAW(x));
    Encoding enc = Rcpp::as<Encoding>(enc_hint);
    Language lang = Rcpp::as<Language>(lang_hint);
    return Rcpp::wrap(detect(s, n, enc, lang));
  }
  if (TYPEOF(x) != STRSXP) {
    Rcpp::stop("'x' must raw or character vector.");
  }
  Rcpp::StringVector xx = Rcpp::as<Rcpp::StringVector>(x);
  std::vector<Encoding> enc = Rcpp::as<std::vector<Encoding>>(enc_hint);
  size_t e_n = enc.size();
  if (e_n != 1 && e_n != n) {
    Rcpp::stop("'enc_hint' must be length one or equal length of the input character vector.");
  }
  std::vector<Language> lang = Rcpp::as<std::vector<Language>>(lang_hint);
  size_t l_n = lang.size();
  if (l_n != 1 && l_n != n) {
    Rcpp::stop("'lang_hint' be length one or equal length of the input character vector.");
  }
  Rcpp::StringVector res = Rcpp::no_init(n);
  for (size_t i = 0; i < n; ++i) {
    if (Rcpp::StringVector::is_na(xx[i])) {
      res[i] = NA_STRING;
    } else {
      std::string s = Rcpp::as<std::string>(xx[i]);
      res[i] = detect(s.c_str(), s.size(), enc[i % e_n], lang[i % l_n]);
    }
  }
  if (xx.hasAttribute("names")) {
    res.names() = xx.names();
  }
  return wrap(res);
}
