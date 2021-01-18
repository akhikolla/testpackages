#include <Rcpp.h>
#include <chunkR.h>
using namespace _chunkR;
using namespace Rcpp;


// [[Rcpp::export]]
RcppExport SEXP chunker__new_data_frame(SEXP path_, SEXP sep_, 
                                        SEXP quoted_, SEXP has_colnames_, 
                                        SEXP has_rownames_, SEXP chunksize_, 
                                        SEXP column_types_) {
  
  std::string path = Rcpp::as < std::string > (path_);
  char sep = Rcpp::as<char>(sep_);
  bool quoted = Rcpp::as<bool>(quoted_);
  bool has_colnames = Rcpp::as<bool>(has_colnames_);
  bool has_rownames = Rcpp::as<bool>(has_rownames_);
  unsigned int chunksize = Rcpp::as<unsigned int>(chunksize_);
  StringVector column_types = Rcpp::as<StringVector>(column_types_);
  Rcpp::XPtr < chunker > ptr(new chunker(path, sep, quoted,
                                         has_colnames, has_rownames, 
                                         chunksize,  column_types), 
                                         true);
  return ptr;
}

// [[Rcpp::export]]
RcppExport SEXP chunker__new_matrix(SEXP path_, SEXP sep_,  
                                    SEXP quoted_, SEXP has_colnames_, 
                                    SEXP has_rownames_, SEXP chunksize_) {
  
  std::string path = Rcpp::as < std::string > (path_);
  char sep = Rcpp::as<char>(sep_);
  bool quoted = Rcpp::as<bool>(quoted_);
  bool has_colnames = Rcpp::as<bool>(has_colnames_);
  bool has_rownames = Rcpp::as<bool>(has_rownames_);
  unsigned int chunksize = Rcpp::as<unsigned int>(chunksize_);
  Rcpp::XPtr < chunker > ptr(new chunker(path, sep, quoted, 
                                         has_colnames,  has_rownames, 
                                         chunksize), true);
  return ptr;
}

// [[Rcpp::export]]
RcppExport SEXP chunker__next_chunk(SEXP ptr) {
  Rcpp::XPtr < chunker > data(ptr);
  return wrap(data->next_chunk());
}

// [[Rcpp::export]]
RcppExport SEXP chunker__set_colnames(SEXP ptr) {
  Rcpp::XPtr < chunker > data(ptr);
  data->set_colnames();
  return wrap(true);
}

// [[Rcpp::export]]
RcppExport SEXP chunker__get_matrix(SEXP ptr) {
	Rcpp::XPtr < chunker > data(ptr);
	return data->get_matrix();
}

// [[Rcpp::export]]
RcppExport SEXP chunker__get_dataframe(SEXP ptr) {
  Rcpp::XPtr < chunker > data(ptr);
  return data->get_dataframe();
}

// [[Rcpp::export]]
RcppExport SEXP chunker__get_colnames(SEXP ptr) {
  Rcpp::XPtr < chunker > data(ptr);
  return data->get_colnames();
}


// [[Rcpp::export]]
RcppExport SEXP chunker__get_total(SEXP ptr) {
  Rcpp::XPtr < chunker > data(ptr);
  return wrap(data->get_total());
}

// [[Rcpp::export]]
RcppExport SEXP chunker__get_completed(SEXP ptr) {
	Rcpp::XPtr < chunker > data(ptr);
	return wrap(data->get_completed());
}

// [[Rcpp::export]]
RcppExport SEXP chunker__get_type(SEXP ptr) {
  Rcpp::XPtr < chunker > data(ptr);
  return wrap(data->get_type());
}
