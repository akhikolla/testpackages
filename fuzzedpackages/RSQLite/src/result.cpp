#include "pch.h"
#include "workarounds/XPtr.h"
#include "RSQLite_types.h"

//#include "DbResult.h"


// [[Rcpp::export]]
XPtr<DbResult> result_create(XPtr<DbConnectionPtr> con, std::string sql) {
  (*con)->check_connection();
  DbResult* res = SqliteResult::create_and_send_query(*con, sql);
  return XPtr<DbResult>(res, true);
}

// [[Rcpp::export]]
void result_release(XPtr<DbResult> res) {
  res.release();
}

// [[Rcpp::export]]
bool result_valid(XPtr<DbResult> res_) {
  DbResult* res = res_.get();
  return res != NULL && res->is_active();
}

// [[Rcpp::export]]
List result_fetch(DbResult* res, const int n) {
  return res->fetch(n);
}

// [[Rcpp::export]]
void result_bind(DbResult* res, List params) {
  res->bind(params);
}

// [[Rcpp::export]]
bool result_has_completed(DbResult* res) {
  return res->complete();
}

// [[Rcpp::export]]
int result_rows_fetched(DbResult* res) {
  return res->n_rows_fetched();
}

// [[Rcpp::export]]
int result_rows_affected(DbResult* res) {
  return res->n_rows_affected();
}

// [[Rcpp::export]]
List result_column_info(DbResult* res) {
  return res->get_column_info();
}

// [[Rcpp::export]]
CharacterVector result_get_placeholder_names(SqliteResult* res) {
  return res->get_placeholder_names();
}

namespace Rcpp {

template<>
DbResult* as(SEXP x) {
  DbResult* result = (DbResult*)(R_ExternalPtrAddr(x));
  if (!result)
    stop("Invalid result set");
  return result;
}

template<>
SqliteResult* as(SEXP x) {
  SqliteResult* result = (SqliteResult*)(R_ExternalPtrAddr(x));
  if (!result)
    stop("Invalid result set");
  return result;
}

}
