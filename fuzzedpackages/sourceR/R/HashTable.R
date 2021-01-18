#' @importFrom methods new
HashTable = function(keys=character(0), vals=numeric(0)) {
  new("Rcpp_HashTable", keys, vals)
}
