loadModule("Rcpp_s2net_export", TRUE)

setClass("Rcpp_s2net")
setMethod("predict", signature("Rcpp_s2net"), predict_Rcpp_s2net)

