



#' An lolog plug-in for easy C++ prototyping and access
#' @param ... plug-in arguments
#' @details
#' The lolog Rcpp plugin allows for the rapid prototyping of compiled code.
#' new functions can be registered and exposed using \code{\link{cppFunction}}
#' and new statistics can be compiled and registered using \code{\link{sourceCpp}}.
#' @examples
#' \dontrun{
#' # This creates a function in C++ to create an empty network of size n
#' # and expose it to R.
#' src <- "
#' lolog::BinaryNet<lolog::Directed> makeEmptyNetwork(const int n){
#' Rcpp::IntegerMatrix tmp(0,2);
#' lolog::BinaryNet<lolog::Directed> net(tmp, n);
#' return net;
#' }
#' "
#' Rcpp::registerPlugin("lolog",inlineLologPlugin)
#' emptyNetwork <- cppFunction(src,plugin="lolog")
#' net <- emptyNetwork(10L)
#' net[1:10,1:10]
#' 
#' }
#' @seealso \code{\link{cppFunction}}, \code{\link{sourceCpp}}, \code{\link{cppFunction}}
inlineLologPlugin <- Rcpp::Rcpp.plugin.maker(
  include.after = "#include <lolog.h>",
  LinkingTo = unique(c("lolog", "BH", "Rcpp")),
  Depends = unique(c("lolog", "BH", "Rcpp")),
  Imports = unique(c("lolog", "BH", "Rcpp")),
  package        = "lolog"
)
