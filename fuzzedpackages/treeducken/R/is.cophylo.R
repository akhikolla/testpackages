#' Test for the cophygenetic set object
#' @description Tests if an object is of class `cophy`
#'
#' @param cophy an object to test to see if it is of class `cophy`
#' @details Checks that an object is of class `cophy`. For multicophy checks that
#' the class is `multiCophylo` and that each element is of class `cophy`.
#' @return A logical vector
#' @seealso as.cophy
#' @examples
#' h_lambda <- 1.0
#' h_mu <- 0.3
#' c_lambda <- 0.0
#' s_lambda <- 1.0
#' s_mu <- 0.3
#' s_her <- 0.0
#' host_symb_sets <- sim_cophylo_bdp(hbr = h_lambda,
#'                                   hdr = h_mu,
#'                                   sbr = s_lambda,
#'                                   cosp_rate = c_lambda,
#'                                   sdr = s_mu,
#'                                   host_exp_rate = s_her,
#'                                   time_to_sim = 2.0,
#'                                   numbsim = 1)
#' is.cophy(host_symb_sets[[1]])
#' is.multiCophylo(host_symb_sets)
#' @export
is.cophy <- function(cophy){
    inherits(cophy, "cophy")
}
#' @describeIn is.cophy Tests for `multiCophylo` composed of `cophy` objects
#' @param multiCophy an object to test for multiCophy
#' @export
is.multiCophylo <- function(multiCophy){
    t <- inherits(multiCophy, "multiCophy")
    if(t){
        tt <- sapply(unclass(multiCophy), inherits, what = "cophy")
    }
    all(c(t, tt))
}