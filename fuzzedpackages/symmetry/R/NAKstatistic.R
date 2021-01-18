#' @template test_stat
#' @templateVar name Allison \eqn{T_2} statistic
#' @templateVar cite Klar2012
#' @templateVar formula \deqn{ \sup_{t>0}\left|\frac{1}{{n\choose{k}}} \sum_{\mathcal{I}_{k}} I\{|X_{(1),X_{i_1},\ldots,X_{i_{k}}}| < t\}- I\{|X_{(k),X_{i_1},\ldots,X_{i_{k}}}| < t \}\right| }
#' @export
NAK <- function(X, k) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  NAK_Cpp(X, k)
}

