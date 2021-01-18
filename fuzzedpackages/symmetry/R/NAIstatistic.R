#' @template test_stat
#' @templateVar name Allison \eqn{T_2} statistic
#' @templateVar cite Klar2012
#' @templateVar formula \deqn{ \frac{1}{n {n\choose{k}}} \sum_{\mathcal{I}_{k}} \sum_{i_{k+1}=1}^n I\{|X_{(1),X_{i_1},\ldots,X_{i_{k}}}| < |X_{i_{k+1}}|\}- I\{|X_{(k),X_{i_1},\ldots,X_{i_{k}}}| < |X_{i_{k+1}}|\} }
#' @export
NAI <- function(X, k) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  NAI_Cpp(X, k)
}
