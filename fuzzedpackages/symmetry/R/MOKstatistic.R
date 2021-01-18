#' @template test_stat
#' @templateVar name The Mira test statistic
#' @templateVar cite Mira1999
#' @templateVar formula \deqn{ \sup_{t>0}\left|\frac{1}{{n\choose{2k}}} \sum_{\mathcal{I}_{2k}} I\{|X_{(k),X_{i_1},\ldots,X_{i_{2k}}}| < t\}- I\{|X_{(k+1),X_{i_1},\ldots,X_{i_{2k}}}| < t \}\right| }
#' @export
MOK <- function(X, k) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  MOK_Cpp(X, k)
}

