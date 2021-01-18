#' @template test_stat
#' @templateVar name The Miao, Gel and Gastwirth test statistic
#' @templateVar cite Miao2006
#' @templateVar formula \deqn{ \sup_{t>0}\frac{1}{n^2} \left| \sum_{i,j=1}^n I\{|X_i - X_j| < t\}- I\{|X_i + X_j| < t\}\right| }
#' @export
K2 <- function(X) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  K2_Cpp(X)
}

