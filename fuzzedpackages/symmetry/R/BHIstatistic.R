#' @template test_stat
#' @templateVar name The Litvinova test statistic
#' @templateVar cite BHI
#' @templateVar formula \deqn{ \frac{1}{n {n\choose{2}}} \sum_{\mathcal{I}_{2}} \sum_{i_{3}=1}^n \left( \frac12 I\{|X_{i_1}| < |X_{i_{3}}|\} + \frac12 I\{|X_{i_2}| < |X_{i_{3}}|\} - I\{|X_{(2),X_{i_1},X_{i_2}}| < |X_{i_{3}}|\}\right) }
#' @export
BHI <- function(X) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  BHI_Cpp(X)
}
