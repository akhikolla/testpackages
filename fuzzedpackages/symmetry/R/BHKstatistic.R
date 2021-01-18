#' @template test_stat
#' @templateVar name The Baringhaus and Henze test statistic
#' @templateVar cite BHK
#' @templateVar formula \deqn{ \sup_{t>0}\left|\frac{1}{{n\choose{2}}} \sum_{\mathcal{I}_{2}}\left( \frac12 I\{|X_{i_1}| <t\} + \frac12 I\{|X_{i_2}| <t\} - I\{|X_{(2),X_{i_1},X_{i_2}}| < t\}\right)\right| }
#' @export
BHK <- function(X) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  BHK_Cpp(X)
}

