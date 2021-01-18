paste_na_omit <- function(..., sep='') {
  .Call('_tokenbrowser_no_na_paste', PACKAGE = 'tokenbrowser', list(...), sep)
}

#' Rescale a numeric variable
#'
#' @param x a numeric vector
#' @param new_min The minimum value of the output
#' @param new_max The maximum value of the output
#' @param x_min The lowest possible value in x. By default this is the actual lowest value in x.
#' @param x_max The highest possible value in x. By default this is the actual highest value in x.
#'
#' @return a numeric vector
#' @export
#' @examples
#' rescale_var(1:10)
#' rescale_var(1:10, new_min = -1, new_max = 1)
rescale_var <- function(x, new_min=0, new_max=1, x_min=min(x), x_max=max(x)){
  if (x_min == x_max) return(x)
  x = (x - x_min) / (x_max - x_min) # normalize
  x = x * (new_max-new_min)
  return(x + new_min)
}


