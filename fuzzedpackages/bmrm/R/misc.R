
#' Columun means of a matrix based on a grouping variable
#' 
#' Similar to rowsum, but for mean values.
#' @param x a matrix
#' @param group a factor with one element per row of x
#' @param ... additional arguments are passed to rowsum()
#' @return a matrix containing the means, with one row per level of group.
#' @export
rowmean <- function(x,group,...) {
  n <- rowsum(x,group,...)
  n / as.vector(table(group)[rownames(n)])
}


