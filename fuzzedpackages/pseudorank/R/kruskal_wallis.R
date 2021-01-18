#' Hettmansperger-Norton Trend Test for k-Samples
#'
#' @description This function calculates the Kruskal-Wallis test using pseudo-ranks under the null hypothesis H0F: F_1 = ... F_k.
#' @param data numeric vector containing the data
#' @param group factor specifying the groups
#' @param na.rm a logical value indicating if NA values should be removed
#' @param formula formula object
#' @param pseudoranks logical value indicating if pseudo-ranks or ranks should be used
#' @param ... further arguments are ignored
#' @return Returns a data.frame with the results
#' @example R/example_3.txt
#' @references Brunner, E., Bathke, A.C., and Konietschke, F. (2018a). Rank- and Pseudo-Rank Procedures for Independent Observations in Factorial Designs - Using R and SAS. Springer Series in Statistics, Springer, Heidelberg. ISBN: 978-3-030-02912-8.
#' @references Hettmansperger, T. P., & Norton, R. M. (1987). Tests for patterned alternatives in k-sample problems. Journal of the American Statistical Association, 82(397), 292-299
#' @keywords internal
kruskal_wallis_internal <- function(data, group, na.rm, formula = NULL, pseudoranks = TRUE, ...) {
  
  stopifnot(is.numeric(data), is.factor(group), is.logical(na.rm), is.logical(pseudoranks))

  if(sum(is.na(data)) > 0) {
    if(na.rm) {
      # remove NAs
      nas <- which(is.na(data))
      data <- data[-nas]
      group <- group[-nas]
    } else {
      stop("There are missing values in your data. They can be omitted by choosing 'na.rm = TRUE'.")
    }
  }
  
  n <- as.numeric(as.matrix(table(group)))
  N <- sum(n)
  a <- length(n)
  if(pseudoranks){
    df <- data.frame(pranks = pseudorank(data, group), group = group)
    
  } else {
    df <- data.frame(pranks = pseudorank(data, as.factor(rep(1, length(data)))), group = group)
  }
  df <- df[order(df$group),]
  R_mean <- summaryBy(pranks~group,data=df, FUN = mean)[, 2]
  
  # test statistic
  numerator = sum(n*( R_mean - (N+1)/2)^2)
  denominator = sum(( df$pranks - (N+1)/2 )^2)
  
  Q_N = numerator/denominator*(N - 1)
  
  pValue = 1 - pchisq(Q_N, a - 1)
  
  output <- list()
  output$name <- "Kruskal-Wallis Test"
  output$test <- Q_N
  output$distribution <- "Chisq"
  output$df <- a-1
  output$pValue <- pValue
  output$ss <- n
  output$pHat <- 1/N*(R_mean-1/2)
  output$formula <- formula
  output$pseudoranks <- pseudoranks
  output$alternative <- NULL
  output$trend <- NULL
  class(output) <- "pseudorank"
  
  return(output)
  
}