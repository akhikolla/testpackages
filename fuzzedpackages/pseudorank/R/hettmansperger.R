#' @keywords internal
I <- function(d) {
  return(diag(d))
}

#' @keywords internal
J <- function(d) {
  return(rep(1,d)%*%t(rep(1,d)))
}


#' Hettmansperger-Norton Trend Test for k-Samples
#'
#' @description This function calculates the Hettmansperger-Norton trend test using pseudo-ranks under the null hypothesis H0F: F_1 = ... F_k.
#' @param data numeric vector containing the data
#' @param group ordered factor vector for the groups
#' @param alternative either decreasing or increasing
#' @param na.rm a logical value indicating if NA values should be removed
#' @param formula formula object
#' @param trend custom numeric vector indicating the trend for the custom alternative, only used if alternative = "custom"
#' @param ... further arguments are ignored
#' @return Returns a data.frame with the results
#' @example R/example_2.txt
#' @references Brunner, E., Bathke, A.C., and Konietschke, F. (2018a). Rank- and Pseudo-Rank Procedures for Independent Observations in Factorial Designs - Using R and SAS. Springer Series in Statistics, Springer, Heidelberg. ISBN: 978-3-030-02912-8.
#' @references Hettmansperger, T. P., & Norton, R. M. (1987). Tests for patterned alternatives in k-sample problems. Journal of the American Statistical Association, 82(397), 292-299
#' @keywords internal
hettmansperger_norton_test_internal <- function(data, group, na.rm, alternative = c("decreasing", "increasing", "custom"), formula = NULL, trend = NULL, pseudoranks = TRUE, ...) {

  stopifnot(is.numeric(data), is.factor(group), is.ordered(group), is.logical(na.rm), is.logical(pseudoranks))
  
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
  a <- length(n)
  
  # calculate either pseudo-ranks or ranks
  if(!pseudoranks) {
    df <- data.frame(pranks = rank(data), group = group)
  } else {
    df <- data.frame(pranks = pseudorank.numeric(data, group), group = group)
  }
  
  df <- df[order(df$group),]
  pHat <- 1/sum(n)*(summaryBy(pranks~group,data=df, FUN = mean)[, 2]-1/2)
  n <- summaryBy(pranks~group,data=df, FUN = length)[, 2]
  alternative <- match.arg(alternative)
  w <- rep(1, a)
  switch(alternative,
         decreasing={
           w <- a:1
         },
         increasing={
           w <- 1:a
         },
         custom={
           stopifnot(is.numeric(trend), length(trend)==a)
           w <- trend
         }
  )

  W <- diag(n)%*%(I(a) - 1/sum(n)*J(a)%*%diag(n))
  v2 <- 1/sum(n)^2*1/(sum(n)-1)*sum( (df$pranks - (sum(n)+1)/2 )^2 )
  sigmaHat2 <- sum(n)*v2*t(w)%*%W%*%diag(1/n)%*%W%*%w

  test <- sqrt(sum(n))*t(w)%*%W%*%pHat*1/sqrt(sigmaHat2)

  pValue <- 1 - pnorm(test)

  output <- list()
  output$name <- "Hettmansperger-Norton Trend Test"
  output$test <- test
  output$distribution <- "Standard-Normal"
  output$df <- NULL
  output$pValue <- pValue
  output$ss <- n
  output$pHat <- pHat
  output$alternative <- alternative
  output$formula <- formula
  output$trend <- trend
  output$pseudoranks <- pseudoranks
  class(output) <- "pseudorank"

  return(output)
}
