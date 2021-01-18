#' Kepner-Robinson Test
#'
#' @description This function calculates the Kepner-Robinsin test under the null hypothesis H0F: F_1 = ... F_k.
#' @param data numeric vector containing the data
#' @param group ordered factor vector for the groups
#' @param alternative either decreasing or increasing
#' @param na.rm a logical value indicating if NA values should be removed
#' @param formula formula object
#' @param trend custom numeric vector indicating the trend for the custom alternative, only used if alternative = "custom"
#' @param ... further arguments are ignored
#' @return Returns a data.frame with the results
#' @example R/example_3.txt
#' @references Kepner, J. L., & Robinson, D. H. (1988). Nonparametric methods for detecting treatment effects in repeated-measures designs. Journal of the American Statistical Association, 83(402), 456-461.
#' @keywords internal
kepner_robinson_test_internal <- function(data, time, subject, distribution, na.rm, formula = NULL, ...) {
  
  stopifnot(is.numeric(data), is.factor(time), is.factor(subject), is.logical(na.rm))

  if(sum(is.na(data)) > 0) {
    if(na.rm) {
      # remove NAs
      nas <- which(is.na(data))
      rem <- vector(mode = "numeric", length=0L)
      for(i in 1:length(nas)) {
        rem <- c(rem, which(subject[nas[i]] == subject))
      }
      rem <- unique(rem)
      data <- data[-rem]
      time <- droplevels(time[-rem])
      subject <- droplevels(subject[-rem])
      
    } else {
      stop("There are missing values in your data. They can be omitted by choosing 'na.rm = TRUE'.")
    }
    
  }
  
  n <- length(unique(subject))
  a <- length(unique(time))
  N <- n*a
  
  df <- data.frame(pranks = rank(data), time = time, subject = subject)

  R_dot <- summaryBy(pranks~time,data=df, FUN = mean)[, 2]
  pHat <- 1/N*(R_dot-1/2)
  R_dot_ind <- summaryBy(pranks~subject,data=df, FUN = mean)[, 2]
  den <- 0
  for(i in 1:N) {
    den <- den + (df$pranks[i]-R_dot_ind[df$subject[i]])^2
  }

  test <- sum( (R_dot - (N+1)/2 )^2 ) * n^2*(a-1)*1/den

  pValue <- 1
  df <- a - 1
  if(distribution == "Chisq") {
    pValue <- 1 - pchisq(test, a-1)
  } else if(distribution == "F") {
    test <- test*1/(a-1)
    pValue <- 1 - pf(test, a-1, n*(a-1))
    df <- c(a-1, n*(a-1))
  }

  
  output <- list()
  output$name <- "Kepner-Robinson Test"
  output$test <- test
  output$distribution <- distribution
  output$df <- df
  output$pValue <- pValue
  output$ss <- n
  output$pHat <- pHat
  output$alternative <- "two.sided"
  output$formula <- formula
  output$trend <- NULL
  output$pseudoranks <- FALSE
  class(output) <- "pseudorank"
  
  return(output)
}