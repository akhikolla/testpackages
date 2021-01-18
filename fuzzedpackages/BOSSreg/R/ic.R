#' Calculate an information criterion.
#'
#' Calculate a specified information criterion (IC) for an estimate or a group of estimates.
#' The choices of IC include AIC, BIC, AICc, BICc, GCV and Mallows' Cp.
#'
#' @param y_hat A vector of fitted values with \code{length(y_hat)=length(y)=n}, or
#'   a matrix, with \code{nrow(coef)=length(y)=n} and \code{ncol(y_hat)=m}, containing m different fits.
#' @param y A vector of response variable, with \code{length(y)=n}.
#' @param ic A specified IC to calculate. Default is AICc ('aicc'). Other choices include AIC ('aic'),
#'   BIC ('bic'), BICc ('bicc'), GCV ('gcv') and Mallows' Cp ('cp').
#' @param df A number if y_hat is a vector, or a vector with \code{length(df)=ncol(y_hat)=m} if y_hat is
#'   a matrix. df represents the degrees of freedom for each fit.
#' @param sigma Standard deviation of the error term. It only needs to be specified if the argument \code{ic='cp'}.
#'
#' @return The value(s) of the specified IC for each fit.
#'
#' @details This function enables the computation of various common IC for model fits, which can
#'   further be used to choose the optimal fit. This allows user comparing the effect of different IC.
#'   In order to calculate an IC, degrees of freedoms (df) needs to be specified. To be more specific,
#'   here are the formulas used to calculate each IC:
#'
#'   \deqn{AIC = \log(\frac{RSS}{n}) + 2\frac{df}{n}}{AIC = log(RSS/n) + 2*df/n}
#'   \deqn{BIC = \log(\frac{RSS}{n}) + \log(n)\frac{df}{n}}{BIC = log(RSS/n) + log(n)*df/n}
#'   \deqn{AICc = \log(\frac{RSS}{n}) + 2\frac{df+1}{n-df-2}}{AICc = log(RSS/n) + 2*(df+1)/(n-df-2)}
#'   \deqn{BICc = \log(\frac{RSS}{n}) + \log(n)\frac{df+1}{n-df-2}}{BICc = log(RSS/n) + log(n)*(df+1)/(n-df-2)}
#'   \deqn{GCV = \frac{RSS}{(n-df)^2}}{GCV = RSS/(n-df)^2}
#'   \deqn{Mallows' Cp = RSS + 2\times \sigma^2 \times df}{AIC = RSS + 2*\sigma^2*df}
#'
#' @author Sen Tian
#' @example R/example/eg.ic.R
#' @export
calc.ic <- function(y_hat, y, ic=c('aicc','bicc','aic','bic','gcv','cp'), df, sigma=NULL){
  # match the argument
  ic = match.arg(ic)

  # unify dimensions
  y = matrix(y, ncol=1)
  df = matrix(df, nrow=1)
  if(is.null(dim(y_hat))){
    y_hat = matrix(y_hat, ncol=1)
  }else if(dim(y_hat)[1]==1){
    y_hat = matrix(y_hat, ncol=1)
  }

  # sanity check
  if(ncol(y_hat) != ncol(df)){
    stop('the number of fits does not match the number of df')
  }
  if(ic=='cp' & is.null(sigma)){
    stop("need to specify sigma for Mallow's Cp")
  }

  n = nrow(y)
  nfit = ncol(y_hat)
  # for AICc and BICc df larger than n-2 will cause trouble, round it
  if(ic=='aicc' | ic=='bicc'){
    df[which(df>=n-2)]=n-3
  }

  rss = Matrix::colSums(sweep(y_hat, 1, y, '-')^2)

  if(ic=='aic'){return(log(rss/n)  + 2*df/n)}
  else if(ic=='bic'){return(log(rss/n)  + log(n)*df/n)}
  else if(ic=='aicc'){return(log(rss/n) + 2*(df+1)/(n-df-2))}
  else if(ic=='bicc'){return(log(rss/n) + log(n)*(df+1)/(n-df-2))}
  else if(ic=='gcv'){return(rss / (n-df)^2)}
  else if(ic=='cp'){return(rss + 2*sigma^2*df)}
}
