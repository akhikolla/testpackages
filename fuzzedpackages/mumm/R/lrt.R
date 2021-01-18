
#' Likelihood Ratio Test
#'
#' A function to perform a likelihood ratio test for testing two nested models against each other.
#'
#' @param fit1 a fitted model object of class \code{mumm}.
#'
#' @param fit2 a fitted model object of class \code{mumm}, \code{lm} or \code{merMod}.
#'
#' @details Performs the likelihood ratio test for testing two nested models against each other. The model in
#' \code{fit2} should be nested within the model in \code{fit1}.
#'
#' @return A matrix with the likelihood ratio test statistic and the corresponding p-value.
#'
#' @examples
#' set.seed(100)
#' sigma_e <- 1.5
#' sigma_a <- 0.8
#' sigma_b <- 0.5
#' sigma_d <- 0.7
#' nu <- c(8.2, 6.2, 2.3, 10.4, 7.5, 1.9)
#'
#' nA <- 15
#' nP <- 6
#' nR <- 5
#'
#' a <- rnorm(nA, mean = 0, sd = sigma_a)
#' b <- rnorm(nA, mean = 0, sd = sigma_b)
#' d <- rnorm(nA*nP, mean = 0, sd = sigma_d)
#' e <- rnorm(nA*nP*nR, mean = 0, sd = sigma_e)
#'
#' Assessor <- factor(rep(seq(1,nA),each = (nP*nR)))
#' Product <- factor(rep(rep(seq(1,nP),each = nR), nA))
#' AssessorProduct <- (Assessor:Product)
#'
#' y <- nu[Product] + a[Assessor] + b[Assessor]*(nu[Product]-mean(nu)) + d[AssessorProduct] + e
#'
#' sim_data <- data.frame(y, Assessor, Product)
#'
#' fit <- mumm(y ~ 1 + Product + (1|Assessor) + (1|Assessor:Product) +
#'              mp(Assessor,Product) ,data = sim_data)
#'
#' fit2 <- mumm(y ~ 1 + Product + (1|Assessor) + mp(Assessor,Product) ,data = sim_data)
#' lrt(fit,fit2)
#'
#' @export
lrt <- function(fit1, fit2) {

  cat("Data:", fit1$call$data, "\n")
  cat("Models: \n Object:")

  for(i in 1:length(deparse(fit1$call$formula,width.cutoff = 50))) {
    if(i>1){cat("    ")}
    cat(deparse(fit1$call$formula,width.cutoff = 50)[i],"\n")
  }


  cat("...1   :")

  if(class(fit2)=="mumm") {

    loglik2 = -fit2$objective
    df2 = fit2$df

    for(i in 1:length(deparse(fit2$call$formula,width.cutoff = 50))) {
      if(i>1){cat("    ")}
      cat(deparse(fit2$call$formula,width.cutoff = 50)[i],"\n")
    }
  } else {
    loglikelihood = stats::logLik(fit2, REML = F)
    loglik2 = loglikelihood[1]
    df2 = attr(loglikelihood,"df")

    if(class(fit2)=="lm"){

      for(i in 1:length(deparse(fit2$call$formula,width.cutoff = 50))) {
        if(i>1){cat("    ")}
        cat(deparse(fit2$call$formula,width.cutoff = 50)[i],"\n")
      }

    } else {

      for(i in 1:length(deparse(fit2@call$formula,width.cutoff = 50))) {
        if(i>1){cat("    ")}
        cat(deparse( fit2@call$formula,width.cutoff = 50)[i],"\n")
      }
    }

  }

  #table
  loglik1 = -fit1$objective
  df1 = fit1$df
  df = df1-df2
  chi = 2*(loglik1-loglik2)
  pvalue = 1-stats::pchisq(chi,df  = df)
  table_sd = data.frame(Df = c(df1,df2), logLik = c(loglik1,loglik2), Chisq = c(NA,chi),
                        ChiDf = c(NA,df), pvalue = c(NA,pvalue))
  row.names(table_sd) <- c("Object","..1")
  print(as.matrix(table_sd), right = FALSE, digits = 4, na.print = "")


}
