
#' Structural break test
#' 
#' This function performs structural break test that is based on p-variation. 
#' 
#' Lets \code{x} be a data that should be tested of structural breaks. 
#' Then the p-variation of the \code{BridgeT(x)} with \code{p=4} is the test's statistics.
#' 
#' The quantiles of H0 distribution is based on Monte-Carlo simulation of 140 millions iterations. 
#' The test is reliable then \code{length(x)} is between 100 and 10000.
#' The test might work with other lengths too, but it is not tested well.
#' The test will not compute then \code{length(x)<20}.
#' 
#' @return
#' If \code{FullInfo=TRUE} then function returns an object of the class \code{PvarBreakTest}. 
#' It is the \code{list} that contains:
#' \item{Stat}{a value of statistics (p-variation of transformed data).}
#' \item{CriticalValue}{the critical value of the test according to significant level.} 
#' \item{alpha}{the significant level.} 
#' \item{p.value}{approximate p-value.} 
#' \item{reject}{\code{logical}. If \code{TRUE}, the H0 was rejected.} 
#' \item{dname}{the name of data vector.} 
#' \item{p}{the power in p-variation calculus. The test performs only with the \code{p=4}.} 
#' \item{x}{a vector of original data.}
#' \item{y}{a vector of transformed data (\code{y=BridgeT(x)}).}
#' \item{Timelabel}{time label of \code{x}. Used only for ploting.}
#' \item{BreakPoints}{the indexes of break points suggestion.}
#' \item{Partition}{a vector of indexes that indicates the partition of \code{y} that achieves the p-variation maximum.}
#' @author Vygantas Butkus <Vygantas.Butkus@@gmail.com>
#' @references 
#' The test was proposed by A. Rackaskas. The test is based on the results given in the flowing article
#' 
#' [1] R. Norvaisa, A. Rackauskas. Convergence in law of partial sum processes in p-variation norm. 
#' Lth. Math. J., 2008., Vol. 48, No. 2, 212-227.
#' 
#' @seealso Tests statistics is  \code{\link{pvar}} of the data \code{BridgeT(x)}(see  \code{\link{BridgeT}}) with (p=4).
#' The critical value and the approximate  p-value of the test might by found by functions
#' \code{\link{PvarQuantile}} and  \code{\link{PvarPvalue}}.  
#' 
#' @param x a numeric vector of data values or an object of class \code{pvar}.
#' @param TimeLabel numeric, a time index of \code{x}. Used only for plotting.
#' @param alpha a small number greater then 0. It indicates the significant level of the test.
#' @param FullInfo \code{logical}. If \code{TRUE} (the default) the function will return an object of the class \code{PvarBreakTest} 
#' that saves all useful information. Otherwise only the statistics will by returned.
#' @export
#' @examples
#' set.seed(1)
#' MiuDiff <- 0.3
#' x <- rnorm(250*4, rep(c(0, MiuDiff, 0, MiuDiff), each=250))
#' 
#' plot(x, pch=19, cex=0.5, main='original data, with several shifts of mean')
#' k <- 50
#' moveAvg <- filter(x, rep(1/k, k))
#' lines(time(x), moveAvg, lwd=2, col=2)
#' legend('topleft', c('sample', 'moving average (k='%.%k%.%')'),
#'        lty=c(NA,1), lwd=c(NA, 2), col=1:2, pch=c(19,NA), pt.cex=c(0.7,1)
#'        ,inset = .03, bg='antiquewhite1')
#'
#' xtest <- PvarBreakTest(x)
#' plot(xtest)
PvarBreakTest <- function(x, TimeLabel = as.vector(time(x)), alpha = 0.05, FullInfo = TRUE) {
  
  dname <- deparse(substitute(x))
  NAInd <- is.na(x)
  if (any(NAInd)) {
    warning("NA values was removed.")
    x <- x[!NAInd]
    TimeLabel <- TimeLabel[!NAInd]
  }
  n <- length(x)
  if (n < 20) 
    stop("The size must be greater than 20.")
  if (n < 100) 
    warning("The test might by biased when n<100.")
  CriticalValue <- PvarQuantile(n, prob = 1 - alpha)
  
  if (sd(x) == 0) {
    y <- x
  } else {
    y <- BridgeT(x)
  }
  
  PvarStat <- pvar(y, p = 4, TimeLabel = TimeLabel)
  if (FullInfo) {
    Stat <- unname(PvarStat$value)
    p.value <- PvarPvalue(n, Stat)
    reject <- Stat >= CriticalValue[1]
    if (reject) {
      BreakPos <- abs(PvarStat$partition/n - 0.5)
      accept <- abs(PvarStat$partition/n - 0.5) < 0.4
      if (!any(accept)) {
        accept[which.min[BreakPos]] <- TRUE
      }
      BreakPoints <- PvarStat$partition[accept]
    } else {
      BreakPoints <- NULL
    }
    ans <- list(Stat = c(Statistics = Stat), CriticalValue = c(`Critical Value` = CriticalValue), alpha = c(alpha = alpha), p.value = c(`~p.value` = p.value), 
      reject = reject, dname = dname, p = unname(PvarStat$p), x = x, y = y, TimeLabel = TimeLabel, BreakPoints = BreakPoints, 
      Partition = PvarStat$partition)
    
    class(ans) <- "PvarBreakTest"
    return(ans)
  } else {
    Stat <- PvarStat
    p.value <- PvarPvalue(n, Stat)
    ans <- c(Statistics = Stat, `Critical Value` = CriticalValue, alpha = alpha, `~p.value` = p.value)
    return(ans)
  }
}


#' @rdname PvarBreakTest
#' @param main1 the \code{main} parameter of the data graph.
#' @param main2 the \code{main} parameter of the Bridge transformation graph.
#' @param ylab1 the \code{ylab} parameter of the data graph.
#' @param ylab2 the \code{ylab} parameter of the Bridge transformation graph.
#' @param sub2 the \code{sub} parameter of the Bridge transformation graph. By default it reports the number of break points.
#' @param cex.DP the cex of data points.
#' @param col.PP the color of partition points.
#' @param cex.PP the cex of partition points.
#' @param col.BP the color of break points.
#' @param cex.BP the cex of break points.
#' @param \dots further arguments, passed to \code{print}. 
#' @export
plot.PvarBreakTest <- function(x, main1 = "Data", main2 = "Bridge transformation", ylab1 = x$dname, ylab2 = "BridgeT(" %.% x$dname %.% 
  ")", sub2 = NULL, col.PP = 3, cex.PP = 0.5, col.BP = 2, cex.BP = 1, cex.DP = 0.5, ...) {
  Time <- x$TimeLabel
  op <- par(mfrow = c(2, 1), mar = c(5.1, 4.1, 2.1, 2.1))
  
  plot(Time, x$x, type = "p", pch = 19, cex = cex.DP, ylab = ylab1, main = main1, ...)
  BP <- c(1, x$BreakPoints, length(x$x))
  for (i in 2:length(BP)) {
    Eseg <- mean(x$x[BP[i - 1]:BP[i]])
    colseg <- i%%2 + 2
    graphics::segments(x0 = Time[BP[i - 1]], y0 = Eseg, x1 = Time[BP[i]], y1 = Eseg, lwd = 3, col = colseg)
  }
  
  if (is.null(sub2)) {
    if (x$reject) {
      if (length(x$BreakPoints) == 1) {
        sub2 <- "Program suggests " %.% length(x$BreakPoints) %.% " break point."
      } else {
        sub2 <- "Program suggests " %.% length(x$BreakPoints) %.% " break points."
      }
    } else {
      sub2 <- "Program didn't find structural breaks at the confidence level of " %.% formatC(utils::head(x$alpha, 1)) %.% "."
    }
  }
  
  plot(Time, x$y, type = "l", ylab = ylab2, main = main2, sub = sub2, ...)
  points(x$TimeLabel[x$Partition], x$y[x$Partition], cex = cex.PP, pch = 19, col = col.PP, bg = col.PP)
  points(x$TimeLabel[x$BreakPoints], x$y[x$BreakPoints], cex = cex.BP, pch = 19, col = col.BP, bg = col.BP)
  par(op)
}

#' @rdname PvarBreakTest
#' @param object the object of the class \code{PvarBreakTest}.
#' @method summary PvarBreakTest
#' @export
summary.PvarBreakTest <- function(object, ...) {
  class(object) <- c("summary.PvarBreakTest", "PvarBreakTest")
  object
}

#' @method print PvarBreakTest
#' @export
print.PvarBreakTest <- function(x, ...) {
  cat("       PvarBreakTest \n\n")
  cat("H0: no structural change \n")
  cat("Results: ")
  if (x$reject) {
    cat("H0 is rejected at the confidence level of " %.% formatC(utils::head(x$alpha, 1)) %.% ".\n")
  } else {
    cat("H0 is accepted at the confidence level of " %.% formatC(utils::head(x$alpha, 1)) %.% ".\n")
  }
  cat("Data: " %.% x$dname %.% ", n=" %.% length(x$x) %.% ".\n")
  cat("The output of the test: \n")
  print(c(x$Stat, x$CriticalValue, x$alpha, x$p.value))
}

#' @method print summary.PvarBreakTest
#' @export
print.summary.PvarBreakTest <- function(x, ...) {
  cat("The summary of PvarBreakTest:\n")
  cat("H0: no structural change. \n")
  cat("Results: ")
  if (x$reject) {
    cat("H0 is rejected at the confidence level of " %.% formatC(utils::head(x$alpha, 1)) %.% ".\n")
  } else {
    cat("H0 is accepted at the confidence level of " %.% formatC(utils::head(x$alpha, 1)) %.% ".\n")
  }
  if (x$reject) {
    if (length(x$BreakPoints) > 6) {
      cat("Suggesting " %.% length(x$BreakPoints) %.% " break points: " %.% paste(formatC(utils::head(x$BreakPoints, 6)), collapse = ", ") %.% 
        ", ...\n")
    } else {
      if (length(x$BreakPoints) == 1) {
        cat("Suggesting " %.% length(x$BreakPoints) %.% " break point: " %.% paste(formatC(utils::head(x$BreakPoints, 6)), collapse = ", ") %.% 
          ".\n")
      } else {
        cat("Suggesting " %.% length(x$BreakPoints) %.% " break points: " %.% paste(formatC(utils::head(x$BreakPoints, 6)), collapse = ", ") %.% 
          ".\n")
      }
    }
  }
  cat("Data: " %.% x$dname %.% ", n=" %.% length(x$x) %.% ".\n")
  cat("Test's output: \n")
  print(c(x$Stat, x$CriticalValue, x$alpha, x$p.value))
  cat("p-avriation calculation info:\n")
  print(unlist(x$info)[-length(x$info)])
  if (length(x$x) > 6) {
    cat("\nData vector (n=" %.% length(x$x) %.% "): " %.% paste(formatC(utils::head(x$x, 6)), collapse = ", ") %.% ", ...\n")
  } else {
    cat("\nData vector (n=" %.% length(x$x) %.% "): " %.% paste(formatC(utils::head(x$x, 6)), collapse = ", ") %.% ".\n")
  }
}

############################################################################################################ 

#' Bridge transformation
#' 
#' Transforms data by Bridge transformation.
#' 
#' Let \code{n} denotes the length ox \code{x}.
#' For  each \eqn{m \in [1,n]} bridge transformations \code{BridgeT}
#' is defined as  
#' 
#' \deqn{
#'   BridgeT(m, x) = \left\{ \sum_{i=1}^m x_i - \frac{m}{n} \sum_{i=1}^n x_i  \right\} .
#' }{
#'   BridgeT(m, x) = ( \sum^m x_i - m/n \sum^n x_i ).
#' }
#' 
#' Meanwhile, the transformation with normalization is 
#' 
#' \deqn{
#'   BridgeT(m, x) = \frac{1}{\sqrt{n var(x)}} \left\{ \sum_{i=1}^m x_i - \frac{m}{n} \sum_{i=1}^n x_i  \right\} .
#' }{
#'   BridgeT(m, x) = ( \sum^m x_i - m/n \sum^n x_i  ) / (n var(x))^0.5.
#' }
#' @return A numeric vector.
#' @seealso \code{\link{PvarBreakTest}},  \code{\link{rbridge}}
#' 
#' @param x x a numeric vector of data values.
#' @param normalize \code{logical}, indicating whether the vector should be normalized. 
#' @export
#' @examples
#' x <- rnorm(1000)
#' Bx <- BridgeT(x, FALSE)
#' 
#' op <- par(mfrow=c(2,1),mar=c(4,4,2,1))
#' plot(cumsum(x), type="l")
#' plot(Bx, type="l")
#' par(op)
BridgeT <- function(x, normalize = TRUE) {
  if (normalize) {
    (cumsum(x) - seq_along(x)/length(x) * sum(x))/sqrt(length(x) * var(x))
  } else {
    (cumsum(x) - seq_along(x)/length(x) * sum(x))
  }
}
################################################################################################################ 

#' Quantiles and probabilities of p-variation
#' 
#' The distribution of p-variation of \code{BridgeT(x)} depends on \code{n=length(x)}. 
#' This fact is important for getting appropriate quantiles (or p-value).
#' These functions helps to deal with it. 
#' 
#' The distribution of p-variance is form Monte-Carlo simulation based on 140 millions iterations. 
#' The data frame \code{\link{PvarQuantileDF}} saves the results of Monte-Carlo simulation.
#' 
#' Meanwhile, \code{MeanCoef} and \code{SdCoef} defines the coefficients of functional 
#' form (conditional on \code{n}) of \code{mean} and \code{sd} statistics.
#' 
#' A functional form of \code{mean} and \code{sd} statistics are the same, namely
#' \deqn{
#'   f(n) = b_1 + b_2 n^b_2 .
#' }{
#'   f(n) = b_1 + b_2 * n^b_2 .
#' }
#' 
#' The coefficients \eqn{(b_1, b_2, b_3)} are saved in vectors \code{MeanCoef} and \code{SdCoef}.
#' Those vectors are estimated with \code{nls} function form Monte-Carlo simulation.
#' @rdname StatisticsPvarBreakTest
#' @return Functions \code{PvarQuantile} and \code{PvarPvalue} returns a corresponding value quantile or the probability.
#' Functions \code{getMean} and \code{getSd} returns a corresponding value of \code{mean} and \code{sd} statistics.
#' Function \code{NormalisePvar} returns normalize values. 
#' @note  Arguments \code{n}, \code{stat} and \code{prob} might be vectors,
#' but they can't be vectors simultaneously (at least one of then must be a number).
#' @seealso  \code{\link{PvarBreakTest}}, \code{\link{PvarQuantileDF}},  
#' \code{\link{NormalisePvar}},  \code{\link{getMean}}, \code{\link{getSd}}
#' 
#' @param n a positive integer indicating the length of data vector.
#' @param prob cumulative probabilities of p-variation distribution.
#' @param DF a \code{data.frame} that links \code{prob} and \code{stat} .
#' @export
PvarQuantile <- function(n, prob = c(0.9, 0.95, 0.99), DF = PvarQuantileDF) {
  intervals <- cut(prob, breaks = DF$prob, include.lowest = TRUE, right = TRUE)
  Quant <- DF$Quant[as.numeric(intervals) + 1]
  if (any(is.na(Quant))) 
    warning("`prob` must be between 0 and 1")
  ans <- Quant * getSd(n) + getMean(n)
  unname(ans)
}

#' @rdname StatisticsPvarBreakTest 
#' @param stat a vector of p-variation statistics.
#' @export
PvarPvalue <- function(n, stat, DF = PvarQuantileDF) {
  NormStat <- NormalisePvar(stat, n)
  intervals <- cut(NormStat, breaks = c(-Inf, DF$Quant, Inf), include.lowest = TRUE, right = FALSE)
  ind <- as.numeric(intervals) - 1
  ind[ind == 0] <- 1
  1 - DF$prob[ind]
}

#' @rdname StatisticsPvarBreakTest 
#' @param bMean a coefficient vector that defines a function of the mean of p-variation.
#' @export
### get the mean of p-variation according of BridgeT(x), then H0 is TRUE (according to n).
getMean <- function(n, bMean = MeanCoef) {
  unname(bMean[1] + bMean[2] * n^bMean[3])
}

#' @rdname StatisticsPvarBreakTest 
#' @param bSd a coefficient vector that defines a function of the standard deviation of p-variation.
#' @export
### get the standart deviation of p-variation according of BridgeT(x), then H0 is TRUE (according to n).
getSd <- function(n, bSd = SdCoef) {
  unname(bSd[1] + bSd[2] * n^bSd[3])
}

#' @rdname StatisticsPvarBreakTest 
#' @param x a numeric vector of data values.
#' @export
### Normalyse p-variation according n.
NormalisePvar <- function(x, n, bMean = MeanCoef, bSd = SdCoef) {
  (x - getMean(n, bMean))/getSd(n, bSd)
} 
