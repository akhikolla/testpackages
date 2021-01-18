##################################################################
# gtTest() function                                              #
##################################################################

#' @title Hypothesis test for one proportion in group testing
#' 
#' @description Calculates p-values for hypothesis tests of 
#' single proportions estimated from group testing 
#' experiments against a threshold proportion 
#' in the hypotheses. Available methods include the exact test, 
#' score test, and Wald test.
#' 
#' @param n integer specifying the number of groups.
#' @param y integer specifying the number of positive groups.
#' @param s integer specifying the common size of groups.
#' @param p.hyp the hypothetical threshold proportion against which to test, 
#' specified as a number between 0 and 1.
#' @param alternative character string defining the alternative 
#' hypothesis, either \kbd{"two.sided"}, \kbd{"less"}, or \kbd{"greater"}.
#' @param method character string defining the test method to be 
#' used. Options include "exact" for an exact test corresponding 
#' to the Clopper-Pearson confidence interval, "score" for a score 
#' test corresponding to the Wilson confidence interval, and "Wald" 
#' for a Wald test corresponding to the Wald confidence interval. 
#' The Wald method is not recommended. The "exact" method uses 
#' \code{binom.test{stats}}.
#' 
#' @details This function assumes equal group sizes, no testing error 
#' (i.e., 100 percent sensitivity and specificity) to test the groups, and 
#' individual units randomly assigned to the groups with identical true 
#' probability of success.
#' 
#' @return A list containing:
#' \item{p.value}{the p-value of the test}
#' \item{estimate}{the estimated proportion}
#' \item{p.hyp}{the threshold proportion provided by the user.}
#' \item{alternative}{the alternative provided by the user.}
#' \item{method}{the test method provided by the user.}
#' 
#' @author This function was originally written as \code{bgtTest} by Frank 
#' Schaarschmidt for the \code{binGroup} package. Minor modifications have 
#' been made for inclusion of the function in the \code{binGroup2} package.
#' 
#' @seealso \code{\link{propCI}} for confidence intervals in 
#' group testing and \code{binom.test(stats)} for the 
#' exact test and corresponding confidence interval.
#' 
#' @family estimation functions
#' 
#' @examples
#' # Consider the following the experiment: Tests are 
#' #   performed on n=10 groups, each group has a size
#' #   of s=100 individuals. The aim is to show that 
#' #   less than 0.5 percent (\eqn{p < 0.005}) of the units
#' #   in the population show a detrimental trait (positive test).
#' #   y=1 positive test and 9 negative tests are observed.
#' gtTest(n=10, y=1, s=100, p.hyp=0.005, alternative="less", 
#'         method="exact")
#' 
#' # The exact test corresponds to the 
#' #   limits of the Clopper-Pearson confidence interval
#' #   in the example of Tebbs & Bilder (2004):
#' gtTest(n=24, y=3, s=7, alternative="two.sided",
#'         method="exact", p.hyp=0.0543)
#'         
#' gtTest(n=24, y=3, s=7, alternative="two.sided",
#'         method="exact", p.hyp=0.0038)
#'         
#' # Hypothesis test with a group size of 1.
#' gtTest(n=24, y=3, s=1, alternative="two.sided",
#'         method="exact", p.hyp=0.1)         
#'         
#' # Further methods:
#' gtTest(n=24, y=3, s=7, alternative="two.sided",
#'         method="score", p.hyp=0.0516)
#'         
#' gtTest(n=24, y=3, s=7, alternative="two.sided",
#'         method="Wald", p.hyp=0.0401)

# Brianna Hitt - 02.23.2020
# Changed the capitalization of CI methods to match propCI()

gtTest <-
  function(n, y, s, p.hyp, alternative="two.sided", method="exact"){
    
    if(length(n)!=1 || (n<1 | abs(round(n)-n) > 1e-07)){stop("number of groups n must be specified as a single integer > 0")}
    if(length(s)!=1 || (s<1 | abs(round(s)-s) > 1e-07)){stop("group size s must be specified as a single integer > 0")}
    if(length(s)!=1 || (y<0 | abs(round(y)-y) > 1e-07)){stop("observed number of positive groups y must be specified as a single integer>0")}
    if(y>n) {stop("number of positive tests y can not be greater than number of groups n")}
    if(length(p.hyp)!=1 || p.hyp<0 || p.hyp>1){stop("the proportion in the hypothesis p.hyp must be specified as a single value between 0 and 1")}
    
    method<-match.arg(method, choices=c("exact", "score", "Wald"))
    alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))
    
    estimate=1-(1-y/n)^(1/s)
    
    bgtsumProb<-function(x, n, s, p)
    {
      sumprob=0
      for(i in x)
      {
        sumprob = sumprob + choose(n=n,k=i) * ((1-(1-p)^s)^(i)) * (1-p)^(s*(n-i))
      }
      sumprob
    }
    
    switch(method,
           
           "exact"={
             if(alternative=="less")
             {p.val = bgtsumProb(x=0:y, n=n, s=s, p=p.hyp)}
             if(alternative=="greater")
             {p.val = bgtsumProb(x=y:n, n=n, s=s, p=p.hyp)}
             if (alternative=="two.sided")
             {p.val = min( 2*(bgtsumProb(x=0:y, n=n, s=s, p=p.hyp)), 2*(p.val = bgtsumProb(x=y:n, n=n, s=s, p=p.hyp)), 1)}
           },
           
           
           "Wald"={
             esti = 1-(1-y/n)^(1/s)
             varesti = (1-(1-esti)^s)/(n*(s^2)*(1-esti)^(s-2))  # variance estimator (see Swallow, 1985)
             teststat = (esti-p.hyp)/sqrt(varesti) 
             
             if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}
             if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 
             if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)} 
           },
           
           "score"={
             esti = y/n
             t.hyp = 1-(1-p.hyp)^s
             teststat = (esti-t.hyp)/(sqrt(t.hyp*(1-t.hyp)/n))
             
             if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}
             if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 
             if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)} 
           })
    
    out<-list(p.value=p.val,
              estimate=estimate,
              alternative=alternative,
              p.hyp=p.hyp,
              method=method)
    
    class(out) <- "gtTest"
    out
  }




##################################################################
# print.gtTest() function                                        #
##################################################################

#' @title Print method for objects of class "gtTest"
#' 
#' @description Print method for objects of class "gtTest" created 
#' by the \code{\link{gtTest}} function.
#' 
#' @param x An object of class "gtTest" (\code{\link{gtTest}}).
#' @param ... Additional arguments to be passed to \code{print}. 
#' Currently only \code{digits} to be passed to \code{signif} for 
#' appropriate rounding.
#' 
#' @return A print out of the p-value and point estimate resulting 
#' from \code{\link{gtTest}}.
#' 
#' @author This function was originally written as \code{print.bgtTest} by 
#' Brad Biggerstaff for the \code{binGroup} 
#' package. Minor modifications were made for inclusion of the function in 
#' the \code{binGroup2} package.

"print.gtTest" <-
  function(x, ...)
    
  {
    args <- list(...)
    if (is.null(args$digits)) {digits <- 4}
    else{digits <- args$digits}
    
    if (x$alternative == "two.sided") {alt.hyp = "true proportion is not equal to"}
    if (x$alternative == "less") {alt.hyp = "true proportion is less than"}
    if (x$alternative == "greater") {alt.hyp = "true proportion is greater than"}
    
    cat("\n")
    cat(x$method, "test for one proportion in group testing\n")
    cat("Alternative hypothesis:",alt.hyp,x$p.hyp,"\n", sep = " ")
    cat("p-value","=",signif( x$p.value, digits),"\n", sep = " ")
    cat("point estimate","=", signif( x$estimate, digits),"\n", sep = " ")
    invisible(x)
  }



