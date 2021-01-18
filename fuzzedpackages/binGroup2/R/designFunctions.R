##################################################################
# nDesign() function                                             #
##################################################################

# Brianna Hitt - 02.24.2020
# Changed capitalization of methods to match propCI()

nDesign <-
  function(nmax, s, delta, p.hyp, conf.level = 0.95, power = 0.8, 
           alternative = "two.sided", method = "CP", biasrest = 0.05)
  {
    
    if (length(nmax) < 1 || length(nmax) > 2 || (min(nmax) <= 3 | abs(round(nmax) - nmax) > 1e-07))
    {stop("the maximal number of groups n allowed in calculations must be one or two integer(s) greater than 1")}
    if (length(s) != 1 || (s < 1 | abs(round(s) - s) > 1e-07))
    {stop("group size s must be specified as a single integer>0")}
    if (length(conf.level) != 1 || conf.level < 0 || conf.level > 1)
    {stop("conf.level must be a positive number between 0 and 1")}
    if (length(power) != 1 || power < 0 || power > 1)
    {stop(" desired power must be a positive number between 0 and 1, f.e. 0.8 for rejecting H0 in 80% of the cases")}
    
    method <- match.arg(method, choices = c("CP", "Blaker", "AC", "score", "Wald", "soc"))
    
    alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
    
    if (length(p.hyp) != 1 || p.hyp > 1 || p.hyp < 0)
    {stop("true proportion p.hyp must be specified as a single number between 0 and 1")}
    
    if (length(delta) != 1)
    {stop("delta must be specified as a single number")}
    if (alternative == "less")
    {
      if (p.hyp - delta < 0 || p.hyp - delta > 1 )
      {stop("alternative=less: specify delta as a number between 0 and the threshold p.hyp")}
    }
    
    if (alternative == "greater")
    {
      if (p.hyp + delta < 0 || p.hyp + delta > 1 )
      {stop("alternative=greater: specify delta as a number between the threshold p.hyp and 1")}
    }
    
    if (alternative == "two.sided")
      
    {
      if (p.hyp + delta < 0 || p.hyp + delta > 1 || p.hyp - delta < 0 || p.hyp - delta > 1)
      {stop("alternative=two.sided: specify delta as a number between the threshold p.hyp and 1")}
    }
    
    if (length(biasrest) != 1 || biasrest >= 1 || biasrest < 0)
    {stop("the maximally allowed bias(p) specified in biasrest must be a single number between 0 and 1, usually should be close to 0")}
    
    # # # # # #
    
    
    if (length(nmax) == 1)
    {
      nit <- 4:nmax
      powerit <- numeric(length = length(nit))
      biasit <- numeric(length = length(nit))
      
      for (i in 1:length(nit))
      {
        temp <- bgtPowerI(n = nit[i], s = s, delta = delta, p.hyp = p.hyp, conf.level = conf.level, alternative = alternative, method = method)
        powerit[i] <- temp$power
        biasit[i] <- temp$bias
        if (temp$bias <= biasrest && temp$power >= power)
        {
          out <- list(nout = nit[i], powerout = powerit[i], biasout = temp$bias, power.reached = TRUE, bias.reached = FALSE, biasit = biasit,
                    nit = nit, powerit = powerit, delta = delta, p.hyp = p.hyp, power = power, biasrest = biasrest, alternative = alternative, maxit = i)
          class(out) <- "nDesign"
          return(out)
        }
      }
      
      # # # bias decreases monotone with increasing n: if nmax has bias> biasrest: all designs have bias> biasrest
      lastn <- length(nit)
      if (biasit[lastn] > biasrest)
      {
        
        out <- list(nout = nmax, powerout = powerit[lastn], biasout = biasit[lastn], power.reached = FALSE, bias.reached = TRUE,biasit = biasit,
                  nit = nit, powerit = powerit,delta = delta,p.hyp = p.hyp, power = power, biasrest = biasrest, alternative = alternative, maxit = lastn)
        class(out) <- "nDesign"
        return(out)
      }
      
      npowmax <- nit[which.max(powerit)]
      powout <- powerit[which.max(powerit)]
      biasout <- biasit[which.max(powerit)]
      {
        out <- list(nout = npowmax, powerout = powout, biasout = biasout, power.reached = FALSE, bias.reached = FALSE,biasit = biasit,
                  nit = nit, powerit = powerit,delta = delta,p.hyp = p.hyp, power = power, biasrest = biasrest, alternative = alternative, maxit = length(nit))
        class(out) <- "nDesign"
        return(out)
        
      }
      
    }
    
    
    
    if (length(nmax) == 2) {
      
      nfrom <- min(nmax)
      nto <- max(nmax)
      if (nfrom < 4 && method == "soc") {stop("The 'soc' interval can have bounds NaN for n < 4")}
      
      nit <- nfrom:nto
      powerit <- numeric(length = length(nit))
      biasit <- numeric(length = length(nit))
      
      for (i in 1:length(nit))
      {
        temp <- bgtPowerI(n = nit[i], s = s, delta = delta, p.hyp = p.hyp, conf.level = conf.level, alternative = alternative, method = method)
        powerit[i] <- temp$power
        biasit[i] <- temp$bias
        if (temp$bias <= biasrest && temp$power >= power)
        {
          out <- list(nout = nit[i], powerout = powerit[i], biasout = temp$bias, power.reached = TRUE, bias.reached = FALSE,biasit = biasit,
                    nit = nit, powerit = powerit,delta = delta,p.hyp = p.hyp, power = power, biasrest = biasrest, alternative = alternative, maxit = i )
          class(out) <- "nDesign"
          return(out)
        }
      }
      
      # # # bias decreases monotone with increasing n: if nmax has bias> biasrest: all designs have bias> biasrest
      lastn <- length(nit)
      if (biasit[lastn] > biasrest)
      {
        out <- list(nout = max(nmax), powerout = powerit[lastn], biasout = biasit[lastn], power.reached = FALSE, 
                    biasit = biasit, bias.reached = TRUE,nit = nit, powerit = powerit,delta = delta,p.hyp = p.hyp, 
                    power = power, biasrest = biasrest, alternative = alternative, maxit = lastn)
        class(out) <- "nDesign"
        return(out)
      }
      
      npowmax <- nit[which.max(powerit)]
      powout <- powerit[which.max(powerit)]
      biasout <- biasit[which.max(powerit)]
      {
        out <- list(nout = npowmax, powerout = powout, biasout = biasout, power.reached = FALSE, bias.reached = FALSE, biasit = biasit,
                  nit = nit, powerit = powerit, delta = delta, p.hyp = p.hyp, power = power, biasrest = biasrest, alternative = alternative, maxit = length(nit))
        class(out) <- "nDesign"
        return(out)
      }
    }
  }




##################################################################
# sDesign() function                                             #
##################################################################

sDesign <-
  function(n,smax,delta,p.hyp,conf.level=0.95, power=0.8, alternative="two.sided", method="CP", biasrest=0.05)
    
  {
    
    if (length(smax) != 1 || (smax < 3 | abs(round(smax) - smax) > 1e-07))
    {stop("the maximal group size smax allowed in calculations must be a single integer greater than 0")}
    if (length(n) != 1 || (n <= 1 | abs(round(n) - n) > 1e-07))
    {stop("the number of groups n must be specified as a single integer>1")}
    if (length(conf.level) != 1 || conf.level < 0 || conf.level > 1)
    {stop("conf.level must be a positive number between 0 and 1")}
    if (length(power) != 1 || power < 0 || power > 1)
    {stop("desired power must be a positive number between 0 and 1, f.e. 0.8 for rejecting H0 in 80% of the cases")}
    
    method <- match.arg(method, choices = c("CP", "Blaker", "AC", "score", "Wald", "soc"))
    
    alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
    
    if (length(p.hyp) != 1 || p.hyp > 1 || p.hyp < 0)
    {stop("threshold p.hyp must be specified as a single number between 0 and 1")}
    
    if (length(delta) != 1)
    {stop("delta must be specified as a single number")}
    if (alternative == "less")
    {
      if (p.hyp - delta < 0 || p.hyp - delta > 1)
      {stop("alternative=less: specify delta as a number between 0 and the threshold p.hyp")}
    }
    
    if (alternative == "greater")
    {
      if (p.hyp + delta < 0 || p.hyp + delta > 1)
      {stop("alternative=greater: specify delta as a number between the threshold p.hyp and 1")}
    }
    
    if (alternative == "two.sided")
    {
      if (p.hyp + delta < 0 || p.hyp + delta > 1 || p.hyp - delta < 0 || p.hyp - delta > 1)
      {stop("alternative=two.sided: specify delta as a number between the threshold p.hyp and 1")}
    }
    
    if (length(biasrest) != 1 || biasrest >= 1 || biasrest < 0)
    {stop("the maximally allowed bias(p) specified in biasrest must be a single number between 0 and 1, usually should be close to 0")}
    
    
    # # # Iteration until smax, until either the desired power is reached or biasrestriction is violated
    
    
    if (method == "soc" && n <= 3) {stop("number of groups n<=3 might cause problems in computation of 'soc' interval")}
    
    
    sit <- 2:smax
    powerit <- numeric(length = length(sit))
    biasit <- numeric(length = length(sit))
    
    for (i in 1:length(sit))
    {
      
      temp <- bgtPowerI(n = n, s = sit[i], delta = delta, p.hyp = p.hyp, conf.level = conf.level, alternative = alternative, method = method)
      powerit[i] <- temp$power
      biasit[i] <- temp$bias
      
      if (temp$bias <= biasrest && temp$power >= power)
      {
        out <- list(sout = sit[i], powerout = powerit[i], biasout = biasit[i],
                    power.reached = TRUE, bias.reached = FALSE, powerit = powerit, biasit = biasit, sit = sit, maxit = i,
                    alternative = alternative, p.hyp = p.hyp, delta = delta, biasrest = biasrest, power = power)
        
        class(out) <- "sDesign"
        return(out)
      }
      
      if (temp$bias > biasrest)
      {
        out <- list(sout = sit[which.max(powerit[1:(i - 1)])],
                    powerout = powerit[which.max(powerit[1:(i - 1)])],
                    biasout = biasit[which.max(powerit[1:(i - 1)])],
                    power.reached = FALSE, bias.reached = TRUE,
                    powerit = powerit,biasit = biasit,sit = sit, maxit = i,
                    alternative = alternative, p.hyp = p.hyp, delta = delta, biasrest = biasrest, power = power)
        
        class(out) <- "sDesign"
        return(out)
      }
    }
    ## end of for statement
    
    out <- list(sout = sit[which.max(powerit)],
                powerout = powerit[which.max(powerit)],
                biasout = biasit[which.max(powerit)],
                power.reached = FALSE, bias.reached = FALSE,
                powerit = powerit, biasit = biasit,sit = sit,maxit = length(sit),
                alternative = alternative, p.hyp = p.hyp, delta = delta, biasrest = biasrest, power = power )
    
    class(out) <- "sDesign"
    return(out)
  }




##################################################################
# designPower() function                                         #
##################################################################

#' @title Number of groups or group size needed to achieve a power level in 
#' one parameter group testing
#' 
#' @description For a fixed number of groups (group size), determine the 
#' group size (number of groups) needed to obtain a specified power level to 
#' reject a hypothesis for a proportion in one parameter group testing.
#' 
#' @param n integer specifying the maximum number of groups \kbd{n} allowed 
#' when \kbd{fixed="s"} or the fixed number of groups when \kbd{fixed="n"}. 
#' When \kbd{fixed="s"}, a vector of two integers giving the range of \kbd{n} 
#' which power shall be iterated over is also allowed.
#' @param s integer specifying the fixed group size (number of units per group) 
#' when \kbd{fixed="s"} or the maximum group size allowed in the planning of 
#' the design when \kbd{fixed="n"}.
#' @param fixed character string specifying whether the number of groups 
#' \kbd{"n"} or the group size \kbd{"s"} is to be held at a fixed value.
#' @param delta the absolute difference between the true proportion and the 
#' hypothesized proportion which shall be detectable with the specified power.
#' @param p.hyp the proportion in the hypotheses, specified as a value between 
#' 0 and 1.
#' @param conf.level confidence level of the decision. The default 
#' confidence level is 0.95.
#' @param power level of power to be achieved, specified as a 
#' probability between 0 and 1.
#' @param alternative character string defining the alternative hypothesis, 
#' either \kbd{"two.sided"}, \kbd{"less"}, or \kbd{"greater"}.
#' @param method character string specifying the confidence interval method 
#' (see \code{\link{propCI}}) to be used.
#' @param biasrest a value between 0 and 1, specifying the absolute bias
#' maximally allowed for a point estimate.
#' 
#' @details The power of a hypothesis test performed by a confidence interval 
#' is defined as the probability that a confidence interval excludes the 
#' thresholdparameter (\kbd{p.hyp}) of the hypothesis.
#' 
#' When \kbd{fixed="s"}, this function increases the number of groups until a 
#' pre-specified level of power is reached or the maximum number of groups 
#' \kbd{n} is reached. Since the power does not increase monotonely with 
#' increasing \kbd{n} for single proportions but oscillates between local 
#' maxima and minima, the simple iteration given here will generally result in 
#' selecting \kbd{n} for which the given confidence interval method shows a 
#' local minimum of coverage if the null hypothesis is true. Bias decreases 
#' monotonely with increasing the number of groups (if other parameters are 
#' fixed). The resulting problems of choosing a number of groups which results 
#' in satisfactory power are solved in the following manner:
#' 
#' In the case that the pre-specified power is reached within the given 
#' range of \kbd{n}, the smallest \kbd{n} is returned for which at least 
#' this power is reached, as well as the actual power for this \kbd{n}.
#' 
#' In the case that the pre-specified power is not reached within the given 
#' value, that \kbd{n} is returned for which maximum power is achieved, and 
#' the corresponding value of power.
#' 
#' In the case that the bias restriction is violated even for the largest 
#' \kbd{n} within the given range of \kbd{n}, simply that \kbd{n} will be 
#' returned for which power was largest in the given range.
#' 
#' Especially for large \kbd{n}, the calculation time may become large 
#' (particularly for the Blaker interval). Alternatively, the function 
#' \code{\link{gtPower}} might be used instead to calculate power and bias 
#' only for some particular combinations of the input arguments.
#' 
#' When \kbd{fixed="n"}, this function increases the size of groups until a 
#' pre-specified level of power is reached. Since the power does not increase 
#' monotonely with increasing \kbd{s} for single proportions but oscillates 
#' between local maxima and minima, the simple iteration given here will 
#' generally result in selecting \kbd{s} for which the given confidence 
#' interval method shows a local minimum of coverage if the null hypothesis 
#' is true. Since the positive bias of the estimator in group testing 
#' increases with increasing group size, this function checks whether the bias 
#' is smaller than a pre-specified level (\kbd{bias.rest}). If the bias violates 
#' this restriction for a given combination \kbd{n}, \kbd{s}, and \kbd{delta}, 
#' \kbd{s} will not be further increased and the actual power of the last 
#' acceptable group size \kbd{s} is returned.
#' 
#' @return A list containing:
#' \item{nout}{the number of groups necessary to reach the power with the 
#' specified parameters, when \kbd{fixed="s"} only.}
#' \item{sout}{the group size necessary to meet the conditions, when 
#' \kbd{fixed="n"} only.}
#' \item{powerout}{the power for the specified parameters and the selected 
#' number of groups \kbd{n} when \kbd{fixed="s"} or the selected group 
#' size \kbd{s} when \kbd{fixed="n"}.}
#' \item{biasout}{the bias for the specified parameters and the selected 
#' number of groups \kbd{n} when \kbd{fixed="s"} or the selected group 
#' size \kbd{s} when \kbd{fixed="n"}.}
#' \item{power.reached}{a logical value indicating whether the 
#' specified level of power was reached.}
#' \item{bias.reached}{a logical value indicating whether the maximum 
#' allowed bias was reached.}
#' \item{nit}{the number of groups for each iteration.}
#' \item{sit}{the group size for each iteration.}
#' \item{powerit}{the power achieved for each iteration.}
#' \item{biasit}{the bias for each iteration.}
#' \item{maxit}{the iteration at which the maximum power was reached, 
#' or the total number of iterations.}
#' \item{alternative}{the alternative hypothesis specified by the user.}
#' \item{p.hyp}{the hypothesized proportion specified by the user.}
#' \item{delta}{the absolute difference between the true proportion and the 
#' hypothesized proportion specified by the user.}
#' \item{power}{the desired power specified by the user.}
#' \item{biasrest}{the maximum absolute bias specified by the user.}
#' 
#' @author The \code{nDesign} and \code{sDesign} functions were originally 
#' written by Frank Schaarschmidt for the \code{binGroup} package. Minor 
#' modifications were made for inclusion in the \code{binGroup2} package.
#' 
#' @references 
#' \insertRef{Swallow1985}{binGroup2}
#' 
#' @seealso 
#' \code{\link{gtPower}} for calculation of power and bias depending 
#' on \kbd{n}, \kbd{s}, \kbd{delta}, \kbd{p.hyp}, \kbd{conf.level}, 
#' and \kbd{method}, and \code{\link{designEst}} to choose the group size 
#' \kbd{s} according to the minimal mse of the estimator, as given in 
#' Swallow (1985).
#' 
#' @family estimation functions
#' 
#' @examples 
#' # Assume the objective is to show that a proportion is 
#' #   smaller than 0.005 (i.e. 0.5 percent) with a power 
#' #   of 0.80 (i.e. 80 percent) if the unknown proportion
#' #   in the population is 0.003 (i.e. 0.3 percent);
#' #   thus, a delta of 0.002 shall be detected.
#' # A 95% Clopper Pearson CI shall be used. 
#' # The maximum group size because of limited
#' #   sensitivity of the diagnostic test might be s=20 and we
#' #   can only afford to perform maximally 100 tests:
#' designPower(n=100, s=20, delta=0.002, p.hyp=0.005, fixed="s",
#'              alternative="less", method="CP", power=0.8)
#'         
#' # One might accept to detect delta=0.004,
#' #   i.e. reject H0: p>=0.005 with power 80 percent 
#' #   when the true proportion is 0.001:
#' designPower(n=100, s=20, delta=0.004, p.hyp=0.005, fixed="s",
#'              alternative="less", method="CP", power=0.8)
#'              
#' # Power for a design with a fixed group size of s=1 
#' #   (individual testing).
#' designPower(n=500, s=1, delta=0.05, p.hyp=0.10, 
#'             fixed="s", method="CP", power=0.80)
#'         
#' # Assume that objective is to show that a proportion
#' #   is smaller than 0.005 (i.e. 0.5%) with a 
#' #   power of 0.80 (i.e. 80%) if the unknown proportion
#' #   in the population is 0.003 (i.e. 0.3%); thus, a 
#' #   delta = 0.002 shall be detected. 
#' # A 95% Clopper-Pearson CI shall be used. 
#' # The maximum number of groups might be 30, where the 
#' #   overall sensitivity is not limited until group 
#' #   size s=100.
#' designPower(s=100, n=30, delta=0.002, p.hyp=0.005, fixed="n",
#'              alternative="less", method="CP", power=0.8)
#'         
#' # One might accept to detect delta=0.004,
#' #   i.e. reject H0: p>=0.005 with power 80 percent 
#' #   when the true proportion is 0.001:
#' designPower(s=100, n=30, delta=0.004, p.hyp=0.005, fixed="n",
#'              alternative="less", method="CP", power=0.8)
#' designPower(s=100, n=30, delta=0.004, p.hyp=0.005, fixed="n",
#'              alternative="less", method="score", power=0.8)

designPower <- function(n, s, fixed = "s", delta, p.hyp, 
                         conf.level = 0.95, power = 0.80, 
                         alternative = "two.sided", 
                         method="CP", biasrest = 0.05) {
  
  if (fixed == "s") {
    results <- nDesign(nmax = n, s = s, delta = delta,
                       p.hyp = p.hyp, conf.level = conf.level,
                       power = power, alternative = alternative,
                       method = method, biasrest = biasrest)

  } else if (fixed == "n") {
    results <- sDesign(n = n, smax = s, delta = delta,
                       p.hyp = p.hyp, conf.level = conf.level,
                       power = power, alternative = alternative,
                       method = method, biasrest = biasrest)
  }
  class(results) <- "designPower"
  return(results)
}




##################################################################
# print.designPower() function -                                 #
##################################################################

#' @title Print method for objects of class "designPower"
#' 
#' @description Print method for objects of class "designPower" 
#' created by \code{\link{designPower}}.
#' 
#' @param x an object of class "designPower" created by 
#' \code{\link{designPower}}.
#' @param ... additional arguments to be passed to \code{print}. 
#' Currently only \code{digits} to be passed to \code{signif} for 
#' appropriate rounding.
#' 
#' @return A print out detailing whether or not power was reached in 
#' the range of values (\kbd{n} or \kbd{s}) provided, the maximal 
#' power reached in the range of values, the alternative hypothesis, 
#' and the assumed true proportion.
#' 
#' @author This function was originally written as \code{print.bgtDesign} 
#' by Frank Schaarschmidt for the \code{binGroup} package. Minor 
#' modifications were made for inclusion in the \code{binGroup2} package.

"print.designPower" <- function(x, ...){
  
  args <- list(...)
  if (is.null(args$digits)) {digits <- 4}
  else{digits <- args$digits}
  
  if (x$alternative == "less")
  {alt.hyp <- paste("true proportion is less than",x$p.hyp )
  ptrue <- paste("assumed true proportion","=", x$p.hyp - x$delta)}
  
  if (x$alternative == "greater")
  {alt.hyp <- paste("true proportion is greater than", x$p.hyp )
  ptrue <- paste("assumed true proportion","=", x$p.hyp + x$delta)}
  
  if (x$alternative == "two.sided")
  {alt.hyp <- paste("true proportion is not equal to",x$p.hyp )
  ptrue <- paste("assumed true proportion","=", x$p.hyp - x$delta,"or", x$p.hyp + x$delta)}
  
  if (names(x)[1] == "nout") {
    
    if (x$power.reached == TRUE && x$bias.reached == FALSE)
    {cat("Power was reached without violating bias restriction\n")
      cat("for n","=",x$nout,"with power", "=",signif(x$powerout, digits),"\n")
      cat("and bias","=",signif(x$biasout, digits),"\n")
      
      cat("alternative hypothesis:", alt.hyp ,"\n")
      cat( ptrue, "\n")
    }
    
    if (x$power.reached == FALSE && x$bias.reached == FALSE && x$powerout != 0)
    {cat("Power was not reached in the range of n","=",min(x$nit),"to",max(x$nit),"\n")
      cat("Maximal power was reached for n","=",x$nout,"with power","=",signif(x$powerout, digits),"\n")
      cat("and bias","=",signif(x$biasout, digits),"\n")
      
      cat("alternative hypothesis:", alt.hyp ,"\n")
      cat( ptrue, "\n")
    }
    
    if (x$power.reached == FALSE && x$bias.reached == TRUE)
    {cat("Power can not be reached without violating bias restriction\n")
      cat("in the range of n","=",min(x$nit),"to",max(x$nit),"\n")
      cat("Maximal power was reached for n","=",x$nout,"\n")
      cat("with power","=",signif(x$powerout, digits),"and bias","=",signif(x$biasout, digits),"\n")
      
      cat("alternative hypothesis:", alt.hyp ,"\n")
      cat( ptrue, "\n")
    }
    
    if (x$power.reached == FALSE && x$bias.reached == FALSE && x$powerout == 0)
    {cat("Power can not be reached in the range of n","=",min(x$nit),",",max(x$nit),"\n")
      cat("null hypothesis can not be rejected for any number of groups in this range\n")
      
      cat("alternative hypothesis:", alt.hyp ,"\n")
      cat( ptrue, "\n")
    }
    
  } else if (names(x)[1] == "sout") {
    
    if (x$power.reached == TRUE && x$bias.reached == FALSE)
    {cat("Power was reached without violating bias restriction","\n")
      cat("for s","=",x$sout,"with power","=",signif(x$powerout, digits),"\n")
      cat("and bias","=",signif(x$biasout, digits),"\n")
      
      cat("alternative hypothesis:", alt.hyp ,"\n")
      cat( ptrue, "\n")
    }
    
    if (x$power.reached == FALSE && x$bias.reached == FALSE && x$powerout != 0)
    {cat("Power was not reached in the range of s","=",min(x$sit),",",max(x$sit),"\n")
      cat("Maximal power was reached for s","=",x$sout,"with power","=",signif(x$powerout, digits),"\n")
      cat("and bias","=",x$biasout,"\n")
      
      cat("alternative hypothesis:", alt.hyp ,"\n")
      cat( ptrue, "\n")
    }
    
    if (x$power.reached == FALSE && x$bias.reached == TRUE)
    {cat("Power can not be reached without violating bias restriction","\n")
      cat("Maximal power without violating biasrest","=", x$biasrest,"\n")
      cat("was reached for s","=",x$sout,"\n")
      cat("with power","=",signif(x$powerout, digits),"and bias","=",signif(x$biasout, digits),"\n")
      
      cat("alternative hypothesis:", alt.hyp ,"\n")
      cat( ptrue, "\n")
    }
    
    if (x$power.reached == FALSE && x$bias.reached == FALSE && x$powerout == 0)
    {cat("Power can not be reached in the range of s","=",min(x$sit),",",max(x$sit),"\n")
      cat("null hypothesis can not be rejected for any group size in this range\n")
      
      cat("alternative hypothesis:", alt.hyp ,"\n")
      cat( ptrue, "\n")
    }
  }
 
  invisible(x)
}




# Supporting function for estDesign()
##################################################################
# msep() function -                                              #
##################################################################

"msep" <- function(n,s,p.tr){
    
    expected <- 0
    for (y in 0:n)
    {
      expected <- expected + ((1 - (1 - y/n)^(1/s))*choose(n,y)*((1 - (1 - p.tr)^s)^y)*((1 - p.tr)^(s*(n - y))))
    }
    expected
    
    
    varsum <- 0
    for (y in 0:n)
    {
      varsum <- varsum + (((1 - y/n)^(2/s))*choose(n, n - y)*(((1 - p.tr)^s)^(n - y))*((1 - (1 - p.tr)^s)^y)) 
    }
    varp <- varsum - (1 - expected)^2
    
    expp <- expected
    bias <- expected - p.tr
    mse <- varp + bias^2
    
    
    list(varp = varp,
         mse = mse,
         bias = bias,
         exp = expected)
  }




##################################################################
# estDesign() function                                           #
##################################################################
# Brianna Hitt - 03.07.2020
# Changed "maximal" and "minimal" to "maximum" and "minimum", respectively
# 

"estDesign" <- function(n, smax, p.tr, biasrest = 0.05){
    
    if (length(n) != 1 || (n <= 1 | abs(round(n) - n) > 1e-07)) {stop("number of groups n must be specified as a single integer greater than 1")}
    if (length(p.tr) != 1 || p.tr > 1 || p.tr < 0) {stop("true proportion p.tr must be specified as a single number between 0 and 1")}
    if (length(smax) != 1 || (smax <= 1 | abs(round(smax) - smax) > 1e-07)) {stop("the maximal group size allowed in calculations must be a single integer greater than 1")}
    if (length(biasrest) != 1 || biasrest >= 1 || biasrest < 0) {stop("the maximally allowed bias(p) specified in biasrest must be a single number between 0 and 1, usually should be close to 0")}
    
  for (i in 2:smax) {
    temp <- msep(n = n, p.tr = p.tr, s = i)
    
    if (temp$bias > biasrest) {
      message("maximum group size within bias restriction is s = ", i - 1, "\n")
      res <- msep(n = n, p.tr = p.tr, s = i - 1)
      sout <- i - 1
      break()
    }
    
    if (i >= 2 && temp$mse > msep(n = n, p.tr = p.tr, s = i - 1)$mse) {
      message("minimum mse(p) is achieved with group size s = ", i - 1, "\n")
      res <- msep(n = n, p.tr = p.tr, s = i - 1)
      sout <- i - 1
      break()
    }
    
    if (i == smax && temp$mse <= msep(n = n, p.tr = p.tr, s = i - 1)$mse) {
      message("minimum mse(p) is achieved with group size s >= smax", "\n")
      res <- msep(n = n, p.tr = p.tr, s = i - 1)
      sout <- smax
      break()
    }
  }
  
  c(list("sout" = sout), res)
}




##################################################################
# designEst() function                                           #
##################################################################

#' @title Optimal group size determination based on minimal MSE 
#' when estimating an overall prevalence
#' 
#' @description Find the group size \kbd{s} for a fixed number of 
#' groups \kbd{n} and an assumed true proportion \kbd{p.tr}, for 
#' which the mean squared error (MSE) of the point estimator is 
#' minimal and bias is within a restriction.
#' 
#' @param n integer specifying the fixed number of groups.
#' @param smax integer specifying the maximum group size allowed in the 
#' planning of the design.
#' @param p.tr assumed true proportion of the "positive" trait 
#' in the population, specified as a value between 0 and 1.
#' @param biasrest a value between 0 and 1 specifying the absolute 
#' bias maximally allowed.
#' 
#' @details Swallow (1985) recommends the use of the upper bound of 
#' the expected range of the true proportion \kbd{p.tr} for optimization 
#' of the design. For further details, see Swallow (1985). Note that the 
#' specified number of groups must be less than \eqn{n=1020}.
#' 
#' @return A list containing:
#' \item{sout}{the group size \kbd{s} for which the MSE of the estimator is 
#' minimal for the given \kbd{n} and \kbd{p.tr} and for which the bias 
#' restriction \kbd{biasrest} is not violated. In the case that the minimum 
#' MSE is achieved for a group size \eqn{s>=smax}, the value of \kbd{smax} 
#' is returned.}
#' \item{varp}{the variance of the estimator.}
#' \item{mse}{the mean square error of the estimator.}
#' \item{bias}{the bias of the estimator.}
#' \item{exp}{the expected value of the estimator.}
#' 
#' @author This function was originally written by Frank Schaarschmidt 
#' as the \code{estDesign} function for the \code{binGroup} package.
#' 
#' @references 
#' \insertRef{Swallow1985}{binGroup2}
#' 
#' @seealso \code{\link{designPower}} for choice of the group testing 
#' design according to the power in a hypothesis test.
#' 
#' @family estimation functions
#' 
#' @examples 
#' # Compare to Table 1 in Swallow (1985):
#' designEst(n=10, smax=100, p.tr=0.001)
#' designEst(n=10, smax=100, p.tr=0.01)
#' designEst(n=25, smax=100, p.tr=0.05)
#' designEst(n=40, smax=100, p.tr=0.25)
#' designEst(n=200, smax=100, p.tr=0.30)

designEst <- estDesign

#