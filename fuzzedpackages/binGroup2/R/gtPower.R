##################################################################
# bgtPower and bgtPowerI() functions                             #
##################################################################

# Brianna Hitt - 02.23.2020
# Changed the capitalization of confidence interval methods to 
#   match those in the propCI() function

"bgtPower" <- 
  function(n, s, delta, p.hyp, conf.level = 0.95, 
           method = "CP", alternative="two.sided") {
    
    if ( any(n<=3) )
    {stop("the number of groups n allowed in calculations must be integers greater than 1")}
    
    if ( any(s<1) ) {stop("group size s must be specified as integers > 0")}
    
    if ( length(conf.level)!=1 || conf.level<0 || conf.level>1)
    {stop("conf.level must be a positive number between 0 and 1")}
    
    if ( length(p.hyp)!=1 || p.hyp>1 || p.hyp<0)
    {stop("true proportion p.hyp must be specified as a single number between 0 and 1")}
    
    method<-match.arg(method, choices=c("CP","Blaker","AC","score","Wald","soc"))
    
    alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))
    
    if (alternative=="less")
    {if ( any( p.hyp-delta < 0) || any(p.hyp-delta > 1) )
    {stop("alternative=less: specify delta as a number between 0 and the threshold p.hyp")}
    }
    
    if (alternative=="greater")
    {if ( any( p.hyp+delta < 0) || any(p.hyp+delta > 1) )
    {stop("alternative=greater: specify delta as a number between the threshold p.hyp and 1")}
    }
    
    if (alternative=="two.sided")
    {if ( any(p.hyp+delta < 0) || any(p.hyp+delta > 1) || any(p.hyp-delta < 0) || any(p.hyp-delta > 1))
    {stop("alternative=two.sided: specify delta as a number between the threshold p.hyp and 1")} 
    }
    
    # calculations:
    
    matnsp <- cbind(n,s,delta)
    matnsp <- cbind("ns"=matnsp[,1]*matnsp[,2], matnsp)
    power <- numeric(length=nrow(matnsp))
    bias <- numeric(length=nrow(matnsp))
    
    
    for( i in 1:length(power))
    {
      temp <- bgtPowerI(n=matnsp[i,2], s=matnsp[i,3], delta=matnsp[i,4], 
                        p.hyp=p.hyp, conf.level=conf.level, method=method, 
                        alternative=alternative)
      power[i] <- temp$power
      bias[i] <- temp$bias
      
    }
    return(cbind(matnsp, power, bias))
    
  }




# Brianna Hitt - 02-13-2020
# Changed class from "bgtPower" to "gtPower"
# Brianna Hitt - 02.23.2020
# Changed the capitalization to match

"bgtPowerI" <-
  function(n, s, delta, p.hyp, conf.level, method, alternative) {
    
    
    # 1) P.Ind returns TRUE in case that the CI doesnot contain p.hyp for a certain event Y=y
    
    
    P.Ind <- function(n,y,s,p.hyp,conf.level,method,alternative) {

      if (method == "score") {
        KI.Wilson <- bgtWilson(n = n,y = y,s = s,conf.level = conf.level, alternative = alternative) 
        
        (KI.Wilson[[1]] >= p.hyp || KI.Wilson[[2]] <= p.hyp)
      }
      
      else{if (method == "AC") {
        KI.AC <- bgtAC(n = n,y = y,s = s,conf.level = conf.level, alternative = alternative) 
        
        (KI.AC[[1]] >= p.hyp || KI.AC[[2]] <= p.hyp)
      }
        
        else{if (method == "Wald") {
          KI.Wald <- bgtWald(n = n,y = y,s = s,conf.level = conf.level, alternative = alternative)
          
          (KI.Wald[[1]] >= p.hyp || KI.Wald[[2]] <= p.hyp)
        }
          
          else{if (method == "CP") {
            KI.CP <- bgtCP(n = n,y = y,s = s,conf.level = conf.level, alternative = alternative) 
            
            (KI.CP[[1]] >= p.hyp || KI.CP[[2]] <= p.hyp)
          }
            
            else{if (method == "soc") {
              KI.SOC <- bgtSOC(n = n,y = y,s = s,conf.level = conf.level, alternative = alternative)
              
              (KI.SOC[[1]] >= p.hyp || KI.SOC[[2]] <= p.hyp)
            }
              
              else{if (method == "Blaker") {
                KI.Bl <- bgtBlaker(n = n,y = y,s = s,conf.level = conf.level) 
                
                if (alternative == "two.sided")
                {dec <- (KI.Bl[[1]] >= p.hyp || KI.Bl[[2]] <= p.hyp)}
                
                if (alternative == "less")
                {dec <- (KI.Bl[[2]] <= p.hyp)}
                
                if (alternative == "greater")
                {dec <- (KI.Bl[[1]] >= p.hyp)}
                dec
              }
                
                else{stop("argument method mis-specified")}}}}}}
    }
    
    # end of P.Ind
    
    # # #
    
    # 2) Probability of a certain event Y=y:
    
    bgt.prob <- function(n, y, s, p.tr)
    {
      theta <- 1 - (1 - p.tr)^s
      dbinom(x = y,size = n, prob = theta)
    }
    
    # # # 
    #  3) sum( P.Ind(y=y)*Prob(y=y) ) for all realizations of y:
    
    if (alternative == "less" || alternative == "greater")
    {
      
      if (alternative == "less") {p.tr = p.hyp - delta}
      if (alternative == "greater") {p.tr = p.hyp + delta}
      
      yvec <- 0:n
      
      probvec <- numeric(length = length(yvec))
      powvec <- numeric(length = length(yvec))
      expvec <- numeric(length = length(yvec))
      
      for (i in 1:length(yvec))
      {
        probvec[i] <- bgt.prob(n = n,y = yvec[i],s = s,p.tr = p.tr)
        powvec[i] <- P.Ind(n = n,y = yvec[i],s = s,p.hyp = p.hyp,conf.level = conf.level,method = method, alternative  =  alternative)
        expvec[i] <- (1 - (1 - yvec[i]/n)^(1/s))
      }
      powex <- sum(powvec * probvec) 
      expex <- sum(expvec * probvec)
      bias <- expex - p.tr
      
      out <- list(power = powex, bias = bias, p.tr = p.tr)
      
    }
    
    if (alternative == "two.sided")
    {
      p.trl = p.hyp - delta
      p.trg = p.hyp + delta
      
      yvec <- 0:n
      powvec <- numeric(length = length(yvec))
      expvec <- numeric(length = length(yvec))
      probvecl <- numeric(length = length(yvec))
      probvecg <- numeric(length = length(yvec))
      
      
      
      for (i in 1:length(yvec))
      {
        probvecl[i] <- bgt.prob(n = n,y = yvec[i],s = s,p.tr = p.trl)
        probvecg[i] <- bgt.prob(n = n,y = yvec[i],s = s,p.tr = p.trg)
        
        powvec[i] <- P.Ind(n = n,y = yvec[i],s = s,p.hyp = p.hyp,conf.level = conf.level,method = method, alternative  =  alternative)
        
        expvec[i] <- (1 - (1 - yvec[i]/n)^(1/s))
      }
      
      powexl <- sum(powvec * probvecl) 
      expexl <- sum(expvec * probvecl)
      biasl <- expexl - p.trl
      
      powexg <- sum(powvec * probvecg) 
      expexg <- sum(expvec * probvecg)
      biasg <- expexg - p.trg
      
      out <- list(power = min(powexl,powexg),bias = max(biasl,biasg), p.tr = c(p.trl,p.trg))
      
    }
    
    class(out) <- "gtPower"
    out
    
  }




##################################################################
# gtPower() function                                             #
##################################################################
# This was the "bgtPower" function from Frank Schaarschmidt

#' @title Power to reject a hypothesis for one proportion in group testing
#' 
#' @description This function calculates the power to reject a hypothesis 
#' in a group testing experiment, using confidence intervals for the 
#' decision. This function also calculates the bias of the point estimator 
#' for a given \eqn{n}, \eqn{s}, and true, unknown proportion.
#' 
#' @param n integer specifying the number of groups. A vector of integers is 
#' also allowed.
#' @param s integer specifying the common group size. A vector of integers is 
#' also allowed.
#' @param delta the absolute difference between the true proportion and the 
#' hypothesized proportion. A vector is also allowed.
#' @param p.hyp the proportion in the hypotheses, specified as a value between 
#' 0 and 1.
#' @param conf.level confidence level required for the decision on the 
#' hypotheses.
#' @param method character string specifying the confidence interval method 
#' (see \code{\link{propCI}}) to be used.
#' @param alternative character string defining the alternative hypothesis, 
#' either \kbd{"two.sided"}, \kbd{"less"}, or \kbd{"greater"}.
#' 
#' @details The power of a hypothesis test performed by a confidence 
#' interval is defined as the probability that a confidence interval 
#' excludes the threshold parameter (\kbd{p.hyp}) of the null hypothesis, 
#' as described in Schaarschmidt (2007). Due to discreteness, the power 
#' does not increase monotonely for an increasing number of groups \eqn{n} 
#' or group size \eqn{s}, but exhibits local maxima and minima, depending 
#' on \eqn{n}, \eqn{s}, \kbd{p.hyp}, and \kbd{conf.level}.
#' 
#' Additional to the power, the bias of the point estimator is calculated 
#' according to Swallow (1985). If vectors are specified for \eqn{n}, 
#' \eqn{s}, and (or) delta, a matrix will be constructed and power and 
#' bias are calculated for each line in this matrix.
#' 
#' @return A matrix containing the following columns:
#' \item{ns}{a vector of the total sample size, \eqn{n*s}.}
#' \item{n}{a vector of the number of groups.}
#' \item{s}{a vector of the group sizes.}
#' \item{delta}{a vector of the delta values.}
#' \item{power}{the power to reject the given null hypothesis.}
#' \item{bias}{the bias of the estimator for the specified 
#' \eqn{n}, \eqn{s}, and the true proportion.}
#' 
#' @author This function was originally written as \code{bgtPower} by Frank 
#' Schaarschmidt for the \code{binGroup} package. Minor modifications have 
#' been made for inclusion of the function in the \code{binGroup2} package.
#' 
#' @references 
#' \insertRef{Schaarschmidt2007}{binGroup2}
#' 
#' \insertRef{Swallow1985}{binGroup2}
#' 
#' @seealso \code{\link{propCI}} for confidence intervals and 
#' \code{\link{gtTest}} for hypothesis tests for one proportion from a 
#' group testing experiment.
#' 
#' @family estimation functions
#' 
#' @examples 
#' # Calculate the power for the design
#' #   in the example given in Tebbs and Bilder(2004):
#' #   n=24 groups each containing 7 insects
#' #   if the true proportion of virus vectors
#' #   in the population is 0.04 (4 percent),
#' #   the power to reject H0: p>=0.1 using an
#' #   upper Clopper-Pearson ("CP") confidence interval
#' #   is calculated with the following call:
#' gtPower(n=24, s=7, delta=0.06, p.hyp=0.1,
#'         conf.level=0.95, alternative="less", method="CP")
#'
#' # Explore development of power and bias for varying 
#' #   n, s, delta. How much can we decrease the number of 
#' #   groups (costly tests to be performed) by pooling the same 
#' #   number of 320 individuals to groups of increasing size 
#' #   without largely decreasing power?
#' gtPower(n=c(320,160,80,64,40,32,20,10,5),
#'         s=c(1,2,4,5,8,10,16,32,64), delta=0.01, p.hyp=0.02)
#'                   
#' # What happens to the power for increasing differences
#' #   between the true proportion and the threshold proportion?
#' gtPower(n=50, s=10, delta=seq(from=0, to=0.01, by=0.001),
#'         p.hyp=0.01, method="CP")
#'          
#' # Calculate power with a group size of 1 (individual testing).
#' gtPower(n=100, s=1, delta=seq(from=0, to=0.01, by=0.001),
#'         p.hyp=0.01, method="CP")

gtPower <- bgtPower


#