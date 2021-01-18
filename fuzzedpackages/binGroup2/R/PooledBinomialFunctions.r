##################################################################
# propDiffCI() function                                          #
##################################################################
# Brianna Hitt - 02.10.2020
# Minor changes were made to the capitalization of confidence interval 
#   method names - only the first letter is capitalized for methods 
#   named after people; otherwise, all lowercase
# The class of the returned object was also changed from 
#   "poolbindiff" to "propDiffCI"

# Brianna Hitt - 03.06.2020
# Removed MIR method from the function 
# Removed the "scale" argument from the function

#' @title Confidence intervals for the difference of proportions
#' in group testing
#'
#' @description Calculates confidence intervals for the difference of two 
#' proportions based on group testing data.
#'
#' @param x1 vector specifying the observed number of positive
#' groups among the number of groups tested (\kbd{n1}) in population 1.
#' @param m1 vector of corresponding group sizes in population 1.
#' Must have the same length as \kbd{x1}.
#' @param x2 vector specifying the observed number of positive
#' groups among the number of groups tested (\kbd{n2}) in population 2.
#' @param m2 vector of corresponding group sizes in population 2.
#' Must have the same length as \kbd{x2}.
#' @param n1 vector of the corresponding number of groups with
#' sizes \kbd{m1}.
#' @param n2 vector of the corresponding number of groups with
#' sizes \kbd{m2}.
#' @param pt.method character string specifying the point estimator to 
#' compute. Options include \kbd{"Firth"} for the bias-preventative estimator 
#' (Hepworth & Biggerstaff, 2017), the default \kbd{"Gart"} for the 
#' bias-corrected MLE (Biggerstaff, 2008), \kbd{"bc-mle"} (same as 
#' \kbd{"Gart"} for backward compatibility), and \kbd{"mle"} for the MLE.
#' @param ci.method character string specifying the confidence interval to 
#' compute. Options include \kbd{"skew-score"} for the skewness-corrected, 
#' \kbd{"score"} for the score (the default), \kbd{"bc-skew-score"} for the 
#' bias- and skewness-corrected, \kbd{"lrt"} for the likelihood ratio test, 
#' and \kbd{"Wald"} for the Wald interval. See Biggerstaff (2008) for
#' additional details.
#' @param conf.level confidence level of the interval.
#' @param tol the accuracy required for iterations in internal functions.
#'
#' @details Confidence interval methods include the Wilson score 
#' (\kbd{ci.method="score"}), skewness-corrected score 
#' (\kbd{ci.method="skew-score"}), bias- and skewness-corrected score 
#' (\kbd{ci.method="bc-skew-score"}), likelihood ratio test 
#' (\kbd{ci.method="lrt"}), and Wald (\kbd{ci.method="Wald"}) interval. 
#' For computational details, simulation results, and recommendations on 
#' confidence interval methods, see Biggerstaff (2008).
#' 
#' Point estimates available include the MLE (\kbd{pt.method="mle"}), 
#' bias-corrected MLE (\kbd{pt.method="Gart"} or \kbd{pt.method="bc-mle"}), 
#' and bias-preventative (\kbd{pt.method="Firth"}). For additional details 
#' and recommendations on point estimation, see Hepworth and Biggerstaff 
#' (2017).
#'
#' @return A list containing:
#' \item{d}{the estimated difference of proportions.}
#' \item{lcl}{the lower confidence limit.}
#' \item{ucl}{the upper confidence limit.}
#' \item{pt.method}{the method used for point estimation.}
#' \item{ci.method}{the method used for confidence interval estimation.}
#' \item{conf.level}{the confidence level of the interval.}
#' \item{x1}{the numbers of positive groups in population 1.}
#' \item{m1}{the sizes of the groups in population 1.}
#' \item{n1}{the numbers of groups with corresponding group sizes
#' \kbd{m1} in population 1.}
#' \item{x2}{the numbers of positive groups in population 2.}
#' \item{m2}{the sizes of the groups in population 2.}
#' \item{n2}{the numbers of groups with corresponding group sizes
#' \kbd{m2} in population 2.}
#'
#' @author This function was originally written as the \code{pooledBinDiff} 
#' function by Brad Biggerstaff for the \code{binGroup} package. Minor 
#' modifications were made for inclusion of the function in the 
#' \code{binGroup2} package.
#' 
#' @references 
#' \insertRef{Biggerstaff2008}{binGroup2}
#' 
#' \insertRef{Hepworth2017}{binGroup2}
#' 
#' @seealso \code{\link{propCI}} for confidence intervals for one proportion 
#' in group testing, \code{\link{gtTest}} for hypothesis tests in group 
#' testing, and \code{\link{gtPower}} for power calculations in group testing.
#' 
#' @family estimation functions
#' 
#' @examples
#' # Estimate the prevalence in two populations 
#' #   with multiple groups of various sizes:
#' # Population 1:
#' #   0 out of 5 groups test positively with 
#' #   groups of size 1 (individual testing); 
#' #   0 out of 5 groups test positively with 
#' #   groups of size 5;
#' #   1 out of 5 groups test positively with 
#' #   groups of size 10; and
#' #   2 out of 5 groups test positively with 
#' #   groups of size 50.
#' # Population 2:
#' #   0 out of 5 groups test positively with
#' #   groups of size 1 (individual testing);
#' #   1 out of 5 groups test positively with 
#' #   groups of size 5;
#' #   0 out of 5 groups test positively with 
#' #   groups of size 10; and 
#' #   4 out of 5 groups test positively with 
#' #   groups of size 50.
#' x1 <- c(0,0,1,2)
#' m <- c(1,5,10,50)
#' n <- c(5,5,5,5)
#' x2 <- c(0,1,0,4)
#' propDiffCI(x1=x1, m1=m, x2=x2, m2=m, n1=n, n2=n, 
#'            pt.method="Gart", ci.method="score")
#'
#' # Compare recommended methods:
#' propDiffCI(x1=x1, m1=m, x2=x2, m2=m, n1=n, n2=n,
#'            pt.method="mle", ci.method="lrt")
#'
#' propDiffCI(x1=x1, m1=m, x2=x2, m2=m, n1=n, n2=n,
#'            pt.method="mle", ci.method="score")
#'
#' propDiffCI(x1=x1, m1=m, x2=x2, m2=m, n1=n, n2=n,
#'            pt.method="mle", ci.method="skew-score")

# Brianna Hitt - 04.02.2020
# Changed cat()/print() to warning()

propDiffCI <- 
  function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),
           pt.method = c("Firth", "Gart", "bc-mle", "mle"),
           ci.method = c("skew-score", "bc-skew-score", "score", "lrt", "Wald"),
           conf.level = 0.95, tol = .Machine$double.eps^0.5)
  {
    
    scale <- 1
    alpha <- 1 - conf.level
    
    call <- match.call()
    pt.method <- match.arg(pt.method)
    ci.method <- match.arg(ci.method)
    if (ci.method == "mir" | pt.method == "mir") {
      ci.method <- "mir"
      pt.methid <- "mir"
    }
    # Brianna Hitt - 02-15-2020
    # Changed the capitalization for "Firth" and "Gart"
    if (pt.method == "Gart") pt.method <- "bc-mle" # backward compatability
    
    switch(pt.method,
           "Firth" = {d <- pooledbinom.firth(x1,m1,n1,tol) - pooledbinom.firth(x2,m2,n2,tol)},
           "mle" = {d <- pooledbinom.mle(x1,m1,n1,tol) - pooledbinom.mle(x2,m2,n2,tol)},
           "bc-mle" = {d <- pooledbinom.cmle(x1,m1,n1,tol) - pooledbinom.cmle(x2,m2,n2,tol)},
           "mir" = {d <- pooledbinom.mir(x1,m1,n1) - pooledbinom.mir(x2,m2,n2)}
    )
    if ((pooledbinom.cmle(x1,m1,n1,tol) < 0 | pooledbinom.cmle(x2,m2,n2,tol) < 0 ) & pt.method == "bc-mle") {
      pt.method <- "mle"
      warning("Bias-correction results in a negative point estimate; using MLE\n")
      d <- pooledbinom.mle(x1,m1,n1,tol) - pooledbinom.mle(x2,m2,n2,tol)
    }
    
    switch(ci.method,
           "skew-score" = {ci.d <- pooledbinom.diff.cscore.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "bc-skew-score" = {ci.d <- pooledbinom.diff.bcscore.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "score" = {ci.d <- pooledbinom.diff.score.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "lrt" = {ci.d <- pooledbinom.diff.lrt.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "Wald" = {ci.d <- pooledbinom.diff.wald.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]},
           "mir" = {ci.d <- pooledbinom.diff.mir.ci(x1,m1,x2,m2,n1,n2,alpha)[4:5]}
    )
    # Brianna Hitt - 02-15-2020
    # Changed the class from "poolbindiff" to "propDiffCI"
    # Brianna Hitt - 03-06-2020
    # Removed the "scale" object from the returned list 
    #   and added conf.level in place of alpha
    # structure(list(d = d, lcl = ci.d[1], ucl = ci.d[2], pt.method = pt.method,
    #                ci.method = ci.method, alpha = alpha, scale = scale, 
    #                x1 = x1, m1 = m1, n1 = n1, x2 = x2, m2 = m2, n2 = n2, call = call), 
    #           class = "propDiffCI")
    
    structure(list(d = d, lcl = ci.d[1], ucl = ci.d[2], pt.method = pt.method,
                   ci.method = ci.method, conf.level = conf.level,
                   x1 = x1, m1 = m1, n1 = n1, x2 = x2, m2 = m2, n2 = n2, call = call), 
              class = "propDiffCI")
  }




##################################################################
# print.propDiffCI() function                                    #
##################################################################

# This function was originally written as print.poolbindiff() for binGroup
# Brianna Hitt - 02/10/2020 
# Minor modifications were made to make this output more consistent 
#   with the output provided by print.propCI

#' @title Print method for objects of class "propDiffCI"
#' 
#' @description Print method for objects of class "propDiffCI" 
#' created by the \code{\link{propDiffCI}} function.
#' 
#' @param x An object of class "propDiffCI" (\code{\link{propDiffCI}}).
#' @param ... Additional arguments to be passed to \code{print}.
#' 
#' @return A print out of the point estimate and confidence interval 
#' found with \code{\link{propDiffCI}}.
#' 
#' @author This function was originally written as \code{print.poolbindiff}
#' by Brad Biggerstaff for the \code{binGroup} package. Minor 
#' modifications were made for inclusion of the function in the 
#' \code{binGroup2} package.

# Brianna Hitt - 03.07.2020
# Removed the "scale" argument from the function

"print.propDiffCI" <- function(x, ...){
  scale <- 1
  
  args <- list(...)
  if (is.null(args$digits)) {digits <- 4}
  else{digits <- args$digits}
  
  switch(x$pt.method, 
         "Firth" = PtEstName <- "Firth's Correction", 
         "Gart" = PtEstName <- "Gart's Correction", 
         "bc-mle" = PtEstName <- "Gart's Correction", 
         "mle" = PtEstName <- "Maximum Likelihood", 
         "mir" = PtEstName <- "Minimum Infection Rate")
  
  switch(x$ci.method, 
         "score" = CIEstName <- "Wilson score", 
         "Wald" = CIEstName <- "Wald", 
         "skew-score" = CIEstName <- "Skew-Corrected score (Gart)", 
         "bc-skew-score" = CIEstName <- "Bias- & Skew-Corrected score (Gart)", 
         "lrt" = CIEstName <- "Likelihood Ratio Test Inversion", 
         "mir" = CIEstName <- "Minimum Infection Rate")
  
  cat("\n")
  cat(x$conf.level*100, "percent", CIEstName,
      "confidence interval:\n", sep = " ")
  # cat(x$conf.level*100, "percent confidence interval:\n", sep = " ")
  cat(" [ ", signif(scale*x$lcl, digits), ", ", 
      signif(scale*x$ucl, digits), " ]\n", sep = "")
  cat("Point estimate (", PtEstName, "): ", signif(x$d, digits),
      "\n", sep = "")
  # cat("Point estimate:", signif(x$estimate, digits), "\n", sep = " ")
  # cat("Scale:", scale, sep = " ")
  
  invisible(x)
}



################################################################################
# One-sample functions
################################################################################
"pooledbinom.loglike" <-
function(p,x,m,n=rep(1,length(m)))
{
	if(p > 0)
		sum(x*log((1-(1-p)^m))) + log(1-p)*sum(m*(n-x))
	else 0
}

"pooledbinom.loglike.vec" <-
function(p,x,m,n=rep(1,length(m)))
{
	np <- length(p)
	lik <- vector(length=np)
	for(i in 1:np) lik[i] <- pooledbinom.loglike(p[i],x,m,n)
	lik
}

"score.p" <-
function(p, x, m, n = rep(1,length(m)))
{
	if(sum(x) == 0 & p >= 0 & p < 1) return(-sum(m*n)/(1-p))
	if(p < 0 | p >= 1) return(0)
	sum(m*x/(1-(1-p)^m) - m*n)/(1-p)
}

"pooledbinom.mir" <-
function(x,m,n=rep(1,length(x)))
{
	sum(x)/sum(m*n)
}

"pooledbinom.mle" <-
function(x, m, n = rep(1., length(x)), tol = 1e-008)
{
#
# This is the implementation using Newton-Raphson, as given
# in the Walter, Hildreth, Beaty paper, Am. J. Epi., 1980
#
	if(length(m) == 1.) m <- rep(m, length(x)) else if(length(m) != length(x))
		stop("\n ... x and m must have same length if length(m) > 1")
	if(any(x > n))
		stop("x elements must be <= n elements")
	if(all(x == 0.))
		return(0.)
	if(sum(x) == sum(n)) return(1)
	#p.new <- 1 - (1 - sum(x)/sum(n))^(1/mean(m)) # starting value
        # Changed by FS 01.03.2018 acc. to email of BB:
        p.new <- sum(x)/sum(m*n)
	done <- 0
	N <- sum(n * m)
	while(!done) {
		p.old <- p.new
		p.new <- p.old - (N - sum((m * x)/(1 - (1 - p.old)^m)))/
					sum((m^2 * x * (1 - p.old)^(m - 1))/(1 - (1 -
					p.old)^m)^2)
		if(abs(p.new - p.old) < tol)
			done <- 1
	}
	p.new
}


"pooledbinom.bias" <-
function(p,m,n=rep(1,length(m)))
{
	if(p > 0)
		sum((m-1)*m^2*n*(1-p)^(m-3)/(1-(1-p)^m)) *
				pooledbinom.mle.var(p,m,n)^2 / 2
	else 0
}

# c is bias-Corrected
"pooledbinom.cmle" <-
function(x, m, n = rep(1, length(m)), tol= 1e-8)
{
	phat <- pooledbinom.mle(x,m,n)
	bias <- pooledbinom.bias(phat,m,n)
	phat - bias
}

"pooledbinom.mle.var" <-
function(p, m, n = rep(1, length(m)))
{
	if(p > 0 & p < 1)
	1/sum((m^2 * n * (1 - p)^(m - 2))/(1 - (1 - p)^m))
	else 0 #1/sum(m*n)
}

"pooledbinom.I" <-
function(p,m,n=rep(1,length(m)))
{
	if(p > 0)
	sum((m^2 * n * (1 - p)^(m - 2))/(1 - (1 - p)^m))
	else 0 # Not sure whether to set this to 0 or not
}

"pooledbinom.mu3" <-
function(p, m, n = rep(1,length(m)))
{
	if(p > 0)
		sum((m^3 * n * (1-p)^(m-3) * (2 * (1-p)^m - 1))/(1-(1-p)^m)^2)
	else 0
}
#
#
# "pooledbinom.firth" <-
#   function(x, m, n = rep(1., length(x)), tol = 1e-008, rel.tol = FALSE)
#   {
#     #
#     # This is the implementation using Newton-Raphson
#     #
#
#     if(length(m) == 1.) m <- rep(m, length(x)) else if(length(m) != length(x))
#       stop("\n ... x and m must have same length if length(m) > 1")
#     if(any(x > n))
#       stop("x elements must be <= n elements")
#     if(all(x == 0.)) return(0.)
#     #if(sum(x) == sum(n)) return(1.)
#
#     # assign a convergence criterion function to avoid repeated
#     # checks of rel.tol during iteration
#
#     if(rel.tol) "converge.f" <- function(old,new) abs((old-new)/old)
#     else "converge.f" <- function(old,new) abs(old-new)
#     N <- sum(n * m)
#
#     # use proportion of positive groups, sum(x)/sum(n)
#     # and inverse of average pool size, sum(n)/N
#     # ...but this doesn't work when all groups are positive, while the
#     # ...MIR does, so use MIR in that case
#
#     if(sum(x) == sum(n)){
#       p.new <- sum(x)/N
#     } else {
#       p.new <- 1 - (1 - sum(x)/sum(n))^(sum(n)/N)
#     }
#     done <- FALSE
#     while(!done) {
#       p.old <- p.new
#       vi.p <- m^2 * n * (1-p.old)^(m-2) / (1 - (1-p.old)^m)
#       wi.p <- vi.p / sum(vi.p)
#       vi.p.prime <- m^2 * n * (1-p.old)^(m-3) * (2 * (1-(1-p.old)^m) - m) / (1-(1-p.old)^m)^2
#       wi.p.prime <- (vi.p.prime * sum(vi.p) - vi.p * sum(vi.p.prime)) / sum(vi.p)^2
#       p.new <- p.old +  (sum(m*x/(1-(1-p.old)^m) - m*wi.p/2) - (N-1/2)) /
#         sum(m^2*x*(1-p.old)^(m-1) / (1-(1-p.old)^m)^2 + m * wi.p.prime / 2)
#       # put in a catch for some wacky cases
#       if(p.new < 0)
#         p.new <- (1 - (1 - sum(x)/sum(n))^(sum(n)/N)) / 2 #sum(x)/N  # 0.5 * MIR
#       if(converge.f(p.old, p.new) < tol) done <- TRUE
#     }
#     p.new
#   }


"pooledbinom.firth" <-
function(x, m, n = rep(1., length(x)), tol = 1e-008, rel.tol = FALSE)
{
   #
   # This is the implementation using Newton-Raphson
   #
   if(length(m) == 1.) m <- rep(m, length(x)) else if(length(m) != length(x))
      stop("\n ... x and m must have same length if length(m) > 1")
   if(any(x > n))
      stop("x elements must be <= n elements")
   if(all(x == 0.)) return(0.)
   #if(sum(x) == sum(n)) return(1.)

   # assign a convergence criterion function to avoid repeated
   # checks of rel.tol during iteration
   if(rel.tol) "converge.f" <- function(old,new) abs((old-new)/old)
   else "converge.f" <- function(old,new) abs(old-new)

   N <- sum(n * m)
   # use proportion of positive groups, sum(x)/sum(n)
   # and inverse of average pool size, sum(n)/N
   # ...but this doesn't work when all groups are positive, while the
   # ...MIR does, so use MIR in that case
   if(sum(x) == sum(n)){
      p.new <- sum(x)/N
   } else {
      p.new <- 1 - (1 - sum(x)/sum(n))^(sum(n)/N)
   }

   done <- FALSE
   while(!done) {
      p.old <- p.new
      vi.p <- m^2 * n * (1-p.old)^(m-2) / (1 - (1-p.old)^m)
      wi.p <- vi.p / sum(vi.p)
      vi.p.prime <- m^2 * n * (1-p.old)^(m-3) * (2 * (1-(1-p.old)^m) - m) / (1-(1-p.old)^m)^2
      wi.p.prime <- (vi.p.prime * sum(vi.p) - vi.p * sum(vi.p.prime)) / sum(vi.p)^2
      p.new <- p.old +  (sum(m*x/(1-(1-p.old)^m) - m*wi.p/2) - (N-1/2)) /
         sum(m^2*x*(1-p.old)^(m-1) / (1-(1-p.old)^m)^2 + m * wi.p.prime / 2)

      if(converge.f(p.old, p.new) < tol) done <- TRUE
   }
   p.new
}

"pooledbinom.score.ci" <-
function(x, m, n = rep(1., length(x)), tol = .Machine$double.eps^0.75, alpha =
	0.05)
{
	f <- function(p0, x, mm, n, alpha)
	{
		# NOTE: I use the square of the score statistic and chisq
		score.p(p0, x, mm, n)^2 * pooledbinom.mle.var(p0, mm, n) -
			qchisq(1 - alpha, 1)
	}
	# The MLE is just used for a sensible cut-value for the search ...
	# it's not needed in the computation
	p.hat <- pooledbinom.mle(x, m, n, tol)
	if(sum(x) == 0) {
		lower.limit <- 0
	}
	else {
		root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd =
			0., upper = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n,
			alpha = alpha)
		lower.limit <- uniroot(f, lower = root.brak[1], upper =
			root.brak[2], tol = .Machine$double.eps^0.75, x = x,
			mm = m, n = n, alpha = alpha)$root
	}
	if(sum(x) == sum(n)) {
		upper.limit <- 1
	}
	else {
		root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10,
			lbnd = p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x,
			mm = m, n = n, alpha = alpha)
		upper.limit <- uniroot(f, lower = root.brak[1], upper =
			root.brak[2], tol = .Machine$double.eps^0.75, x = x,
			mm = m, n = n, alpha = alpha)$root
	}
	c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha
		)
}


"pooledbinom.cscore.ci"<-
function(x, m, n = rep(1., length(x)), tol = 1e-008, alpha = 0.05)
{
	# use gamm not gam b/c of the function gam()
	f <- function(p, x, mm, n, alpha)
	{
		gamm <- function(p0, mm, n = rep(1, length(mm)))
		{
			pooledbinom.mu3(p0, mm, n) * pooledbinom.mle.var(p0,
				mm, n)^(3/2)
		}
		# NOTE: I use the square of the score statistic and chisq
		(score.p(p, x, mm, n) * sqrt(pooledbinom.mle.var(p, mm, n)) -
			(gamm(p, mm, n) * (qchisq(1 - alpha, 1) - 1))/6)^2 -
			qchisq(1. - alpha, 1)
	}
	# The MLE is just used for a sensible cut-value for the search ...
	# it's not needed in the computation
	p.hat <- pooledbinom.mle(x, m, n, tol) # used to be cmle
	if(sum(x) == 0)
		return(pooledbinom.score.ci(x, m, n, tol, alpha))
	# Effectively "else"
	root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd = 0., upper
		 = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n, alpha = alpha)
	lower.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
		2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
		 = alpha)$root
	# use mm because m is confused with maxiter
	root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10, lbnd =
		p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x, mm = m, n = n,
		alpha = alpha)
	upper.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
		2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
		 = alpha)$root
	c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha
		)
}


"pooledbinom.bcscore.ci"<-
function(x, m, n = rep(1., length(x)), tol = 1e-008, alpha = 0.05)
{
	f <- function(p, x, mm, n, alpha)
	{
		gamm <- function(p0, mm, n = rep(1, length(mm)))
		{
			pooledbinom.mu3(p0, mm, n) * pooledbinom.mle.var(p0,
				mm, n)^(3/2)
		}
		# NOTE: I use the square of the score statistic and chisq
		(score.p(p, x, mm, n) * sqrt(pooledbinom.mle.var(p, mm, n)) -
			pooledbinom.bias(p, mm, n) - (gamm(p, mm, n) * (qchisq(
			1 - alpha, 1) - 1))/6)^2 - qchisq(1. - alpha, 1)
	}
	# The MLE is just used for a sensible cut-value for the search ...
	# it's not needed in the computation
	p.hat <- pooledbinom.cmle(x, m, n, tol)
	p.hat <- pooledbinom.cmle(x, m, n, tol)
	if(sum(x) == 0)
		return(pooledbinom.score.ci(x, m, n, tol, alpha))
	# Effectively "else"
	root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd = 0., upper
		 = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n, alpha = alpha)
	lower.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
		2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
		 = alpha)$root
	# use mm because m is confused with maxiter
	root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10, lbnd =
		p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x, mm = m, n = n,
		alpha = alpha)
	upper.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
		2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, alpha
		 = alpha)$root
	c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha
		)
}


"pooledbinom.lrt.ci"<-
function(x, m, n = rep(1., length(x)), tol = 1e-008, alpha = 0.05)
{
	p.hat <- pooledbinom.mle(x, m, n, tol)
	f <- function(p, x, mm, n, f.tol, alpha, p.hat)
	{
		if(p.hat == 0.)
			return( - Inf)
		else -2. * (pooledbinom.loglike(p, x, mm, n) -
				pooledbinom.loglike(p.hat, x, mm, n)) - qchisq(
				1. - alpha, 1.)
	}
	if(sum(x) == 0)
		return(pooledbinom.score.ci(x, m, n, tol, alpha))
	# Effectively "else"
	root.brak <- bracket.bounded.root(f, lower = p.hat/10, lbnd = 0., upper
		 = p.hat/2, ubnd = p.hat, x = x, mm = m, n = n, f.tol = tol,
		alpha = alpha, p.hat = p.hat)
	lower.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
		2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, f.tol
		 = tol, alpha = alpha, p.hat = p.hat)$root
	# use mm because m is confused with maxiter
	#print(lower.limit)
	root.brak <- bracket.bounded.root(f, lower = (p.hat + 1)/10, lbnd =
		p.hat, upper = (p.hat + 1)/2, ubnd = 1., x = x, mm = m, n = n,
		f.tol = tol, alpha = alpha, p.hat = p.hat)
	upper.limit <- uniroot(f, lower = root.brak[1], upper = root.brak[
		2], tol = .Machine$double.eps^0.75, x = x, mm = m, n = n, f.tol
		 = tol, alpha = alpha, p.hat = p.hat)$root
	c(p = p.hat, lower = lower.limit, upper = upper.limit, alpha = alpha)
}



"pooledbinom.wald.ci"<-
function(x, m, n = rep(1, length(x)), tol = 1e-008, alpha = 0.05)
{
	if(length(m) == 1)
		m <- rep(m, length(x))
	p.hat <- pooledbinom.mle(x, m, n, tol)
	p.stderr <- sqrt(pooledbinom.mle.var(p.hat, m, n))
	z <- qnorm(1 - alpha/2)
	c(p = p.hat, lower = max(0, p.hat - z * p.stderr), upper = min(1, p.hat +
		z * p.stderr), alpha = alpha)
}

"pooledbinom.mir.ci" <-
function(x, m, n = rep(1, length(x)), tol = 1e-008, alpha = 0.05)
{
	N <- sum(m*n)
	mir <- sum(x)/N
	mir.stderr <- sqrt(mir*(1-mir)/N)
	z <- qnorm(1-alpha/2)
	c(p=mir,lower=max(0, mir - z * mir.stderr), upper = min(1, mir + z * mir.stderr),alpha=alpha)
}

#################################################################################
# two-sample functions
#################################################################################
"pooledbinom.diff.loglike" <-
function(d,s,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
	pooledbinom.loglike((s+d)/2,x1,m1,n1) + pooledbinom.loglike((s-d)/2,x2,m2,n2)


"score.p" <-
function(p, x, m, n = rep(1,length(m)))
{
	if(sum(x) == 0 & p > 0) return(-sum(m*n)/(1-p))
	if(p < 0 | p >= 1) return(0)
	sum(m*x/(1-(1-p)^m) - m*n)/(1-p)
}


"score2" <-
function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
	0.05*(score.p((s+d)/2,x1,m1,n1) + score.p((s-d)/2,x2,m2,n2))
}

"score2.vec" <-
function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
	ns <- length(s)
	s2 <- vector(length=ns)
	for(i in 1:ns)
		s2[i] <- score2(s[i],d,x1,m1,x2,m2,n1,n2)
	s2
}

# VITAL here to have tol = something very small!!!!
"Shatd" <-
function(d, x1, m1, x2, m2, n1, n2)
{
	#if(sum(x1)==0 & sum(x2)==0){
	#	N1 <- sum(m1*n1)
	#	N2 <- sum(m2*n2)
	#	return(2 - abs(N1-N2)*abs(d)/(N1+N2))
	#}
	if(sum(x1) == 0){
		if(d < 0){
			# g is value of dScore(s,d)/ds on the edge s = -d
			# must be vectorized for uniroot()
			g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
					{
						nd <- length(d)
						gv <- vector(length=nd)
						for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
						gv
					}
			brkt <- bracket.bounded.root(g,lower=-0.5,lbnd=-1,upper=-0.001,ubnd=0,
								x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
			dstar <- uniroot(g,lower = brkt[1],upper = brkt[2], tol = .Machine$double.eps^0.75,
								x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
#			dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
#							tol = .Machine$double.eps^0.75,
#						x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
			if(d <= dstar){
				return(-d)
			} else {

			brkt <- bracket.bounded.root(score2,lower= -0.5,lbnd=-d,upper=0,ubnd=2+d,
								d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
				ans <- uniroot(score2,lower = brkt[1],upper = brkt[2],
						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root

#				ans <- uniroot(score2,lower = -d + .Machine$double.eps^0.5,upper = 2 + d - .Machine$double.eps^0.5,
#						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
				return(ans)
		}
		}
		else {

			brkt <- bracket.bounded.root(score2,lower= d+0.001,lbnd=d,upper=2-d-0.001,ubnd=2-d,
								d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
			ans <- uniroot(score2,lower= brkt[1],upper = brkt[2],
						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root

#			ans <- uniroot(score2,lower= d + .Machine$double.eps^0.5,upper = 2 - d -.Machine$double.eps^0.5,
#						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
		    return(ans)
		}
	}
	if(sum(x2) == 0){
		if(d > 0){
			# g is value of dScore(s,d)/ds on the edge s = d
			# must be vectorized for uniroot()
			g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
					{
						nd <- length(d)
						gv <- vector(length=nd)
						for(i in 1:nd) gv[i] <- score2(d[i],d[i],x1,m1,x2,m2,n1,n2)
						gv
					}
			brkt <- bracket.bounded.root(g,lower=0.001,lbnd=0,upper=0.5,ubnd=1,
								x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
				dstar <- uniroot(g,lower = brkt[1],upper = brkt[2],
						tol = .Machine$double.eps^0.75,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
#				dstar <- uniroot(g,lower = 0 + .Machine$double.eps^0.5,upper = 1 -.Machine$double.eps^0.5,
#						tol = .Machine$double.eps^0.75,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
			if(d >= dstar)
				return(d)
			else {

			brkt <- bracket.bounded.root(score2,lower=d+0.001,lbnd=d,upper=2-d-0.001,ubnd=2-d,
								d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
			ans <- uniroot(score2,lower=brkt[1],upper = brkt[2],
						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root

#				ans <- uniroot(score2,lower= d + .Machine$double.eps,upper = 2 - d - .Machine$double.eps,
#						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
				return(ans)
			}
		}
		else {
			brkt <- bracket.bounded.root(score2,lower=-d+0.001,lbnd=-d,upper=2+d-0.001,ubnd=2+d,
								d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
			ans <- uniroot(score2,lower= brkt[1],upper=brkt[2],
						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
#			ans <- uniroot(score2,lower= -d+.Machine$double.eps^0.5,upper=2+d-.Machine$double.eps^0.5,
#						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
		    return(ans)
		}
	}
	if(d < 0){
			brkt <- bracket.bounded.root(score2,lower=-d+0.001,lbnd=-d,upper=2+d-0.001,ubnd=2+d,
								d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
		ans <- uniroot(score2, lower =  brkt[1],
					upper = brkt[2],
					tol = .Machine$double.eps^0.75, d = d, x1 = x1,
					m1 = m1, x2 = x2, m2 = m2, n1 = n1,
					n2 = n2)$root

#		ans <- uniroot(score2, lower =  - d + .Machine$double.eps^0.5,
#					upper = 2 + d - .Machine$double.eps^0.5,
#					tol = .Machine$double.eps^0.75, d = d, x1 = x1,
#					m1 = m1, x2 = x2, m2 = m2, n1 = n1,
#					n2 = n2)$root
	}
	else {
			brkt <- bracket.bounded.root(score2,lower=d+0.025,lbnd=d,upper=2-d-0.05,ubnd=2-d,
								d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
		ans <- uniroot(score2, lower = brkt[1], upper = brkt[2],
					tol = .Machine$double.eps^0.75,d = d, x1 = x1,	m1 = m1, x2 = x2, m2 = m2, n1 = n1,	n2 = n2)$root
#		ans <- uniroot(score2, lower = d + .Machine$double.eps^0.5, upper = 2 - d - .Machine$double.eps^0.5,
#					tol = .Machine$double.eps^0.75,d = d, x1 = x1,	m1 = m1, x2 = x2, m2 = m2, n1 = n1,	n2 = n2)$root
	}
	ans
}

"Shatd.vec" <-
function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
	nd <- length(d)
	shat <- vector(length=nd)
	dhat <- pooledbinom.diff.mle(x1,m1,x2,m2,n1,n2)
	for(i in 1:nd){
		 shat[i] <- Shatd(d[i],x1,m1,x2,m2,n1,n2)
	 }
	shat
}

# This is the profile likelihood ratio interval
"pooledbinom.diff.lrt.ci" <-
function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
{
	p1 <- pooledbinom.mle(x1, m1, n1)
	p2 <- pooledbinom.mle(x2, m2, n2)
	dhat <- p1 - p2
	shat <- p1 + p2

	if(sum(x1)==0 & sum(x2) == 0){
		# This is from Newcombe, p. 889
		# NOTE: exp(-3.84/2) = 0.1465
		N1 <- sum(m1*n1)
		N2 <- sum(m2*n2)
		g <- exp(-qchisq(1-alpha,1)/2)
		lcl <- -1 + g^(1/N1)
		ucl <- 1 - g^(1/N2)
		return(c(p1=p1,p2=p2,diff=dhat,lcl=lcl,ucl=ucl,alpha=alpha))
	}

	loglike.p <- function(p,x,m,n){
		if(p > 0) sum(x*log((1-(1-p)^m))) + log(1-p)*sum(m*(n-x))
		else 0
	}

	loglike.ds <- function(d,s,x1,m1,x2,m2,n1,n2)
		loglike.p((s+d)/2,x1,m1,n1) + loglike.p((s-d)/2,x2,m2,n2)

	lrt.stat <- function(d0,dhat,shat,x1,m1,x2,m2,n1,n2)
				{
					shat0 <- Shatd(d0,x1,m1,x2,m2,n1,n2)
					2*(loglike.ds(dhat,shat,x1,m1,x2,m2,n1,n2) -
						loglike.ds(d0,shat0,x1,m1,x2,m2,n1,n2))
				}

	f <- function(d0,dhat,shat,x1,m1,x2,m2,n1,n2,alpha)
					lrt.stat(d0,dhat,shat,x1,m1,x2,m2,n1,n2) - qchisq(1-alpha,1)

  	f.dhat <- f(dhat,dhat,shat,x1,m1,x2,m2,n1,n2,alpha)
	stepsize <- 10
 	for(i in 1:stepsize){
	  	if(dhat - i/stepsize <= -1){
		  	 lower.lim <- -1 + .Machine$double.eps
		  	 break
	  	 }
  		if(f(dhat - i/stepsize,dhat,shat,x1,m1,x2,m2,n1,n2,alpha)*f.dhat < 0){
	  		 lower.lim <- dhat - i/stepsize
	  		 break
  		 }
	}
	# should have upper = dhat.mpl

	lower <- uniroot(f, lower = lower.lim, upper = dhat , dhat=dhat,shat=shat,x1 = x1, m1 =
				m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,alpha=alpha)$root

  	for(i in 1:stepsize){
	  	if(dhat + i/stepsize >= 1){
		  	 upper.lim <- 1 - .Machine$double.eps
		  	 break
	  	 }
  		if(f(dhat + i/stepsize,dhat,shat,x1,m1,x2,m2,n1,n2,alpha)*f.dhat < 0){
	  		 upper.lim <- dhat + i/stepsize
	  		 break
  		 }
	}
	# should have lower = dhat.mpl
 	upper <- uniroot(f, lower = dhat, upper = upper.lim, dhat=dhat,shat=shat,x1 = x1, m1 =
 				m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,alpha=alpha)$root

	c(p1=p1,p2=p2,diff=dhat,lower=lower,upper=upper,alpha=alpha)
}

"pooledbinom.diff.lrtstat" <-
function(d,x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
{
	p1 <- pooledbinom.mle(x1, m1, n1)
	p2 <- pooledbinom.mle(x2, m2, n2)
	dhat <- p1 - p2
	shat <- p1 + p2

	loglike.p <- function(p,x,m,n){
		if(p > 0) sum(x*log((1-(1-p)^m))) + log(1-p)*sum(m*(n-x))
		else 0
	}

	loglike.ds <- function(d,s,x1,m1,x2,m2,n1,n2)
		loglike.p((s+d)/2,x1,m1,n1) + loglike.p((s-d)/2,x2,m2,n2)


	lrt.stat <- function(d0,dhat,shat,x1,m1,x2,m2,n1,n2)
				{
					shat0 <- Shatd(d0,x1,m1,x2,m2,n1,n2)
					2*(loglike.ds(dhat,shat,x1,m1,x2,m2,n1,n2) -
						loglike.ds(d0,shat0,x1,m1,x2,m2,n1,n2))
				}
	lrt.stat(d,dhat,shat,x1,m1,x2,m2,n1,n2)
}

"pooledbinom.diff.lrtstat.vec" <-
function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
	n <- length(d)
	lmp <- vector(length=n)
	for(i in 1:n){
		lmp[i] <- pooledbinom.diff.lrtstat(d[i],x1,m1,x2,m2,n1,n2)
	}
	lmp
}


# maximum likelihood estimate
"pooledbinom.diff.mle" <-
function(x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
	pooledbinom.mle(x1,m1,n1) - pooledbinom.mle(x2,m2,n2)

"pooledbinom.diff.newcombe.ci" <-
function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha = 0.05)
{
	s1 <- as.vector(pooledbinom.score.ci(x1, m1, n1, alpha = alpha))
	# as.vector to drop names
	p1 <- s1[1]
	L1 <- s1[2]
	U1 <- s1[3]
	s2 <- as.vector(pooledbinom.score.ci(x2, m2, n2, alpha = alpha))
	p2 <- s2[1]
	L2 <- s2[2]
	U2 <- s2[3]
	thetahat <- p1 - p2
	delta <- sqrt((p1 - L1)^2 + (U2 - p2)^2)
	epsilon <- sqrt((U1 - p1)^2 + (p2 - L2)^2)
	c(p1 = p1, p2 = p2, diff = thetahat, lower = thetahat - delta,
		upper = thetahat + epsilon, alpha = alpha)
}

"pooledbinom.diff.newcombecscore.ci" <-
function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha = 0.05)
{
	s1 <- as.vector(pooledbinom.cscore.ci(x1, m1, n1, alpha = alpha))
	# as.vector to drop names
	p1 <- s1[1]
	L1 <- s1[2]
	U1 <- s1[3]
	s2 <- as.vector(pooledbinom.cscore.ci(x2, m2, n2, alpha = alpha))
	p2 <- s2[1]
	L2 <- s2[2]
	U2 <- s2[3]
	thetahat <- p1 - p2
	delta <- sqrt((p1 - L1)^2 + (U2 - p2)^2)
	epsilon <- sqrt((U1 - p1)^2 + (p2 - L2)^2)
	c(p1 = p1, p2 = p2, diff = thetahat, lower = thetahat - delta,
		upper = thetahat + epsilon, alpha = alpha)
}

"pooledbinom.diff.wald.ci" <-
function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha = 0.05)
{
	p1 <- pooledbinom.mle(x1, m1, n1)
	p1.var <- pooledbinom.mle.var(p1, m1, n1)
	p2 <- pooledbinom.mle(x2, m2, n2)
	p2.var <- pooledbinom.mle.var(p2, m2, n2)
	diff <- p1 - p2
	diff.var <- p1.var + p2.var
	z <- qnorm(1 - alpha/2)
	c(p1 = p1, p2 = p2, diff = diff, lower = diff - z * sqrt(diff.var),
		upper = diff + z * sqrt(diff.var), alpha = alpha)
}


#######################################################################################################################
#######################################################################################################################
#######################################################
# Score methods
#######################################################
#######################################################################################################################
#######################################################################################################################
Z <- function(d, s, x1, m1, x2, m2, n1=rep(1,length(x1)), n2=rep(1,length(x2)))
	 {
		0.5*(score.p((d + s)/2, x1, m1, n1) - score.p((s -
		d)/2, x2, m2, n2)) * sqrt((pooledbinom.mle.var(
		(d + s)/2, m1, n1) + pooledbinom.mle.var((
		s - d)/2, m2, n2)))
	 }

"score.p" <-
function(p, x, m, n = rep(1,length(m)))
{
	if(sum(x) == 0 & p >= 0) return(-sum(m*n)/(1-p))
	if(p <= 0 | p >= 1) return(0)
	sum(m*x/(1-(1-p)^m) - m*n)/(1-p)
}

"score2" <-
function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
	0.5*(score.p((s+d)/2,x1,m1,n1) + score.p((s-d)/2,x2,m2,n2))
}

"score2.vec" <-
function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
	ns <- length(s)
	s2 <- vector(length=ns)
	for(i in 1:ns) s2[i] <- score2(s[i],d,x1,m1,x2,m2,n1,n2)
	s2
}

"pooledbinom.diff.scorestat" <-
function(d, x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)))
{
	dhat <- pooledbinom.diff.mle(x1,m1,x2,m2,n1,n2)

	Z <- function(d, shat, x1, m1, x2, m2, n1, n2)
		 {
			if(sum(x1) == 0){
				g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
						{
							nd <- length(d)
							gv <- vector(length=nd)
							for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
							gv
						}
				dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
								tol = .Machine$double.eps,
							x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
				N1 <- sum(m1*n1)

				if(d <= dstar){
					0.5*(-sum(m1*n1)/(1-(shat+d)/2) * as.numeric(shat + d > 0)  - score.p((shat - d)/2, x2, m2, n2)) *
					sqrt(((1-(shat+d)/2)^2 / N1)  * as.numeric(shat+d > 0) +
						pooledbinom.mle.var((shat - d)/2, m2, n2))

				} else {
				0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
				sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2))
				}
			} else {
				#
				# This is the REGULAR (easy) case -- x <> 0
				#
			p1 <- (shat+d)/2
			p2 <- (shat-d)/2
			I1 <- pooledbinom.I(p1,m1,n1)
			I2 <- pooledbinom.I(p2,m2,n2)
			R <- I1 / (I1 + I2)
			V1 <- 1/I1
			V2 <- 1/I2
			S1 <- score.p(p1,x1,m1,n1)
			S2 <- score.p(p2,x2,m2,n2)
			ans <- (R*S1 - (1-R)*S2) * sqrt(V1 + V2)
			#cprint(ans)
			return(ans)
			#	0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
			#	sqrt((pooledbinom.mle.var((shat + d)/2, m1, n1) + pooledbinom.mle.var((
	 		#	shat - d)/2, m2, n2)))
			}
		 }

	shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
	Z(d, shat, x1, m1, x2, m2, n1, n2)  # - qnorm(1-0.05/2)
}


"pooledbinom.diff.scorestat.vec" <-
function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
{
	n <- length(d)
	sc <- vector(length=n)
	for(i in 1:n) {
		sc[i] <- pooledbinom.diff.scorestat(d[i],x1,m1,x2,m2,n1,n2)
	}
	sc
}

"pooledbinom.diff.score.ci" <-
function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
{
	p1 <- pooledbinom.mle(x1, m1, n1)
	p2 <- pooledbinom.mle(x2, m2, n2)
	dhat <- p1 - p2
	shat <- p1 + p2

	#if(p1 == 0 & p2 == 0){
	if(sum(x1) == 0 & sum(x2)==0){
		# Not quite sure of this...MISSING 1/2!?!
		N1 <- sum(m1*n1)
		N2 <- sum(m2*n2)
		z2 <- qnorm(1-alpha/2)^2
		#lower <- -z2/(N2 + z2)
		lower <- -as.vector(pooledbinom.score.ci(x2,m2,n2,alpha=alpha)[3]) # -upper limit from pop'n 2
		#upper <- z2/(N1 + z2)
		upper <- as.vector(pooledbinom.score.ci(x1,m1,n1,alpha=alpha)[3]) # upper limit from pop'n 1
		return(c(p1=0,p2=0,diff=0,lower=lower,upper=upper,alpha=alpha))
	}
	if(sum(x1) == sum(n1) | sum(x2) == sum(n2))
		return(c(p1=pooledbinom.mle(x1,m1,n1),p2=pooledbinom.mle(x2,m2,n2),diff=0,lower=-1,upper=1))

	Z <- function(d, shat, x1, m1, x2, m2, n1, n2)
		 {
			if(sum(x1) == 0){
				g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
						{
							nd <- length(d)
							gv <- vector(length=nd)
							for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
							gv
						}
				dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
								tol = .Machine$double.eps,
							x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
				N1 <- sum(m1*n1)
				if(d <= dstar){
					return(0.5*(score.p((shat+d)/2, x1, m1, n1) * as.numeric(shat + d > 0)  - score.p((shat - d)/2, x2, m2, n2)) *
					sqrt(((1-(shat+d)/2)^2 / N1)  * as.numeric(shat+d > 0) +
						pooledbinom.mle.var((shat - d)/2, m2, n2)))

				} else {
				return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
					sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
				}
			}
			if(sum(x2) == 0){
				g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
						{
							nd <- length(d)
							gv <- vector(length=nd)
							for(i in 1:nd) gv[i] <-  score2(d[i],d[i],x1,m1,x2,m2,n1,n2)
							gv
						}
				dstar <- uniroot(g,lower = 0 + .Machine$double.eps^0.5,upper = 1 -.Machine$double.eps^0.5,
								tol = .Machine$double.eps,
							x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
				N2 <- sum(m2*n2)
				if(d >= dstar){
					return(0.5*(score.p((shat+d)/2, x1, m1, n1)  -
						score.p((shat - d)/2, x2, m2, n2)*as.numeric(shat - d > 0)) *
					sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + (1-(shat-d)/2)^2/N2 * as.numeric(shat-d>0)))

				} else {
					return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
						sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
				}
			}
			#
			# This is the REGULAR (easy) case -- x <> 0
			#
			0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
			sqrt((pooledbinom.mle.var((shat + d)/2, m1, n1) + pooledbinom.mle.var((
 			shat - d)/2, m2, n2)))
		 }

 	f.lower <- function(d, x1, m1, x2, m2, n1, n2,dhat)
		 {
			shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
			z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - qnorm(1-alpha/2)
			z
  		 }
 	f.upper <- function(d, x1, m1, x2, m2, n1, n2,dhat)
		 {
			shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
			z <- Z(d,shat,x1,m1,x2,m2,n1,n2) + qnorm(1-alpha/2)
			z
  		 }

  	# Bound the root, in a rather slow but safe way
  	f.dhat <- f.lower(dhat,x1,m1,x2,m2,n1,n2,dhat)

	stepsize <- 10
 	for(i in 1:stepsize){
	  	if(dhat - i/stepsize <= -1){
		  	 lower.lim <- -1 + .Machine$double.eps
		  	 break
	  	 }
  		if(f.lower(dhat - i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
	  		 lower.lim <- dhat - i/stepsize
	  		 break
  		 }
	}
	lower <- uniroot(f.lower, lower = lower.lim, upper = dhat , x1 = x1, m1 =
				m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root
  	f.dhat <- f.upper(dhat,x1,m1,x2,m2,n1,n2,dhat)
  	for(i in 1:stepsize){
	  	if(dhat + i/stepsize >= 1){
		  	 upper.lim <- 1 - .Machine$double.eps
		  	 break
	  	 }
  		if(f.upper(dhat + i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
	  		 upper.lim <- dhat + i/stepsize
	  		 break
  		 }
	}
	#if(!exists("upper.lim")) upper.lim <- 1-.Machine$double.eps^0.5

 	upper <- uniroot(f.upper, lower = dhat, upper = upper.lim,x1 = x1, m1 =
 				m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root
	c(p1 = p1, p2 = p2, diff = dhat, lower = lower, upper = upper, alpha=alpha)
}


"pooledbinom.diff.cscore.ci" <-
function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
{
	p1 <- pooledbinom.mle(x1, m1, n1)
	p2 <- pooledbinom.mle(x2, m2, n2)
	dhat <- p1 - p2
	shat <- p1 + p2
	chi2 <- qnorm(1 - alpha/2)^2
	if(sum(x1) == 0 & sum(x2)==0){
		# Not quite sure of this...MISSING 1/2!?!
		N1 <- sum(m1*n1)
		N2 <- sum(m2*n2)
		z2 <- qnorm(1-alpha/2)^2
		lower <- -z2/(N2 + z2)
		upper <- z2/(N1 + z2)
		return(c(p1=0,p2=0,diff=0,lower=lower,upper=upper,alpha=alpha))
	}
	if(sum(x1) == sum(n1) | sum(x2) == sum(n2))
		return(c(p1=pooledbinom.mle(x1,m1,n1),p2=pooledbinom.mle(x2,m2,n2),diff=0,lower=-1,upper=1))

	Z <- function(d, shat, x1, m1, x2, m2, n1, n2)
		 {
			if(sum(x1) == 0){
				g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
						{
							nd <- length(d)
							gv <- vector(length=nd)
							for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
							gv
						}
				dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
								tol = .Machine$double.eps,
							x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
				N1 <- sum(m1*n1)
				if(d <= dstar){
					return(0.5*(score.p((shat+d)/2, x1, m1, n1) * as.numeric(shat + d > 0)  - score.p((shat - d)/2, x2, m2, n2)) *
					sqrt(((1-(shat+d)/2)^2 / N1)  * as.numeric(shat+d > 0) +
						pooledbinom.mle.var((shat - d)/2, m2, n2)))

				} else {
				return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
					sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
				}
			}
			if(sum(x2) == 0){
				g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
						{
							nd <- length(d)
							gv <- vector(length=nd)
							for(i in 1:nd) gv[i] <-  score2(d[i],d[i],x1,m1,x2,m2,n1,n2)
							gv
						}
				dstar <- uniroot(g,lower = 0 + .Machine$double.eps^0.5,upper = 1 -.Machine$double.eps^0.5,
								tol = .Machine$double.eps,
							x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
				N2 <- sum(m2*n2)
				if(d >= dstar){
					return(0.5*(score.p((shat+d)/2, x1, m1, n1)  -
						score.p((shat - d)/2, x2, m2, n2)*as.numeric(shat - d > 0)) *
					sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + (1-(shat-d)/2)^2/N2 * as.numeric(shat-d>0)))

				} else {
					return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
						sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
				}
			}
			#
			# This is the REGULAR (easy) case -- x <> 0
			#
			0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
			sqrt((pooledbinom.mle.var((shat + d)/2, m1, n1) + pooledbinom.mle.var((
 			shat - d)/2, m2, n2)))
		 }

	mu3 <- function(d,shat,m1,m2,n1,n2)
		   {
			   I1 <- pooledbinom.I((shat+d)/2,m1,n1)
			   I2 <- pooledbinom.I((shat-d)/2,m2,n2)
			   R <- I1/(I1+I2)
			   (1-R)^3*pooledbinom.mu3((shat+d)/2,m1,n1) - R^3*pooledbinom.mu3((shat-d)/2,m2,n2)
		   }
    gamma1 <- function(d,shat,m1,m2,n1,n2)
              {
	              mu3(d,shat,m1,m2,n1,n2)*
	               (pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat-d)/2,m2,n2))^(3/2)
              }

 	f.lower <- function(d, x1, m1, x2, m2, n1, n2,dhat)
		 {
			shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
			z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - gamma1(d,shat,m1,m2,n1,n2)*(chi2 - 1)/6 - qnorm(1-alpha/2)
			z
  		 }
 	f.upper <- function(d, x1, m1, x2, m2, n1, n2,dhat)
		 {
			shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
			z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - gamma1(d,shat,m1,m2,n1,n2)*(chi2 - 1)/6 + qnorm(1-alpha/2)
			z
  		 }

  	# Bound the root, in a rather slow but safe way
  	f.dhat <- f.lower(dhat,x1,m1,x2,m2,n1,n2,dhat)

	stepsize <- 10
 	for(i in 1:stepsize){
	  	if(dhat - i/stepsize <= -1){
		  	 lower.lim <- -1 + .Machine$double.eps
		  	 break
	  	 }
  		if(f.lower(dhat - i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
	  		 lower.lim <- dhat - i/stepsize
	  		 break
  		 }
	}
	lower <- uniroot(f.lower, lower = lower.lim, upper = dhat , x1 = x1, m1 =
				m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root

  	f.dhat <- f.upper(dhat,x1,m1,x2,m2,n1,n2,dhat)
  	for(i in 1:stepsize){
	  	if(dhat + i/stepsize >= 1){
		  	 upper.lim <- 1 - .Machine$double.eps
		  	 break
	  	 }
  		if(f.upper(dhat + i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
	  		 upper.lim <- dhat + i/stepsize
	  		 break
  		 }
	}

 	upper <- uniroot(f.upper, lower = dhat, upper = upper.lim,x1 = x1, m1 =
 				m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root


	c(p1=p1,p2=p2,diff = dhat, lower = lower, upper = upper,alpha=alpha)
}


"pooledbinom.diff.bcscore.ci" <-
function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
{
	p1 <- pooledbinom.mle(x1, m1, n1)
	p2 <- pooledbinom.mle(x2, m2, n2)
	dhat <- p1 - p2
	shat <- p1 + p2
	chi2 <- qnorm(1 - alpha/2)^2
		if(sum(x1) == 0 & sum(x2)==0){
		# Not quite sure of this...MISSING 1/2!?!
		N1 <- sum(m1*n1)
		N2 <- sum(m2*n2)
		z2 <- qnorm(1-alpha/2)^2
		lower <- -z2/(N2 + z2)
		upper <- z2/(N1 + z2)
		return(c(p1=0,p2=0,diff=0,lower=lower,upper=upper,alpha=alpha))
	}
	if(sum(x1) == sum(n1) | sum(x2) == sum(n2))
		return(c(p1=pooledbinom.mle(x1,m1,n1),p2=pooledbinom.mle(x2,m2,n2),diff=0,lower=-1,upper=1))

	Z <- function(d, shat, x1, m1, x2, m2, n1, n2)
		 {
			if(sum(x1) == 0){
				g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
						{
							nd <- length(d)
							gv <- vector(length=nd)
							for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
							gv
						}
				dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
								tol = .Machine$double.eps,
							x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
				N1 <- sum(m1*n1)
				if(d <= dstar){
					return(0.5*(score.p((shat+d)/2, x1, m1, n1) * as.numeric(shat + d > 0)  - score.p((shat - d)/2, x2, m2, n2)) *
					sqrt(((1-(shat+d)/2)^2 / N1)  * as.numeric(shat+d > 0) +
						pooledbinom.mle.var((shat - d)/2, m2, n2)))

				} else {
				return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
					sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
				}
			}
			if(sum(x2) == 0){
				g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
						{
							nd <- length(d)
							gv <- vector(length=nd)
							for(i in 1:nd) gv[i] <-  score2(d[i],d[i],x1,m1,x2,m2,n1,n2)
							gv
						}
				dstar <- uniroot(g,lower = 0 + .Machine$double.eps^0.5,upper = 1 -.Machine$double.eps^0.5,
								tol = .Machine$double.eps,
							x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
				N2 <- sum(m2*n2)
				if(d >= dstar){
					return(0.5*(score.p((shat+d)/2, x1, m1, n1)  -
						score.p((shat - d)/2, x2, m2, n2)*as.numeric(shat - d > 0)) *
					sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + (1-(shat-d)/2)^2/N2 * as.numeric(shat-d>0)))

				} else {
					return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
						sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
				}
			}
			#
			# This is the REGULAR (easy) case -- x <> 0
			#
			0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) *
			sqrt((pooledbinom.mle.var((shat + d)/2, m1, n1) + pooledbinom.mle.var((
 			shat - d)/2, m2, n2)))
		 }
	mu3 <- function(d,shat,m1,m2,n1,n2)
		   {
			   I1 <- pooledbinom.I((shat+d)/2,m1,n1)
			   I2 <- pooledbinom.I((shat-d)/2,m2,n2)
			   R <- I1/(I1+I2)
			   (1-R)^3*pooledbinom.mu3((shat+d)/2,m1,n1) - R^3 * pooledbinom.mu3((shat-d)/2,m2,n2)
		   }
    gamma1 <- function(d,shat,m1,m2,n1,n2)
              {
	              mu3(d,shat,m1,m2,n1,n2)*
	               (pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat-d)/2,m2,n2))^(3/2)
              }
	bias <- function(d,shat,m1,m2,n1,n2)
	        {
			   I1 <- pooledbinom.I((shat+d)/2,m1,n1)
			   I2 <- pooledbinom.I((shat-d)/2,m2,n2)
			   R <- I1/(I1+I2)
		       R^(3/2) * sqrt(I2) * pooledbinom.bias((shat+d)/2,m1,n1) -
		         (1-R)^(3/2) * sqrt(I1) * pooledbinom.bias((shat-d)/2,m2,n2)
	        }

  	f.lower <- function(d, x1, m1, x2, m2, n1, n2,dhat)
		 {
			shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
			z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - bias(d,shat,m1,m2,n1,n2) - gamma1(d,shat,m1,m2,n1,n2)*(chi2 - 1)/6 - qnorm(1-alpha/2)
			z
  		 }
 	f.upper <- function(d, x1, m1, x2, m2, n1, n2,dhat)
		 {
			shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
			z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - bias(d,shat,m1,m2,n1,n2) - gamma1(d,shat,m1,m2,n1,n2)*(chi2 - 1)/6 + qnorm(1-alpha/2)
			z
  		 }

  	# Bound the root, in a rather slow but safe way
  	f.dhat <- f.lower(dhat,x1,m1,x2,m2,n1,n2,dhat)

	stepsize <- 10
 	for(i in 1:stepsize){
	  	if(dhat - i/stepsize <= -1){
		  	 lower.lim <- -1 + .Machine$double.eps
		  	 break
	  	 }
  		if(f.lower(dhat - i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
	  		 lower.lim <- dhat - i/stepsize
	  		 break
  		 }
	}
	lower <- uniroot(f.lower, lower = lower.lim, upper = dhat , x1 = x1, m1 =
				m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root

  	f.dhat <- f.upper(dhat,x1,m1,x2,m2,n1,n2,dhat)
  	for(i in 1:stepsize){
	  	if(dhat + i/stepsize >= 1){
		  	 upper.lim <- 1 - .Machine$double.eps
		  	 break
	  	 }
  		if(f.upper(dhat + i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
	  		 upper.lim <- dhat + i/stepsize
	  		 break
  		 }
	}

 	upper <- uniroot(f.upper, lower = dhat, upper = upper.lim,x1 = x1, m1 =
 				m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root

	c(p1=p1,p2=p2,diff = dhat, lower = lower, upper = upper,alpha=alpha)
}


########################
# MIR
########################

"pooledbinom.diff.mir.ci" <-
function(x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)),alpha=0.05)
{
	N1 <- sum(m1*n1)
	N2 <- sum(m2*n2)
	mir1 <- sum(x1)/N1
	mir2 <- sum(x2)/N2
	mir1.var <- mir1*(1-mir1)/N1
	mir2.var <- mir2*(1-mir2)/N2
	mir.diff <- mir1 - mir2
	mir.diff.stderr <- sqrt(mir1.var + mir2.var)
	z <- qnorm(1-alpha/2)
	c(p1=mir1,p2=mir2,diff=mir.diff,lower=mir.diff - z * mir.diff.stderr,upper=mir.diff + z*mir.diff.stderr,alpha=alpha)
}


#########################################################################
# Some utility functions
#########################################################################
"bracket.root" <-
function(f, lower=0, upper=1, ..., scaling = 1.6, direction = c("both", "up",
	"down"), max.iter = 50)
{
	direction <- match.arg(direction)
	if(lower == upper)
		stop("lower cannot equal upper")
	f1 <- f(lower, ...)
	f2 <- f(upper, ...)
	done <- F
	num.iter <- 0
	while(!done) {
		#cprint(lower)
		#cprint(upper)
		num.iter <- num.iter + 1
		if(num.iter >= max.iter)
			stop(paste("root not bracketed in", max.iter,
				"iterations\n"))
		if(f1 * f2 < 0) {
			done <- T
			ans <- c(lower, upper)
		}
		switch(direction,
			both = if(abs(f1) < abs(f2)) f1 <- f(lower <- lower +
					scaling * (lower - upper), ...) else f2 <-
					f(upper <- upper + scaling * (upper -
					lower), ...),
			down = f1 <- f(lower <- lower + scaling * (lower -
				upper), ...),
			up = f2 <- f(upper <- upper + scaling * (upper - lower),
				...))
	}
	sort(ans)
}



# Use this function to bracket a root by lower and upper
# that is bounded below by lbnd and above by ubnd
# idea is to start with (lower,upper) as brackets, then
# step toward the appropriate bound, lbnd or ubnd,
# by moving the the midpoint between lower and lbnd
# (upper and ubnd) updating lower (upper) each time
"bracket.bounded.root" <-
function(f,lower=0,upper=1,lbnd=0,ubnd=1,...,slicing=2,max.iter=50)
{
	if(lower == upper) stop("lower cannot equal upper")
	f1 <- f(lower,...)
	f2 <- f(upper,...)
	done <- F
	num.iter <- 0
	while(!done){
		num.iter <- num.iter + 1
		if(num.iter >= max.iter) stop(paste("root not bracketed in",max.iter,"iterations\n"))
		if(f1*f2 < 0) {
			done <- T
			ans <- c(lower,upper)
		}
		# march toward the bound in slicing steps
		if(abs(f1) < abs(f2))
			f1 <- f(lower <- (lower+lbnd)/slicing,...)
		else
			f2 <- f(upper <- (upper+ubnd)/slicing,...)
	}
	sort(ans)
}


"printc" <-
function(x,digits=options()$digits)
{
	namex <- deparse(substitute(x))
	cat(paste(namex,"=",round(x,digits=digits),"\n"))
}

