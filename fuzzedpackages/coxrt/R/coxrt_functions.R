
# assumes that t is a scalar
# get S(t)=Pr(T>t) from (time,surv)
.get.S <- function(t, time, surv)	{
  n <- length(time)
  if (t<time[1] && abs(t-time[1])>100*.Machine$double.eps) return(1)
  if (t>time[n] && abs(t-time[n])>100*.Machine$double.eps) return(surv[n])
  return( surv[sum(time<=t | abs(time-t)<100*.Machine$double.eps)] )
} # due to the necessarily limited precision of the underlying element types;


.get_est <- function(X, T, Z, delta)
{
  # compute weights using NA estimate:
  fit <- summary(survfit(Surv(time=X, time2=T, event=delta)~1))
  ti <- fit$time # time of events
  died <- fit$n.event # number of events
  n.ar <- fit$n.risk # number at risk
  H <- cumsum(died/n.ar)
  S.NA <- exp(-H)
  w.simple <- sapply(X, .get.S, time=ti, surv=S.NA)
  # estimate beta
  srv.trunc <- Surv(X, event=delta)
  sol <- try(coxph(srv.trunc ~ Z +offset(-log(w.simple))) )
  if (class(sol)!="coxph") {
    cat("coxph.RT: Error occurred. ", sol,"\n")
    return(NULL)
  }
  return(list(est=sol$coef, w=w.simple))
}


#' Fits Cox Regression Model Using Right Truncated Data
#'
#'
#' Estimates covariate effects in a Cox proportional hazard regression
#' from right-truncated survival data assuming positivity, that is
#' \code{P(lifetime>max(right) | Z=0)=0}.
#'
#' When positivity does not hold, the estimator of regression coefficients
#' will be biased.
#' But if all the covariates are independent in the population,
#' the Wald test performed by this function is still valid and can be used
#' for testing partial hypotheses about regression coefficients
#' even in the absence of positivity. If the covariates are not independent and
#' positivity does not hold, the partial tests cannot guarantee the correct
#' level of type I error.
#'
#' @param formula a formula object, with the response on the left of a ~ operator, and covariates on the right.
#' The response is a target lifetime variable.
#' @param right a right truncation variable.
#' @param data a data frame that includes the variables used in both sides of \code{formula}
#' and in \code{right}.
#' The observations with missing values in one of the variables are dropped.
#' @param bs logical value: if \code{TRUE}, the bootstrap esimator of standard error,
#' confidence interval,
#' and confidence upper and lower limits for one-sided confidence intervals
#' based on the bootstrap distribution are calculated. The default value is \code{FALSE}.
#' @param nbs.rep number of bootstrap replications. The default number is 500.
#' @param conf.int The confidence level for confidence intervals and hypotheses tests.
#' The default level is 0.95.
#'
#' @return  A list with components:
#' \tabular{llr}{
#' \code{coef} \tab an estimate of regression coefficients \tab   \cr
#' \code{var} \tab covariance matrix of estimates of regression coefficients based on the analytic formula\tab  \cr
#' \code{n} \tab the number of observations used to fit the model \tab   \cr
#' \code{summary} \tab a data frame with a summary of fit: \tab  \cr}
#' \itemize{
#' \item{\code{coef}} a vector of coefficients
#' \item{\code{exp.coef}} exponent of regression coefficients (=hazard ratio)
#' \item{\code{SE}} asymptotic standard error estimate based on the analytic formula derived in Vakulenko-Lagun et al. (2018)
#' \item{\code{CI.L}} lower confidence limit for two-sided hypothesis \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>} = 0}{\eqn{\beta_i=0}}
#' \item{\code{CI.U}} upper confidence limit for two-sided hypothesis \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>} = 0}{\eqn{\beta_i=0}}
#' \item{\code{pvalue}} p-value from a Wald test for a two-sided hypothesis
#' \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:}  \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>} = 0}{\eqn{\beta_i=0}}
#' \item{\code{pvalue.H1.b.gr0}} p-value from the Wald test for a one-sided
#' partial hypothesis \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:}  \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>}\eqn{\le 0}}{\eqn{\beta_i\le 0}}
#' based on the analytical asymptotic standard error estimate
#' \item{\code{pvalue.H1.b.le0}} p-value from the Wald test a for one-sided
#' partial hypothesis \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>}\eqn{\ge 0}}{\eqn{\beta_i\ge 0}}
#' based on the analytical asymptotic standard error estimate }
#' \tabular{ll}{
#' \code{bs } \tab if the input argument \code{bs} was TRUE, then an output list also includes an element \code{bs} with\cr
#'  \tab statistics from the bootstrap distribution of estimated
#'  coefficients:\cr}
#'  \itemize{
#' \item{\code{num.bs.rep}}
#' {number of bootsrap replications used to obtain the sample distribution}
#' \item{\code{var}} {estimated variance}
#' \item{\code{summary}} {a data frame with a summary
#' of bootstrap distribution that includes:
#'  \code{SE}, a bootstrap estimated standard error;
#'  \code{CI.L},  a quantile estimated lower confidence limit
#'  for two-sided hypothesis \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>} = 0}{\eqn{\beta_i=0}};
#'  \code{CI.U}, a  quantile estimated upper confidence limit for two-sided hypothesis
#' \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:}  \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>} = 0}{\eqn{\beta_i=0}};
#' \code{CI.L.H1.b.gr0},
#' a quantile estimated the limit for one-sided hypothesis
#' \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>}\eqn{\le 0}}{\eqn{\beta_i\le 0}};
#' \code{CI.U.H1.b.le0}, a
#' quantile estimated the limit for one-sided hypothesis
#' \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>}\eqn{\ge 0}}{\eqn{\beta_i\ge 0}}.}
#'}
#'
#'
#' @seealso \code{\link{coxph.RT.a0}}, \code{\link{coxrt}}, \code{\link[survival]{coxph}}
#' @examples
#' # loading AIDS data set
#' library(gss)
#' data(aids)
#' all <- data.frame(age=aids$age, ageg=as.numeric(aids$age<=59), T=aids$incu,
#' R=aids$infe, hiv.mon =102-aids$infe)
#' all$T[all$T==0] <- 0.5 # as in Kalbfeisch and Lawless (1989)
#' s <- all[all$hiv.mon>60,] # select those who were infected in 1983 or later
#' # analysis assuming positivity
#' # we request bootstrap SE estimate as well:
#' sol <- coxph.RT(T~ageg, right=R, data=s, bs=FALSE)
#' sol
#' sol$summary # print the summary of fit based on the analytic Asymptotic Standard Error estimate
#'
#' @export
coxph.RT <- function(formula, right, data, bs=FALSE, nbs.rep=500, conf.int=0.95) {
  Call <- match.call()
  mnames <- c("", "formula", "data", "right")
  cnames <- names(Call) #  formula right data bs
  cnames <- cnames[match(mnames, cnames, 0)] # formula data right
  mcall <- Call[cnames]
  mcall[[1]] <- as.name("model.frame")
  obj <- eval(mcall, parent.frame())
  d <- data[, c(names(obj)[-length(names(obj))], as.character(Call[cnames[length(cnames)]]))]
  X <- d[,1]
  T <- d[,ncol(d)]
  if ((ncol(d)-2)!=1)
    Z <- as.matrix( d[,-c(1, ncol(d))])
  else
    Z <- matrix( d[,-c(1, ncol(d))], nrow=nrow(d), ncol=ncol(d)-2)

  covNames <- colnames(d)[-c(1,ncol(d))]
  colnames(Z) <- covNames

  # removing missing observations:
  na.i <- ( is.na(X) | is.na(T) | apply(is.na(Z), 1, any))
  if (sum(na.i)>0)
  {
    cat(sum(na.i), " observations were deleted due to missing values.\n")
    X <- X[!na.i]
    T <- T[!na.i]
    Z <- matrix( Z[!na.i,], nrow=length(X), ncol= ncol(d)-2)
  }

  n <- length(X)
  delta <- rep(1, n)

  sol <- .get_est(X, T, Z, delta)

  if (is.null(sol)) {
    return(NULL)
  }
  est <- matrix(sol$est, nrow=length(covNames), ncol=1)
  rownames(est) <- covNames

  # estimate asymptotic variance(beta):
  exp.bZ <- exp(Z %*% est)
  Var <- getVar(exp.bZ, X, T,  Z,  sol$w)
  rownames(Var) <- colnames(Var) <- covNames

  # estimate variance(beta) by simple bootstrap:
  if (bs) {
    b.b <-  NULL
    for (b in 1:nbs.rep)
    {
       samp.b <- sample(n, size=n, replace=TRUE)
       est.b <- .get_est(X[samp.b], T[samp.b], Z[samp.b,], delta)$est
       if (is.null(est.b)) next
       b.b <- rbind(b.b, est.b)
    } # bs loop
    colnames(b.b) <- covNames
    se.bs <- apply(b.b, 2, sd)
    var.bs <- var(b.b)
    left.lim <- apply(b.b, 2, quantile, prob=(1-conf.int)/2, na.rm=TRUE)
    right.lim <- apply(b.b, 2, quantile, prob=1-(1-conf.int)/2, na.rm=TRUE)
    # the next row is a CI of [b.bar-z_alpha*SE, +infty)  is for H0:b<=0, H1:b>0
    left.H1.b.gr0 <- apply(b.b, 2, quantile, prob=(1-conf.int), na.rm=TRUE)
    # the next row is a CI of (-infty, b.bar+z_alpha*SE]  is for H0:b>=0, H1:b<0
    right.H1.b.le0 <- apply(b.b, 2, quantile, prob=conf.int, na.rm=TRUE)
  }
  names(sol$est) <- covNames
  out <- list()
  se <- sqrt(diag(Var))
  qn <- qnorm(1-(1-conf.int)/2)
  stat <- data.frame(coef=sol$est,
                     exp.coef=exp(sol$est),
                     SE=se,
                     CI.L=sol$est - qn*se,
                     CI.U=sol$est + qn*se,
                     pvalue=2*pnorm(-abs(sol$est)/se),
                     pvalue.H1.b.gr0=pnorm(sol$est/se, lower.tail=FALSE),
                     pvalue.H1.b.le0=pnorm(sol$est/se, lower.tail=TRUE)
  )# summary table
  rownames(stat) <- covNames
  out$Call <- Call
  out$summary <- stat
  out$coef <- sol$est
  out$var <- Var
  out$n <- n
  class(out) <- "coxph.RT"
  if (bs){
    out$bs$num.bs.rep <- nrow(b.b)
    out$bs$var <- var.bs
    out$bs$summary <- data.frame(SE=se.bs, CI.L=left.lim, CI.U=right.lim, CI.L.H1.b.gr0=left.H1.b.gr0,
                                 CI.U.H1.b.le0=right.H1.b.le0)
  } else {
    out$bs <- NULL
  }

  return(out)
}


# the equation that we want to solve:
.EE.RT <- function(par, T, Z, R, a0, W)
{
  n <- length(T)
  p <- length(par)
  par <- as.matrix(par, p,1)
  Ti <- matrix(T, n, n, byrow=TRUE)
  Tj <- matrix(T, n, n, byrow=FALSE)
  Wm <- matrix(W, n, n, byrow=FALSE)
  Cj <- matrix( c(a0^exp(Z %*% par)), n,n, byrow=FALSE)
  tmp <- matrix( c(exp(Z %*% par)), n,n, byrow=FALSE)/Wm*((Tj>=Ti) + Cj/(1-Cj))
  p2 <- NULL
  for (q in 1:p)
  {
    Zm <- matrix( Z[,q], n, n, byrow=FALSE)
    p2 <- c(p2, sum(colSums(Zm*tmp)/colSums(tmp)) )
  }
  colSums(Z) - p2
}

#' Fits Cox Regression with Adjustment for the Lack of Positivity
#'
#' Estimates covariate effects in a Cox proportional hazard regression from
#' right truncated survival data for a given value of \code{a0=P(lifetime>max(right) | Z=0)}.
#'   This probability reflects the chance of falling to the right of the upper bound
#'   of the support of the right truncation variable in the reference stratum
#'   where all the covariates are zero. Right truncation might result
#'   in a completely unobserved right tail of the distribution of the target lifetime.
#'   That means that it can happen there will be no "representatives" in a sample
#'   from the right tail.  Ignoring this and using \code{\link{coxph.RT}} in general
#'   will result in biased estimation of regression coefficients, whereas
#'   \code{coxph.RT.a0} allows to account for this violation.
#'
#' @param formula a formula object, with the response on the left of a ~ operator, and covariates on the right.
#' The response is a target lifetime variable.
#' @param right a right truncation variable.
#' @param data a data frame that includes the variables used in \code{formula} and in \code{right}.
#' @param a0 probability of falling into the unobservable region in the stratum of \code{Z=0},
#'  i.e. \code{P(lifetime>max(right) | Z=0)}. By default \code{a0=0}, which is equivalent to assuming positivity.
#' @param bs logical value: if TRUE, the bootstrap esimator of standard error, confidence interval,
#' and confidence upper and lower limits for one-sided confidence intervals
#' based on the bootstrap distribution are calculated. The default value is FALSE.
#' @param nbs.rep number of bootstrap replications. The default number is 200.
#' @param conf.int The confidence level for confidence intervals and hypotheses tests.
#' The default level is 0.95.
#'
#' @return  a list with components:
#' \tabular{ll}{
#' \code{convergence} \tab convergence code as returned by \code{\link[BB:BBsolve]{BBsolve}}. \cr
#' \tab \code{convergence} > 0 means that the algorithm diverged and a solution was not found. \cr
#' \tab \code{BBsolve} is used with a default parameters setting. \cr
#' \code{coef} \tab a vector of estimated regression coefficients. \cr
#' \code{var} \tab covariance matrix of regression coefficients, if the input argument \code{bs} was \code{TRUE};\cr
#' \tab  \code{NULL}, otherwise. \cr
#' \code{n} \tab the number of observations used to fit the model. \cr
#' \code{a0} \tab plugged-in value of \code{a0}. \cr
#' \code{bs} \tab if the input argument \code{bs} was TRUE, then an output list also includes an element \code{bs}\cr
#'  \tab  with statistics from the bootstrap distribution of estimated coefficients:\cr}
#' \itemize{
#' \item{\code{num.bs.rep}}
#' number of successful bootsrap replications used to obtain the sample distribution
#' \item{\code{var}} estimated variance of regression coefficients
#' \item{\code{summary}} a data frame with a summary
#' of bootstrap distribution that includes:
#'  \code{coef}, a vector of estimated regression coefficients;
#'  \code{exp.coef}, an exponent of regression coefficients (=hazard ratio);
#'  \code{SE}, a bootstrap estimated standard error;
#'  \code{CI.L}, a quantile estimated lower confidence limit
#'  for two-sided hypothesis \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>} = 0}{\eqn{\beta_i=0}};
#'  \code{CI.U}, a quantile estimated upper confidence limit for two-sided hypothesis
#' \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>} = 0}{\eqn{\beta_i=0}};
#' \code{CI.L.H1.b.gr0}, a
#' quantile estimated the limit for one-sided hypothesis
#' \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>}\eqn{\le 0}}{\eqn{\beta_i\le 0}};
#' \code{CI.U.H1.b.le0}, a
#' quantile estimated the limit for one-sided hypothesis
#' \ifelse{html}{\out{H<sub>0</sub>:}}{\eqn{H_0}:} \ifelse{html}{\eqn{\beta}\out{<sub>i</sub>}\eqn{\ge 0}}{\eqn{\beta_i\ge 0}}. }
#'
#'
#'
#' @seealso \code{\link{coxph.RT}}, \code{\link[BB:BBsolve]{BBsolve}}
#'
#' @examples
#' # loading AIDS data set
#' library(gss)
#' data(aids)
#' all <- data.frame(age=aids$age, ageg=as.numeric(aids$age<=59), T=aids$incu,
#' R=aids$infe, hiv.mon =102-aids$infe)
#' all$T[all$T==0] <- 0.5 # as in Kalbfeisch and Lawless (1989)
#' s <- all[all$hiv.mon>60,] # select those who were infected in 1983 or later
#'
#' # analysis using adjusted estimating equations for a0=0.2
#' sol.02 <- try(coxph.RT.a0(T~ageg, right=R, data=s, a0=0.2, bs=FALSE))
#' sol.02
#
#' # for a0=0
#' sol <- try(coxph.RT(T~ageg, right=R, data=s, bs=FALSE) )
#' sol$summary # print the summary of fit based on the asymptotic SE estimate
#'
#'
#' # sensitivity analysis for different values of a0
#' a_ <- seq(0.05, 0.55, by=0.05)
#' est <- NULL
#'
#' for(q in 1:length(a_))
#' {
#'   sol.a <- try(coxph.RT.a0(T~ageg, right=R, data=s, a0=a_[q], bs=FALSE))
#'   if (sol.a$convergence!=0)
#'   {
#'     cat("a0=", a_[q], ". Error occurred in BBsolve.\n")
#'   } else
#'   {
#'     cat("a=", a_[q]," ", " IPW.adj.est=", sol.a$coef, "\n")
#'     est <- c(est, sol.a$coef)
#'   }
#' }
#' require(ggplot2)
#' res.d <- data.frame(a0=c(0, a_), beta=c(sol$coef, est))
#'
#' p <- ggplot(res.d, aes(x=a0, y=beta)) +
#'   geom_line() + geom_point() +
#'   geom_hline(yintercept=0)
#' p + xlab(expression( paste(a[0], "=P(T>", r['*']," | z=0)" , sep="")) )+
#'   ylab(expression( paste(hat(beta), "(", a[0], ")" , sep="")) ) +
#'   scale_x_continuous(breaks=res.d$a0, labels=res.d$a0) +
#'   theme(axis.text.x = element_text(face="bold", angle=45),
#'         axis.text.y = element_text(face="bold"))
#'
#' @export
coxph.RT.a0 <-function(formula, right, data, a0=0, bs=FALSE, nbs.rep=200, conf.int=0.95)
{
  Call <- match.call()
  mnames <- c("", "formula", "data", "right")
  cnames <- names(Call) #  formula right data bs
  cnames <- cnames[match(mnames, cnames, 0)] # formula data right
  mcall <- Call[cnames]
  mcall[[1]] <- as.name("model.frame")
  obj <- eval(mcall, parent.frame())
  d <- data[, c(names(obj)[-length(names(obj))], as.character(Call[cnames[length(cnames)]]))]
  X <- d[,1]
  T <- d[,ncol(d)]
  if ((ncol(d)-2)!=1)
    Z <- as.matrix( d[,-c(1, ncol(d))])
  else
    Z <- matrix( d[,-c(1, ncol(d))], nrow=nrow(d), ncol=ncol(d)-2)
  covNames <- colnames(d)[-c(1,ncol(d))]
  colnames(Z) <- covNames

  # removing missing observations:
  na.i <- ( is.na(X) | is.na(T) | apply(is.na(Z), 1, any))
  if (sum(na.i)>0)
  {
    cat(sum(na.i), " observations were deleted due to missing values.\n")
    X <- X[!na.i]
    T <- T[!na.i]
    Z <- matrix( Z[!na.i,], nrow=length(X), ncol= ncol(d)-2)
  }

  n <- length(X)
  delta <- rep(1, n)


  # compute weights using NA estimate:
  fit <- summary(survival::survfit(survival::Surv(time=X, time2=T, event=delta)~1))
  ti <- fit$time # time of events
  died <- fit$n.event # number of events
  n.ar <- fit$n.risk # number at risk
  H <- cumsum(died/n.ar)
  S.NA <- exp(-H)
  w.simple <- sapply(X, .get.S, time=ti, surv=S.NA)
  # get the initial estimate of beta for BBsolve
  srv.trunc <- survival::Surv(X, event=delta)
  sol <- try(survival::coxph(srv.trunc ~ Z + offset(-log(w.simple)), data = obj) )
  if (class(sol)!="coxph") {
    cat("coxph.RT: Error occurred. ", sol,"\n")
  }
  names(sol$coef) <- covNames

  out <- list()
  out$Call <- Call
  out$coef <- NULL
  out$var <- NULL
  out$n <- n
  out$a0 <- a0
  class(out) <- "coxph.RT.a0"

  sol.a0 <- try(BBsolve(par=sol$coef, fn=.EE.RT, T=X, Z=Z, R=T, a0=a0, W=w.simple, quiet=TRUE))
  if (sol.a0$convergence!=0)
  {
    cat("RT sensitivity: Error occurred in BBsolve for a0=",a0 , "\n")
    out$convergence <- sol.a0$convergence
    return(out)
  }

  if (bs)
  {
    act.num <- 0
    b.b <-  NULL
    for (b in 1:nbs.rep)
    {
      samp.b <- sample(n, size=n, replace=TRUE)
      fit <- summary(survival::survfit(survival::Surv(time=X[samp.b], time2=T[samp.b], event=delta)~1))
      ti <- fit$time # time of events
      died <- fit$n.event # number of events
      n.ar <- fit$n.risk # number at risk
      H <- cumsum(died/n.ar)
      S.NA <- exp(-H)
      re.w.simple <- sapply(X[samp.b], .get.S, time=ti, surv=S.NA)
      srv.trunc.re <- survival::Surv(X[samp.b], event=delta)
      est.b <- try(survival::coxph(srv.trunc.re ~ Z[samp.b,] + offset(-log(re.w.simple))) )
      if (class(est.b)!="coxph")
      {
        cat("bs rep=",b," Error occurred in bootstrap resampling: ", est.b,"\n")
        next
      }

      sol.b <- try(BBsolve(par=est.b$coef, fn=.EE.RT , T=X[samp.b], Z=matrix(Z[samp.b,], n, ncol(Z)),
                           R=T[samp.b], a0=a0, W=re.w.simple, quiet = TRUE))
      if (class(sol.b)!="list" || sol.b$convergence!=0)
      {
        cat("bs rep=", b, "Error occurred in BBsolve. " , "\n")
        next
      }
      else
      {
        b.b <- rbind(b.b, sol.b$par)
        act.num <- act.num +1
      }
    }
    colnames(b.b) <- covNames
    se.bs <- apply(b.b, 2, sd)
    var.bs <- var(b.b)
    left.lim <- apply(b.b, 2, quantile, prob=(1-conf.int)/2, na.rm=TRUE)
    right.lim <- apply(b.b, 2, quantile, prob=1-(1-conf.int)/2, na.rm=TRUE)
    # the next row is a CI of [b.bar-z_alpha*SE,infty)  is for H0:b<=0, H1:b>0
    left.H1.b.gr0 <- apply(b.b, 2, quantile, prob=(1-conf.int), na.rm=TRUE)
    # the next row is a CI of ([)-infty, b.bar+z_alpha*SE]  is for H0:b>=0, H1:b<0
    right.H1.b.le0 <- apply(b.b, 2, quantile, prob=conf.int, na.rm=TRUE)
  }

  qn <- qnorm(1-(1-conf.int)/2)
  out$coef <- sol.a0$par
  out$convergence <- sol.a0$convergence

  if (bs)
  {
    stat <- data.frame(coef=sol.a0$par,
                       exp.coef=exp(sol.a0$par),
                       SE=se.bs,
                       CI.L=left.lim, CI.U=right.lim,
                       CI.L.H1.b.gr0=left.H1.b.gr0,
                       CI.U.H1.b.le0=right.H1.b.le0,
                       pvalue=2*pnorm(-abs(sol.a0$par)/se.bs),
                       pvalue.H1.b.gr0=pnorm(sol.a0$par/se.bs, lower.tail=FALSE),
                       pvalue.H1.b.le0=pnorm(sol.a0$par/se.bs, lower.tail=TRUE))
    out$bs$nbs.rep <- act.num
    out$bs$summary <- stat
    out$var <- var.bs
  }
  else {
    out$bs <- NULL
  }

  return(out)
}












