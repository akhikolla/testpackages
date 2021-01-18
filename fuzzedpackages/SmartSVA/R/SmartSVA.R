#' @importFrom utils getFromNamespace
# Accessing 'mono' function from 'SVA' package
getFromNamespace("mono", "sva")

#' @importFrom "stats" "dnorm" "qnorm" "density" "smooth.spline" "predict"
# This is the edge.lfdr function from sva package without any modifications.
edge.lfdr <- function (p, trunc = TRUE, monotone = TRUE,
                       transf = c("probit","logit"), adj = 1.5, eps = 10^-8,
                       lambda = 0.8, ...)
{
  pi0 <- mean(p >= lambda)/(1 - lambda)
  pi0 <- min(pi0, 1)
  n = length(p)
  transf = match.arg(transf)
  if (transf == "probit") {
    p = pmax(p, eps)
    p = pmin(p, 1 - eps)
    x = qnorm(p)
    myd = density(x, adjust = adj)
    mys = smooth.spline(x = myd$x, y = myd$y)
    y = predict(mys, x)$y
    lfdr = pi0 * dnorm(x)/y
  }
  if (transf == "logit") {
    x = log((p + eps)/(1 - p + eps))
    myd = density(x, adjust = adj)
    mys = smooth.spline(x = myd$x, y = myd$y)
    y = predict(mys, x)$y
    dx = exp(x)/(1 + exp(x))^2
    lfdr = pi0 * dx/y
  }
  if (trunc) {
    lfdr[lfdr > 1] = 1
  }

  # Get 'mono' function from 'SVA' package.
  mono <- getFromNamespace("mono", ns="sva")

  if (monotone) {
    lfdr = lfdr[order(p)]
    lfdr = mono(lfdr)
    lfdr = lfdr[rank(p)]
  }
  return(lfdr)
}


#' @importFrom "stats" "cor" "pf"
#'
#' @import sva
#' @import isva
# Computing p-value from an F-test.
f.pval <- function (dat, orth11, orth01, y.norm, rss00, df00)  {

  n <- dim(dat)[2]

  df11 <- dim(orth11)[2]
  df01 <- dim(orth01)[2]

  prj11 <- dat %*% orth11
  prj01 <- dat %*% orth01

  rss11 <- y.norm - rowSums(prj11 * prj11)
  rss01 <- y.norm - rowSums(prj01 * prj01)

  fstats <- ((rss01 - rss11)/(df11 - df01))/(rss11/(n - df11))
  p1 <- 1 - pf(fstats, df1 = (df11 - df01), df2 = (n - df11))

  fstats <- ((rss00 - rss01)/(df01 - df00))/(rss01/(n - df01))
  p2 <- 1 - pf(fstats, df1 = (df01 - df00), df2 = (n - df01))

  return(list(p1=p1, p2=p2))
}


#' Performs a significantly more efficient surrogate variable analysis for EWAS.
#'
#' @param dat the measurement matrix (preferably M-values)
#' where rows are variables and columns are samples.
#'
#' @param mod the model matrix used to fit data.
#'
#' @param mod0 contains the nul model matrix.
#'
#' @param n.sv number of svs. The use of random matrix theory is recommended to estimate
#' n.sv. See the example for more details.
#'
#' @param B iteration number which is typically varies between 20 and 100.
#'
#' @param alpha determines the initial point for optimization
#' which affects the convergence rate.
#'
#' @param epsilon specifies the convergence threshold. the spearman
#'  correlation between posterior probabilities of consecutive iterations
#'  of the algorithm is compared to epsilon. Empirical evaluation on several
#'  dataset revealed epsilon=0.005 gives very reasonable results.
#'  However, we suggest epsilon=1e-3 as a conservative threshold.
#'
#' @param VERBOSE a logical variable. If TRUE, prints some details about
#'  iterative progress of the algorithm.
#'
#' @return Returs a list containing the surrogate variables and some meta data
#' about the convergence criterion.
#'
#' @examples
#' ## Methylation M values (CpG by Sample)
#' Y <- matrix(rnorm(20*1000), 1000, 20)
#' df <- data.frame(pred=gl(2, 10))
#'
#' ## Determine the number of SVs
#' Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
#' ## Add one to compensate potential loss of 1 degree of freedom
#' ##  in confounded scenarios
#' n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
#' mod <- model.matrix( ~ pred, df)
#' sv.obj <- smartsva(Y, mod, mod0=NULL, n.sv=n.sv)
#'
#' @export
smartsva <-  function(dat, mod, mod0 = NULL, n.sv, B = 100,
		alpha=0.25, epsilon=1e-3, VERBOSE = F) {
	if (is.null(mod0)) {
		mod0 <- mod[, 1]
	}
	
	qr.obj <- qr(mod)
	orth1 <- qr.Q(qr.obj)
	uu <- eigen(crossprod(dat - tcrossprod(dat %*% orth1, orth1)),
			 symmetric=TRUE)$vectors[, 1:n.sv, drop=F]
	
	# Precompute the quantites
	y.norm <- rowSums(dat * dat)
	mod00 <- cbind(mod0)
	orth00 <- qr.Q(qr(mod00))
	prj00 <- dat %*% orth00
	rss00 <- y.norm - rowSums(prj00 * prj00)
	df00 <- dim(orth00)[2]
	
	if (VERBOSE)
		cat(paste("Iteration (out of", B, "):\n"))
	
	i = 0
	rho = 0
	
	while (i < B && rho < 1 - epsilon) {
		i <- i + 1
		mod11 <- cbind(mod, uu)
		mod01 <- cbind(mod0, uu)
		
		orth11 <- qr.Q(qr(mod11))
		orth01<- qr.Q(qr(mod01))
		
		ptmp <- f.pval(dat, orth11, orth01, y.norm, rss00, df00)
		
		if (i == 1) {
			pprob.b <- (1 - edge.lfdr(ptmp[['p1']])^alpha)
		} else {
			pprob.b <- (1 - edge.lfdr(ptmp[['p1']]))
		}
		
		pprob.gam <- (1 - edge.lfdr(ptmp[['p2']]))
		pprob <- pprob.gam * (1 - pprob.b)
		
		uu <- eigen(crossprod(dat * pprob - rowMeans(dat * pprob)),
				symmetric=TRUE)$vectors[, 1:n.sv, drop=F]
		# Update spearman Rho.
		if (i > 1) {
			rho <- cor(x=pprob, y=p.prev, use="pairwise.complete.obs",
					method="spearman")
			p.prev <- pprob
		}else{
			p.prev <- pprob
		}
		if (VERBOSE)
			cat(paste(i, " ", rho, "\n"))
	}
	
	sv <- uu[, 1:n.sv, drop=F]
	retval <- list(sv = sv, n.sv = n.sv, pprob.gam = pprob.gam, pprob.b = pprob.b, rho = rho, iter = i)
	return(retval)
}

#' @importFrom "stats" "cor" "pf"
#'
#' @import sva
#' @import isva
# Computing p-value from an F-test (EigenCpp accelerated).
f.pval.cpp <- function (dat, orth11, orth01, y.norm, rss00, df00)  {
	
	n <- dim(dat)[2]
	
	df11 <- dim(orth11)[2]
	df01 <- dim(orth01)[2]
	
	prj11 <- prodCpp(dat, orth11)
	prj01 <- prodCpp(dat, orth01)
	
	rss11 <- y.norm - rowSums(prj11 * prj11)
	rss01 <- y.norm - rowSums(prj01 * prj01)
	
	fstats <- ((rss01 - rss11)/(df11 - df01))/(rss11/(n - df11))
	p1 <- 1 - pf(fstats, df1 = (df11 - df01), df2 = (n - df11))
	
	fstats <- ((rss00 - rss01)/(df01 - df00))/(rss01/(n - df01))
	p2 <- 1 - pf(fstats, df1 = (df01 - df00), df2 = (n - df01))
	
	return(list(p1=p1, p2=p2))
}

#' Performs a significantly more efficient surrogate variable analysis for EWAS (CppEigen accelaration).
#'
#' @param dat the measurement matrix (preferably M-values)
#' where rows are variables and columns are samples.
#'
#' @param mod the model matrix used to fit data.
#'
#' @param mod0 contains the nul model matrix.
#'
#' @param n.sv number of svs. The use of random matrix theory is recommended to estimate
#' n.sv. See the example for more details.
#'
#' @param B iteration number which is typically varies between 20 and 100.
#'
#' @param alpha determines the initial point for optimization
#' which affects the convergence rate.
#'
#' @param epsilon specifies the convergence threshold. the spearman
#'  correlation between posterior probabilities of consecutive iterations
#'  of the algorithm is compared to epsilon. Empirical evaluation on several
#'  dataset revealed epsilon=0.005 gives very reasonable results.
#'  However, we suggest epsilon=1e-3 as a conservative threshold.
#'
#' @param VERBOSE a logical variable. If TRUE, prints some details about
#'  iterative progress of the algorithm.
#'
#' @return Returs a list containing the surrogate variables and some meta data
#' about the convergence criterion.
#'
#' @examples
#' ## Methylation M values (CpG by Sample)
#' Y <- matrix(rnorm(20*1000), 1000, 20)
#' df <- data.frame(pred=gl(2, 10))
#'
#' ## Determine the number of SVs
#' Y.r <- t(resid(lm(t(Y) ~ pred, data=df)))
#' ## Add one to compensate potential loss of 1 degree of freedom
#' ##  in confounded scenarios
#' n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
#' mod <- model.matrix( ~ pred, df)
#' sv.obj <- smartsva.cpp(Y, mod, mod0=NULL, n.sv=n.sv)
#'
#' @export

smartsva.cpp <-  function(dat, mod, mod0 = NULL, n.sv, B = 100,
		alpha=0.25, epsilon=1e-3, VERBOSE = F) {
	if (is.null(mod0)) {
		mod0 <- mod[, 1]
	}
	
	qr.obj <- qr(mod)
	orth1 <- qr.Q(qr.obj)
	uu <- eigs_sym(crossprodCpp(dat - tcrossprodCpp(prodCpp(dat, orth1), orth1)),
			k=n.sv)$vectors[, 1:n.sv, drop=F]
	
	# Precompute the quantites
	y.norm <- rowSums(dat * dat)
	mod00 <- cbind(mod0)
	orth00 <- qr.Q(qr(mod00))
	prj00 <- prodCpp(dat, orth00)
	rss00 <- y.norm - rowSums(prj00 * prj00)
	df00 <- dim(orth00)[2]
	
	if (VERBOSE)
		cat(paste("Iteration (out of", B, "):\n"))
	
	i = 0
	rho = 0
	
	while (i < B && rho < 1 - epsilon) {
		i <- i + 1
		mod11 <- cbind(mod, uu)
		mod01 <- cbind(mod0, uu)
		
		orth11 <- qr.Q(qr(mod11))
		orth01<- qr.Q(qr(mod01))
		
		ptmp <- f.pval.cpp(dat, orth11, orth01, y.norm, rss00, df00)
		
		if (i == 1) {
			pprob.b <- (1 - edge.lfdr(ptmp[['p1']])^alpha)
		} else {
			pprob.b <- (1 - edge.lfdr(ptmp[['p1']]))
		}
		
		pprob.gam <- (1 - edge.lfdr(ptmp[['p2']]))
		pprob <- pprob.gam * (1 - pprob.b)
		
		uu <- eigs_sym(crossprodCpp(dat * pprob - rowMeans(dat * pprob)),
				k=n.sv)$vectors[, 1:n.sv, drop=F]
		# Update spearman Rho.
		if (i > 1) {
			rho <- cor(x=pprob, y=p.prev, use="pairwise.complete.obs",
					method="spearman")
			p.prev <- pprob
		}else{
			p.prev <- pprob
		}
		if (VERBOSE)
			cat(paste(i, " ", rho, "\n"))
	}
	
	sv <- uu[, 1:n.sv, drop=FALSE]
	retval <- list(sv = sv, n.sv = n.sv, pprob.gam = pprob.gam, pprob.b = pprob.b, rho = rho, iter = i)
	return(retval)
}



