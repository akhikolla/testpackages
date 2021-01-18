#' Profile likelihood based confidence intervals for GPD
#' 
#' Calculates profile likelilhood based confidence intervals for a given fitted GPD model -- this is only implemented for two parameter GPD with no covariates in the model.
#' @param z a fitted \code{evmOpt} object
#' @param m return period : units are number of observations
#' @param xmax point estimate of the return level, this is used to bracket the roots of the equation used to calculate the ends of the profile likelihood based confidence interval.  The value need not be exact.
#' @param xlow value lower than the lower end of the confidence interval, for bracketing in root finding
#' @param conf confidence level, defaults to 0.95
#' @param nint used for plotting if required, number of points at which to calculate the profile likelihood for plotting, defaults to 50
#' @param PlotIt logical, whether or not to plot the profile likelihood, defaults to \code{FALSE}
#' @param mult used to calculate the starting point for the root finding for solving to find the upper end of the confidence interval.  The starting point is \code{mult*xmax} minus the lower end point. If this starting point is beyond the estimated upper endpoint of the fitted distribution then this can cause an error, and the value of \code{mult} should be reduced
#' @param priorParameters optional, value of prior/penalty parameters used for penalised likelihood estimation, default to NULL
#' @return Numeric vector of length two, with lower and upper ends of the estiamted confidence intervals respectively.
#' @export
gpd.prof <- 
	function (z, m, xmax, xlow, conf = 0.95, nint = 50, PlotIt=FALSE,mult=2,priorParameters=NULL) 
	{
		if(class(z) != "evmOpt"){
			stop("Please pass an object of class evmOpt in gpd.prof")
		}
		if(length(z$par) >2 ){
			stop("gpd.prof only implemented for 2 parameter GPD models, no covariates in the model.")
		}
		xdat <- z$data[[1]]
		u <- z$threshold
		la <- z$rate
		
		sol <- z$par[2]
		prior <- .make.quadratic.penalty(priorParameters)
		gpd.plik <- function(a,xp) {
			if (m != Inf) 
				sc <- (a * (xp - u))/((m * la)^a - 1)
			else sc <- (u - xp)/a
			if (abs(a) < 10^(-4)) 
				l <- length(xdat) * log(sc) + sum(xdat - u)/sc
			else {
				y <- (xdat - u)/sc
				y <- 1 + a * y
				if (any(y <= 0) || sc <= 0) 
					l <- 10^6
				else l <- length(xdat) * log(sc) + sum(log(y)) * 
					(1/a + 1) + prior(c(log(sc),a))
			}
			l
		}
		plikInt <- function(r){
			f <- optim(z$par[2], gpd.plik, method = "BFGS",xp=r)$value 
			f + z$loglik - dev
		}

		dev <- 0.5 * stats::qchisq(conf, 1)
		Low <- uniroot(plikInt,c(xlow,xmax))$root
		Upp <- .secant.method(plikInt,xmax,mult*xmax-Low)
		
		if(PlotIt){
			v <- numeric(nint)
			x <- seq(Low-0.1, Upp+0.1, length = nint)
			
			for (i in 1:nint) {
				opt <- optim(sol, gpd.plik, method = "BFGS",xp=x[i])
				sol <- opt$par
				v[i] <- opt$value
			}
			plot(x, -v, type = "l", xlab = "Return Level", ylab = "Profile Log-likelihood",xlim=c(Low-0.1,Upp+0.1))
			ma <- z$loglik
			abline(h = ma)
			abline(h = ma - dev)
			abline(v=Low)
			abline(v=Upp)
		}
		result <- c(Low,Upp)
		names(result) <-  c(paste("Lower",conf,sep=""),
							paste("Upper",conf,sep=""))					
		result
	}

.secant.method <- function(f, x0, x1, tol = 1e-9, n = 1000) {
	for (i in 1:n) {
		x2 <- x1 - f(x1) / ((f(x1) - f(x0)) / (x1 - x0)) # Calculate the new x value
		if (abs(x2 - x1) < tol) { # If the difference between the new value and the previous value is small enough, end iteration and output root.
			return(x2)
		}
		# If the root was not determined in the previous iteration, update the values and proceed to the next iteration.
		x0 <- x1
		x1 <- x2 
	}
}

.make.quadratic.penalty <- function(priorParameters=NULL) {
	
	if(is.null(priorParameters)){
		f <- function(param) 0
	} else {
		centre <- priorParameters[[1]]
		cov <- priorParameters[[2]]
	
		factors <- svd(as.matrix(cov))
		if (min(factors$d) == 0) {
			stop("Singular covariance matrix: quadratic penalty impossible")
			}
	
		stacked <- rbind(t(factors$u), t(factors$v))
		n.prior <- length(centre)
		
		f <- function(param) {
			delta <- param - centre
			prod <- as.numeric(stacked %*% delta)
			first <- prod[1:n.prior]
			last <- prod[(1 + n.prior):(2*n.prior)]
			as.numeric(first %*% (last / factors$d))
		}
	}
	f
}
