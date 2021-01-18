#' Provide full marginal reference distribution for for maringal transformation
#'
#' This gives the option of providing a set of reference marginal distributions to use for marginal transformation if the data's own marginal distribution is not appropriate (for instance if only data for which one variable is large is available, the marginal distributions of the other variables will not be represented by the available data).  In such situations, the user can supply the full marginal information of the non-thresholded variables which are necessary to transform these variables correctly from the original margins to Gumbel/Laplace for estimation of dependence model parameters.
#' 
#' @param x output from migpd fit to the original data which does not represent at least one marginal distribution
#' @param r output from migpd fit to the reference data which does represent the correct marginal distribution of the variable with incomplete representation in \code{x}
#' @param whichNoChange Margins which are not to use the supplied reference distribution \code{r} have numeric indices, giving column numbers in original dataframe, listed in \code{whichNoChange}.
#' @return An object of class "migpd".
#' @export
makeReferenceMarginalDistribution <- function(x,r,whichNoChange=NULL){ 
	# x,r both output from migpd. Margins which are not to use reference distribution r have numeric indices listed in whichNoChange
	r$transData <- r$models
	x$transData <- x$models
	for(i in 1:length(x$transData)){
		x$transData[[i]] <- x$data[,i]
		r$transData[[i]] <- r$data[,i]
	}
	for(i in whichNoChange){
		r$models[[i]] <- x$models[[i]]
		r$transData[[i]] <- x$transData[[i]]
		r$mth[i] <- x$mth[i]
		r$mqu[i] <- x$mqu[i]
		r$priorParameters[[i]] <- x$priorParameters[[i]]
	}
	r
}


mexTransform <-
    function(x, margins, r=NULL, method = "mixture", divisor = "n+1", na.rm = TRUE){ # if required, r output from makeReferenceMarginalDistribution
        
        if (!is.element(method, c("mixture", "empirical")))
            stop("method should be either 'mixture' or 'empirical'")
        if (!is.element(divisor, c("n", "n+1")))
            stop("divisor can be 'n' or 'n+1'")
        
        if (is.null(r)){
            r <- x
            r$transData <- lapply(1:dim(x$data)[2],function(i)x$data[,i])
        }
        
        transFun <- function(i, x, r, mod, th, divisor, method){
            x <- x[, i]
            r <- r[[i]]
            mod <- mod[[i]]
            th <- th[i]
            
            if (divisor == "n") divisor <- length(r)
            else if (divisor == "n+1") divisor <- length(r) + 1
            
            ox <- order(x)
            
            r <- sort(r)
            run <- rle(r)
            p <- cumsum(run$lengths) / divisor
            p <- rep(p, run$lengths) # this calculated from r
            
            Femp <- p[sapply(x,function(y) which.min(abs(r-y)))]
            if (method == "mixture"){
                sigma <- exp(mod$coefficients[1])
                xi <- mod$coefficients[2]
                
                Para <- (1 + xi * (x - th) / sigma) ^ (-1 / xi) # this calculated from model fitted to r but data values are x
                Para <- 1 - mean(r > th) * Para
                res <- ifelse(x <= th, Femp, Para)
            }
            else res <- Femp
            
            res[ox] <- sort(res)
            res
        } # Close transfun
        
        res <- sapply(1:ncol(x$data), transFun,
                      x = x$data, r = r$transData, mod = r$models, th = r$mth,
                      divisor = divisor, method=method)
        
        dimnames(res) <- list(NULL, names(r$models))
        
        x$transformed <- margins$p2q(res)
        
        invisible(x)
    }
