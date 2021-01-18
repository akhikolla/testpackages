pvalues <- function(object)
{
	if(attr(object,'standardized'))
	 	return(2*pnorm(-abs(object$J)))
	else
		 stop("This function is only to be used with standardized results.")
}

