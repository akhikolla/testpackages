# Small variant of the pd.solve function of Adelchi Azzalini's 'mnormt' package
# that specifies the tolerance for recognising unsymmetric matrics as an
# aditional argument

pdwt.solve <- function(x, silent=FALSE, log.det=FALSE, onlylogdet=FALSE, symtol=.Machine$double.eps)
{
  if(is.null(x))  if(silent) return (NULL) else stop("pdwt.solve invoked with a NULL first argument")
  if(any(is.na(x))) if(silent) return (NULL) else stop("NA's in x") 
  if(length(x) == 1) x <- as.matrix(x)
  if(!is.matrix(x)) if(silent) return(NULL) else stop("x is not a matrix")
  maxsymdif <- max(abs(x - t(x))) 
  if(maxsymdif > symtol) if(silent) return (NULL) else stop("x appears to be not symmetric") 
  if (maxsymdif > 0.) x <- (x + t(x))/2
  u <- try(chol(x, pivot = FALSE), silent = silent)
  if(class(u)[1] == "try-error") if(silent) return(NULL) else stop("x appears to be not positive definite")
  if(log.det || onlylogdet) logdet <- 2*sum(log(diag(u)))
  if(onlylogdet) return(logdet)
  inv <- chol2inv(u)
  if(log.det) attr(inv, "log.det") <- logdet
	
  inv
}

# Safe variant of pdwt.solve that only accepts well-conditioned p.d. matrices
# (i.e. p.d. matrices with log-determinants above the ldetmax argument)

safepdwt.solve <- function(x, silent=TRUE, ldetmax=-10., onlylogdet=FALSE, symtol=.Machine$double.eps)
{
  inv <- pdwt.solve(x,silent=TRUE,log.det=TRUE,onlylogdet=onlylogdet,symtol=symtol)
  if (is.null(inv) || attr(inv,"log.det") < ldetmax)
  if (!silent) stop("Matrix appears no to be positive definite\n")
  else return(NULL)
  inv
}


