msnCP.dev <- function(param, y, grpind, Config, n=ifelse(is.matrix(y),nrow(y),1),
                     p=ifelse(is.matrix(y),ncol(y),length(y)), k=1,	
                     limlnk2=log(1e8), trace=FALSE, c2tol=1e-3, ldRtol=-500,
                     PenF=1e12, PenC=0., nopenalty=FALSE, Srpar=TRUE)
{  
  if (!all(is.finite(param))) return(.Machine$double.xmax)
  res <- .Call( "msnCP_dev", param, y, grpind, Config, n, p, k, limlnk2, trace,  
    c2tol, ldRtol, PenF, PenC, nopenalty, .Machine$double.eps, Srpar, PACKAGE = "MAINT.Data" )
  if (!is.finite(res))  return(.Machine$double.xmax/10)
  res
}

