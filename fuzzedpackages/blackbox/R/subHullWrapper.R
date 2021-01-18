############################### 'equality' is a structure which each element is [a unique named value, compatible with unlist()]
# equality is a list(varname1=value1,varname2=value2,...)
# FR 05/2018: it's not clear whether I wanted to use geometry rather than rcdd, or the reverse.
#   subHullWrapper() itself uses rcdd; but calls redundant.addVeq() whose 'preferences' are unclear:
#   it uses rcdd if (include.Hrepr), otherwise depends only on geometry::convhulln() except when the latter fails 
#   in which case rcdd::redundant() is called. If we used redundant() consistently in redundant.addVeq(), 
#   the whole subHullWrapper() code would ot depend on geometry::. 
subHullWrapper <- function(vertices, equality, include.Hrepr=FALSE,precision="rational") {
  tmp <- vertices
  for (ii in seq_len(length(equality))) {
    value <- unlist(equality[ii]) ## named scalar
    tmp <- redundant.addVeq(vertices=tmp, value=value)
  }
  if ( ! is.null(tmp) && include.Hrepr ) {
    if (precision=="double") {
      Hrepr <- try(scdd(cbind(0, 1, tmp) , representation="V", roworder="maxcutoff")$output,silent=TRUE)
      if (inherits(Hrepr,"try-error")) {
        Hrepr <- NULL
      } else {
        verif <- apply(Hrepr[,-c(1,2)] %*% t(tmp) + Hrepr[,2],2L,range) ## matrix elements must be >= 0
        verif <- apply(verif,2L,range) ## first row must be ~0 if points are vertices 
        ##                                (second row must be positive but first test is sufficient) 
        if (max(abs(verif[1,]))>1e-8) Hrepr <- NULL ## and then next tests -> use rational 
      }
    }
    if (precision=="rational" || is.null(Hrepr)) { ## FR->FR bottleneck
      Hrepr <- q2d(scdd(d2q( cbind(0, 1, tmp) ), representation="V", roworder="maxcutoff")$output) 
    }
    if (any(abs(Hrepr)>1e16)) { ## Problems sometimes happen despite a rational computation
      ## patch deduced from a single case and clearly heuristic:
      Hrepr <- Hrepr[ (apply(abs(Hrepr),1,max)<1e16),,drop=FALSE] ## remove suspect constraints
    }
    colnames(Hrepr) <- c("eq", "b", colnames(tmp) ) ## Hrep should never be NULL here
  } else Hrepr <- NULL
  return(list(vertices=tmp, Hrepr=Hrepr))
}
####### Hrepr used in profileBySubHull -> generateInitpts -> rhull -> scdd.addHline; and in constroptim calls within profile computations.
