volTriangulation <- function(vertices) { ##
  if (is.data.frame(vertices)) vertices <- as.matrix(vertices)
  ## probleme with repeated vertices occurs sometimes:
  vertices <- unique(vertices)
  ## ...otherwise it is possible that a vertex (from $vertices) is later selected, which is not in the $simplicesTable
  tc <- delaunayn(vertices,"Pp") ## triangulation by simplices
  # cf blckbox::volTriangulationWrapper for cases where delaunayn fails
  pmul <- cbind(-1,diag(ncol(vertices)))
  factorialdim <- factorial(ncol(vertices)) ## which is not always try if only a subspace is sampled
  vb <- apply(tc,1,function(v){
    simplex <- vertices[v,,drop=FALSE]
    # pmul %*% simplex is simplex[-1,]-simplex[1,] for each simplex
    # volume = abs(det())/ dim!
    vol <- abs(det(pmul %*% simplex))/factorialdim
    bary <- colMeans(simplex)
    c(vol,bary) ## volume, and barycenter \n\
  })
  resu <- list(vol=vb[1,], ## a vector
               bary=t(vb[-1,,drop=FALSE]), ## matrix
               vertices=vertices,simplicesTable=tc)
  class(resu) <- c("volTriangulation",class(resu))
  resu
}
#vT <- volTriangulation(mvH)


## 'simplices' are indices kept, or subtracted if negative
# originally named subset.volTriangulation but subset is a generic with a different usage
subsimplices.volTriangulation <- function(x,simplices,...) {
  vT <- x
  vT$vol <- vT$vol[simplices]
  vT$bary <- vT$bary[simplices,]
  vT$simplicesTable <- vT$simplicesTable[simplices,,drop=FALSE]
  # vT$simplicesTable refers to original point indices => $vertices is unchanged
  vT
}


## to sample uniformly wrt volume within a simplex, not uniformly wrt perimeter
# old complicated code with additional argument 'u' removed 04/2016
#
rsimplex <- function(simplex,
                     expand=NULL, ## *volume* expansion
                     bary=NULL ## used for expansion
) {
  d <- NROW(simplex)-1 ## 2 pour triangle
  if(NCOL(simplex)!= d) {stop("(!) From 'rsimplex': simplex should have one more row than columns.")}
  weights <- diff(sort(c(0, runif(n=d), 1))) ## Wilks 1962, p.238 in Rubin 1981
  # http://cs.stackexchange.com/questions/3227/uniform-sampling-from-a-simplex for some more discussion
  if ( ! is.null(expand)) {
    if (is.null(bary)) bary <- colMeans(simplex) ## OK pour bary de d-dim simplex avec densité uniforme
    simplex <- expand^(1/d) * sweep(simplex,2L,bary,`-`)
    simplex <- sweep(simplex,2L,bary,`+`)
  }
  ws <- sweep(simplex,1L,weights,`*`)
  return(colSums(ws))
}

## ideally we would have a rvolume S3 generic with .volTriangulation and .data.frame methods
rvolTriangulation <- function(n=1,volTriangulationObj,replace=TRUE,expand=NULL) {
  if (! inherits(volTriangulationObj,"volTriangulation")) stop("(!) From 'volTriangulation': input 'volTriangulationObj' is not a volTriangulation object.")
  simplexProbs <- volTriangulationObj$vol
  simplexProbs <- simplexProbs/sum(simplexProbs)
  whichSimplex <- sample(seq_len(length(simplexProbs)),n,prob=simplexProbs,replace=replace)
  vertices <- volTriangulationObj$vertices
  vI <- volTriangulationObj$simplicesTable
  resu <- sapply(whichSimplex,function(idx) {rsimplex(vertices[vI[idx,],,drop=FALSE],expand=expand)})
  if (ncol(vertices)==1L) {
    resu <- as.matrix(resu)
  } else resu <- t(resu)
  rownames(resu) <- NULL ## bc all rownames=parm when ncol=1
  colnames(resu) <- colnames(vertices,do.NULL=FALSE)
  return(resu)
} ## end def rvolTriangulation


# old.nextPoints <- function(n=1,optr,replace=TRUE) { ## random sampling of volume defined from previous fit
#   uP <- upperPoints(optr$predictions) ## indices
#   uP <- optr$predictions[uP,attr(optr$predictions,"fittedPars")]
#   uP <- rbind(uP,optr$par) ## not sure this is useful for volumetric sampling
#   erV <- elim.redundant.V(uP)
#   vT <- volTriangulation(erV)
#   rvolTriangulation(n,vT,replace=replace)
# }

## FR->FR same as in blackbox,; (F I X M E ?) see spaMM:::.locate_in_tv() for an optimized version for 2D only but several points
locatePointinvT <- function(point, ## numeric (not matrix or data frame: see use of 'point' below)
                            vT,fallback=TRUE) { ## in which simplex ? with fall back if 'numerically outside' vT (but quite distinct from minimal distance)
  pmul <- cbind(-1,diag(ncol(vT$vertices)))
  minw <- apply(vT$simplicesTable[vT$vol>0,,drop=FALSE],1,function(v) {
    simplex <- vT$vertices[v,,drop=FALSE]
    vM <- pmul %*% simplex # is simplex[-1,]-simplex[1,] for each simplex
    vWeights <- try(solve(t(vM),point-simplex[1,]),silent=TRUE) ## as.numeric(point) would be required if point were a (1-row) matrix
    # problem may occur if volume of simplex is nearly zero
    if (inherits(vWeights,"try-error")) {vWeights <- ginv(t(vM)) %*% (point-simplex[1,])}
    vWeights <- c(1-sum(vWeights),vWeights) ## weights for all vertices
    min(vWeights) ## if the point is within/marginally outside/clearly outside the vT, this will return a positive/small negative/large negative value
  })
  resu <- which(minw>0)
  if (length(resu)==0L && fallback) resu <- which.max(minw) ## if numerically outside...
  return(resu)
}
