#####################################################################
####
#### Conversion between different parametrizations of mixed models
####
#####################################################################

#' @title Conversion between different parametrizations of mixed models
#'
#' @description Functions to convert between canonical parametrization
#'     (g,h,K), moment parametrization (p,m,S) and mixed
#'     parametrization (p,h,K).
#'
#' @name parm-conversion
#' 
#' @param parms Parameters of a mixed interaction model
#' @param type Output parameter type; either "ghk" or "pms".
#' @param SS List of moment parameters.
#' @return Parameters of a mixed interaction model.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @keywords utilities
#' 

#' @rdname parm-conversion
#' @export
parm_pms2ghk <- function(parms){
    ##cat("pms2ghkParms\n")
    c(.Call("C_pms2ghk", parms), parms[-(1:4)])
}

#' @rdname parm-conversion
#' @export
parm_ghk2pms <- function(parms){
    ##cat("ghk2pmsParms\n")
    c(.Call("C_ghk2pms", parms), parms[-(1:4)])
}

#' @rdname parm-conversion
#' @export
parm_pms2phk <- function(parms){
  ##parms <- unclass(parms)  
  switch(parms[['gentype']],
         "discrete"={
           res    <- list(p=parms[["p"]],h=NULL, K=NULL, gentype="discrete")
         },
         "mixed"={
           KK     <- solveSPD(parms[["Sigma"]])           
           hh     <- KK %*% parms[["mu"]] # h = Sigma.inv %*% mu
           res    <-list(p=parms[['p']], h=hh, K=KK, gentype="mixed")
         },
         "continuous"={
           KK     <- solveSPD(parms[["Sigma"]])
           hh     <- KK %*% parms[["mu"]] # h = Sigma.inv %*% mu           
           res    <- list(p=parms[['p']], h=hh, K=KK,gentype="continuous")
         })
  
  val <- c(res, parms[-(1:4)])
  ##class(val)<- c("phk","MIparms")
  val
}

#' @rdname parm-conversion
#' @export
parm_phk2ghk <- function(parms){
  ##parms <- unclass(parms)    
  switch(parms[['gentype']],
         "discrete"={
           res    <- list(g=log(parms[["p"]]),h=NULL, K=NULL, gentype="discrete")
         },
         "mixed"={
           KK     <- parms[["K"]]           
           Q      <- nrow(KK)
           detSig <- 1/det(KK)
           hh     <- parms[["h"]]
           mu     <- solveSPD(KK) %*% hh # K.inv %*% h
           quad   <- colSumsPrim(hh * mu)

           gg     <- log(parms[["p"]]) + (- log(detSig) - Q*log(2*pi) - quad) / 2
           res    <- list(g=gg, h=hh, K=KK, gentype="mixed")
         },
         "continuous"={
           KK     <- parms[["K"]]           
           Q      <- nrow(KK)
           detSig <- 1/det(KK)
           mu     <- solveSPD(KK) %*% hh # K.inv %*% h
           quad   <- colSumsPrim(hh * mu)

           gg     <- (- log(detSig) - Q*log(2*pi) - quad) / 2
           res    <-list(g=gg, h=hh, K=KK, gentype="continuous")
         })
  
  val <- c(res, parms[-(1:4)])
  ##class(val)<- c("ghk","MIparms")
  val
}

#' @rdname parm-conversion
#' @export
parm_phk2pms <- function(parms){
  ##parms <- unclass(parms)  
  switch(parms[['gentype']],
         "discrete"={
           res <- list(p=parms[['g']], mu=NULL, Sigma=NULL, gentype="discrete")
         },
         "mixed"={
           Sigma     <- solveSPD(parms[['K']])
           mu        <- Sigma %*% parms[['h']] # Kinv %*% h
           res       <- list(p=parms[['p']], mu=mu, Sigma=Sigma, gentype="mixed")  
         },
         "continuous"={
           Sigma     <- solveSPD(parms[['K']])
           mu        <- Sigma %*% parms[['h']]
           res       <- list(p=1, mu=mu, Sigma=Sigma, gentype="continuous")
         })

  val <- c(res, parms[-(1:4)])
  ##class(val)<- c("pms","MIparms")
  val
}

#' @rdname parm-conversion
#' @export
parm_ghk2phk <- function(parms){
  ##parms <- unclass(parms)  
  switch(parms[['gentype']],
         "discrete"={
           zzz <- parms[['g']]
           pp  <- exp(zzz - mean.default(zzz)) 
           res <- list(p=pp / sum(pp), mu=NULL, Sigma=NULL, gentype="discrete")
         },
         "mixed"={
           Sigma     <- solveSPD(parms[['K']])
           hh        <- parms[['h']]
           mu        <- Sigma %*% hh # Kinv %*% h
           g.quad    <- parms[['g']] + colSumsPrim(hh * mu)/2
           pp        <- exp(g.quad - mean.default(g.quad))           
           res       <- list(p=pp / sum(pp), h=hh, K=parms[['K']], gentype="mixed")  
         },
         "continuous"={
           Sigma     <- solveSPD(parms[['K']])
           mu        <- Sigma %*% parms[['h']]
           res       <- list(p=1, h=parms[['h']], K=parms[['K']], gentype="continuous")
         })

  val <- c(res, parms[-(1:4)])
  ##class(val)<- c("phk","MIparms")
  val
}

#' @rdname parm-conversion
#' @export
parm_CGstats2mmod <- function(parms, type="ghk"){
  type <- match.arg(type, c("ghk", "pms"))

  ans <- switch(type,
                "pms"={parm_moment2pms(parms)},
                "ghk"={parm_pms2ghk(parm_moment2pms(parms))}
                )                
  ans
}

#' @rdname parm-conversion
#' @export
parm_moment2pms <- function(SS){
  parms <- list(p=SS$n.obs/sum(SS$n.obs), mu=SS$center, Sigma=SS$cov, n.total=sum(SS$n.obs), gentype="mixed")
  ans   <- c(parms, SS[-(1:3)])
  ##class(ans) <- c("pms","MIparms")
  ans
}



### dot-functions below here





.pms2ghkParms <- function(parms){
  ##parms <- unclass(parms)  
  switch(parms[['gentype']],
         "discrete"={
           res    <- list(g=log(parms[["p"]]),h=NULL, K=NULL, gentype="discrete")
         },
         "mixed"={
           KK     <- solveSPD(parms[["Sigma"]])           
           Q      <- nrow(KK)
           #.logdetSig <- -log(det(KK))
           .logdetSig <- -c(determinant.matrix(KK)[['modulus']])

           mu     <- parms[["mu"]]
           hh     <- KK %*% mu # h = Sigma.inv %*% mu
           quad   <- colSumsPrim(hh * mu)
           gg     <- log(parms[["p"]]) + (- .logdetSig - Q*log(2*pi) - quad) / 2
           res    <- list(g=gg, h=hh, K=KK, gentype="mixed")
         },
         "continuous"={
           KK     <- solveSPD(parms[["Sigma"]])
           Q      <- nrow(KK)
           #detSig <- 1/det(KK) ##FIXME: brug det.matrix isf.
           detSig <- 1/c(determinant.matrix(KK, logarithm=FALSE)[['modulus']])
           mu     <- parms[["mu"]]           
           hh     <- KK %*% mu # h = Sigma.inv %*% mu
           quad   <- colSumsPrim(hh * mu)
           gg     <- (- log(detSig) - Q*log(2*pi) - quad) / 2
           res    <- list(g=gg, h=hh, K=KK,gentype="continuous")
         })
  
  val <- c(res, parms[-(1:4)])
  ##class(val)<- c("ghk","MIparms")
  return(val)
}

.ghk2pmsParms<-function(parms){

  ##parms <- unclass(parms)
  switch(parms[['gentype']],
         "discrete"={
           zzz <- parms[['g']]
           pp  <- exp(zzz-mean.default(zzz))
           res <- list(p=pp/sum(pp), mu=NULL, Sigma=NULL, gentype="discrete")
         },
         "mixed"={
           Sigma     <- solveSPD(parms[['K']])
           hh        <- parms[['h']]
           mu        <- Sigma %*% hh        # Kinv %*% h
           g.quad    <- parms[['g']] + colSumsPrim(hh * mu)/2
           pp        <- exp( g.quad - mean.default(g.quad) )           
           res       <- list(p=pp/sum(pp), mu=mu, Sigma=Sigma, gentype="mixed")  
         },
         "continuous"={
           Sigma     <- solveSPD(parms[['K']])
           mu        <- Sigma %*% parms[['h']]
           res       <- list(p=1, mu=mu, Sigma=Sigma, gentype="continuous")
         })

  val <- c(res, parms[-(1:4)])
  ##class(val)<- c("pms","MIparms")
  return(val)
}


### Normalizes ghK representation
###
.normalize.ghkParms <- function(parms){

  K.idx <- 3
  hh        <- parms[['h']]
  mu        <- solveSPD(parms[[K.idx]]) %*% hh
  logdetK   <- .logdet(parms[[K.idx]])
  Q         <- nrow(parms[[K.idx]])

  quad   <- colSumsPrim(hh * mu)
  zzz    <- parms[['g']] + quad / 2
  ppp    <- exp( zzz - mean.default(zzz))
  normconst <- sum(ppp)
  pppn      <- ppp / normconst

  g.new <- log(pppn) + (logdetK - Q*log(2*pi) - quad)/2

  parms[['g']] <- g.new
  parms
}

## sandbox
.normalize.ghkParms <- function(parms){

  K.idx <- 3
  hh        <- parms[['h']]
  mu        <- solveSPD(parms[[K.idx]]) %*% hh
  logdetK   <- .logdet(parms[[K.idx]])
  Q         <- nrow(parms[[K.idx]])

  #cat("logdetK=", logdetK, "\n")
  #print(list(dimh=dim(hh), dimK=dim(parms[[K.idx]])))
  
  #print(hh * mu)
  quad   <- colSumsPrim(hh * mu)
  #print(list(dimh=dim(hh), dimK=dim(parms[[K.idx]]), dimq=dim(quad)))

  zzz    <- parms[['g']] + quad / 2
  ##print(zzz)
  ppp    <- exp( zzz - mean.default(zzz))
  ##print(ppp)
  normconst <- sum(ppp)
  pppn      <- ppp / normconst

  g.new <- log(pppn) + (logdetK - Q*log(2*pi) - quad)/2

  #cat("g.new:\n"); print(g.new)
  parms[['g']] <- g.new
  parms
}





