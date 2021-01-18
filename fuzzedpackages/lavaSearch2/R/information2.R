### information2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 19 2018 (14:17) 
## Version: 
## Last-Updated: jul 31 2020 (10:44) 
##           By: Brice Ozenne
##     Update #: 276
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - information2
#' @title  Extract The full Information Matrix
#' @description  Extract the full information matrix from a Gaussian linear model.
#' @name information2
#'
#' @param object a linear model or a latent variable model
#' @param ... arguments to be passed to \code{vcov2}.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @return A matrix.
#' 
#' @examples
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- Sequence.lvm(0)
#' set.seed(10)
#' d <- lava::sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' info.tempo <- vcov2(e.lm, bias.correct = TRUE)
#' info.tempo[names(coef(e.lm)),names(coef(e.lm))] - vcov(e.lm)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' vcov.tempo <- vcov2(e.lvm, bias.correct = FALSE)
#' round(vcov.tempo %*% information(e.lvm), 5)
#'
#' @concept small sample inference
#' @export
`information2` <-
  function(object, ...) UseMethod("information2")

## * information2.lm
#' @rdname information2
#' @export
information2.lm <- function(object, ...){
    return(solve(vcov2(object, ...)))
}

## * information2.gls
#' @rdname information2
#' @export
information2.gls <- information2.lm

## * information2.lme
#' @rdname information2
#' @export
information2.lme <- information2.lm

## * information2.lvmfit
#' @rdname information2
#' @export
information2.lvmfit <- information2.lm

## * information2.lm2
#' @rdname information2
#' @export
information2.lm2 <- function(object, ...){
    return(solve(vcov2(object, ...)))
}

## * information2.gls2
#' @rdname information2
#' @export
information2.gls2 <- information2.lm2

## * information2.lme2
#' @rdname information2
#' @export
information2.lme2 <- information2.lm2

## * information2.lvmfit
#' @rdname information2
#' @export
information2.lvmfit2 <- information2.lm2

## * .information2
#' @title Compute the Expected Information Matrix From the Conditional Moments
#' @description Compute the expected information matrix from the conditional moments.
#' @name information2-internal
#' 
#' @details \code{.information2} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.information2 <- function(dmu, dOmega,
                          Omega, n.corrected,
                          index.Omega, leverage, n.cluster,
                          grid.meanparam, n.grid.meanparam,
                          grid.varparam, n.grid.varparam,
                          name.param, n.param){

### ** Prepare
    test.global <- is.null(index.Omega)
    if(!test.global){
        OmegaM1 <- lapply(1:n.cluster, function(iC){
            return(chol2inv(chol(Omega[index.Omega[[iC]],index.Omega[[iC]]])))
        })
    }else{
        OmegaM1 <- chol2inv(chol(Omega))
    }
    
    Info <- matrix(0, nrow = n.param, ncol = n.param,
                   dimnames = list(name.param,name.param))
    if(length(dmu)>0){
        index.meanparam <- 1:n.grid.meanparam
    }else{
        index.meanparam <- NULL
    }
    if(length(dOmega)>0){
        index.varparam <- 1:n.grid.varparam
    }else{
        index.varparam <- NULL
    } 

### ** Global    
    if(test.global){
        ## *** Information relative to the mean parameters
        for(iG in index.meanparam){ # iG <- 1
            iP1 <- grid.meanparam[iG,1]
            iP2 <- grid.meanparam[iG,2]

            Info[iP1,iP2] <- Info[iP1,iP2] + sum(dmu[[iP1]] %*% OmegaM1 * dmu[[iP2]])
        }

        ## *** Information realtive to the variance parameters
        for(iG in index.varparam){ # iG <- 1
            iP1 <- grid.varparam[iG,1]
            iP2 <- grid.varparam[iG,2]

            iDiag <- diag(OmegaM1 %*% dOmega[[iP1]] %*% OmegaM1 %*% dOmega[[iP2]])
            Info[iP1,iP2] <- Info[iP1,iP2] + 1/2*sum(iDiag*n.corrected)
        }
    }

### ** Individual specific (missing data)
    if(!test.global){
        ## *** Information relative to the mean parameters
        for(iC in 1:n.cluster){ # iC <- 1
            iIndex <- index.Omega[[iC]]

            for(iG in index.meanparam){ # iG <- 1
                iP1 <- grid.meanparam[iG,1]
                iP2 <- grid.meanparam[iG,2]

                Info[iP1,iP2] <- Info[iP1,iP2] + sum(dmu[[iP1]][iC,iIndex,drop=FALSE] %*% OmegaM1[[iC]] * dmu[[iP2]][iC,iIndex,drop=FALSE])            
            }

            ## *** Information relative to the variance parameters
            for(iG in index.varparam){ # iG <- 1
                iP1 <- grid.varparam[iG,1]
                iP2 <- grid.varparam[iG,2]
                iDiag <- diag(OmegaM1[[iC]] %*% dOmega[[iP1]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[[iC]] %*% dOmega[[iP2]][iIndex,iIndex,drop=FALSE])
                Info[iP1,iP2] <- Info[iP1,iP2] + 1/2 * sum(iDiag * (1 - leverage[iC,iIndex]))
            }
        }        
    }

    ## ** Make Info a symmetric matrix
    Info <- symmetrize(Info, update.upper = NULL)
        
    ## ** export
    return(Info)
}

## * .hessian2
#' @title Compute the Hessian Matrix From the Conditional Moments
#' @description Compute the Hessian matrix from the conditional moments.
#' @name information2-internal
#' 
#' @details \code{.hessian} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.hessian2 <- function(dmu, d2mu, dOmega, d2Omega, 
                      Omega, n.corrected,
                      index.Omega, leverage, n.cluster,
                      grid.meanparam, n.grid.meanparam,
                      grid.varparam, n.grid.varparam,
                      name.param, n.param, residuals){

    ## ** Prepare
    test.global <- is.null(index.Omega)
    if(!test.global){
        OmegaM1 <- lapply(1:n.cluster, function(iC){
            return(chol2inv(chol(Omega[index.Omega[[iC]],index.Omega[[iC]]])))
        })
    }else{
        OmegaM1 <- chol2inv(chol(Omega))
    }
    
    hessian <- array(0, dim = c(n.param, n.param, n.cluster),
                     dimnames = list(name.param,name.param,NULL))
    if(length(dmu)>0){
        index.meanparam <- 1:n.grid.meanparam
    }else{
        index.meanparam <- NULL
    }
    if(length(dOmega)>0){
        index.varparam <- 1:n.grid.varparam
    }else{
        index.varparam <- NULL
    }

    if((n.grid.meanparam>0) && (n.grid.varparam>0)){
        name.meanparam <- unique(unlist(grid.meanparam[,c("Var1","Var2")]))
        name.varparam <- unique(unlist(grid.varparam[,c("Var1","Var2")]))

        grid.hybridparam <- .combination(name.meanparam,name.varparam)
        n.hybridparam <- NROW(grid.hybridparam)
        index.hybridparam <- 1:n.hybridparam
    }else{
        grid.hybridparam <- NULL
        n.hybridparam <- 0
        index.hybridparam <- NULL
    }

    ## ** Global    
    if(test.global){
        ## *** second derivative relative to the mean parameters
        for(iG in index.meanparam){ # iG <- 1
            iP1 <- grid.meanparam[iG,1]
            iP2 <- grid.meanparam[iG,2]

            if(grid.meanparam[iG,"deriv12"]){
                term1 <- rowSums((d2mu[[iP1]][[iP2]] %*% OmegaM1) * residuals)
            }else if(grid.meanparam[iG,"deriv21"]){
                term1 <- rowSums((d2mu[[iP2]][[iP1]] %*% OmegaM1) * residuals)
            }else{
                term1 <- 0
            }
            term2 <- -rowSums((dmu[[iP1]] %*% OmegaM1) * dmu[[iP2]])
            hessian[iP1,iP2,] <- hessian[iP1,iP2,] + term1 + term2
            hessian[iP2,iP1,] <- hessian[iP1,iP2,]
        }

        ## *** second derivative relative to the variance parameters
        for(iG in index.varparam){ # iG <- 1
            iP1 <- grid.varparam[iG,1]
            iP2 <- grid.varparam[iG,2]

            term1a <- - diag(OmegaM1 %*% dOmega[[iP1]] %*% OmegaM1 %*% dOmega[[iP2]])
            term2 <- - rowSums((residuals %*% OmegaM1 %*% dOmega[[iP2]] %*% OmegaM1 %*% dOmega[[iP1]] %*% OmegaM1) * residuals)
            if(grid.varparam[iG,"deriv12"]){
                term1b <- diag(OmegaM1 %*% d2Omega[[iP1]][[iP2]])
                term3 <- 1/2 * rowSums((residuals %*% OmegaM1 %*% d2Omega[[iP1]][[iP2]] %*% OmegaM1) * residuals)
            }else if(grid.varparam[iG,"deriv21"]){
                term1b <- diag(OmegaM1 %*% d2Omega[[iP2]][[iP1]])
                term3 <- 1/2 * rowSums((residuals %*% OmegaM1 %*% d2Omega[[iP2]][[iP1]] %*% OmegaM1) * residuals)
            }else{
                term1b <- 0
                term3 <- 0
            }
            hessian[iP1,iP2,] <- hessian[iP1,iP2,] - 1/2 * rowSums( sweep(1-leverage, FUN = "*", STATS = term1a + term1b, MARGIN = 2) ) + term2 + term3
            hessian[iP2,iP1,] <- hessian[iP1,iP2,]
        }
        ## *** second derivative relative to the mean and variance parameters
        for(iG in index.hybridparam){ # iG <- 1
            iP1 <- grid.hybridparam[iG,1]
            iP2 <- grid.hybridparam[iG,2]

            term1 <- - rowSums((dmu[[iP1]] %*% OmegaM1 %*% dOmega[[iP2]] %*% OmegaM1) * residuals)
            if(!is.null(dmu[[iP2]]) && !is.null(dOmega[[iP1]])){
                term2 <- - rowSums((dmu[[iP2]] %*% OmegaM1 %*% dOmega[[iP1]] %*% OmegaM1) * residuals)
            }else{
                term2 <- 0
            }
            
            hessian[iP1,iP2,] <- hessian[iP1,iP2,] + term1 + term2
            hessian[iP2,iP1,] <- hessian[iP1,iP2,]
        }
    }

    ## ** Individual specific (missing data)
    if(!test.global){

        ## *** Information relative to the mean parameters
        for(iC in 1:n.cluster){ # iC <- 1
            iIndex <- index.Omega[[iC]]
            
            ## *** second derivative relative to the mean parameters
            for(iG in index.meanparam){ # iG <- 1
                iP1 <- grid.meanparam[iG,1]
                iP2 <- grid.meanparam[iG,2]

                if(grid.meanparam[iG,"deriv12"]){
                    term1 <- sum((d2mu[[iP1]][[iP2]][iC,iIndex,drop=FALSE] %*% OmegaM1[[iC]]) * residuals[iC,iIndex,drop=FALSE])
                }else if(grid.meanparam[iG,"deriv21"]){
                    term1 <- sum((d2mu[[iP2]][[iP1]][iC,iIndex,drop=FALSE] %*% OmegaM1[[iC]]) * residuals[iC,iIndex,drop=FALSE])
                }else{
                    term1 <- 0
                }
                term2 <- -sum((dmu[[iP1]][iC,iIndex,drop=FALSE] %*% OmegaM1[[iC]]) * dmu[[iP2]][iC,iIndex,drop=FALSE])
                hessian[iP1,iP2,iC] <- hessian[iP1,iP2,iC] + term1 + term2
                hessian[iP2,iP1,iC] <- hessian[iP1,iP2,iC]
            }

            ## *** second derivative relative to the variance parameters
            for(iG in index.varparam){ # iG <- 1
                iP1 <- grid.varparam[iG,1]
                iP2 <- grid.varparam[iG,2]

                term1a <- - diag(OmegaM1[[iC]] %*% dOmega[[iP1]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[[iC]] %*% dOmega[[iP2]][iIndex,iIndex,drop=FALSE])
                term2 <- - sum((residuals[iC,iIndex,drop=FALSE] %*% OmegaM1[[iC]] %*% dOmega[[iP2]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[[iC]] %*% dOmega[[iP1]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[[iC]]) * residuals[iC,iIndex,drop=FALSE])
                if(grid.varparam[iG,"deriv12"]){
                    term1b <- diag(OmegaM1[[iC]] %*% d2Omega[[iP1]][[iP2]][iIndex,iIndex,drop=FALSE])
                    term3 <- 1/2 * sum((residuals[iC,iIndex,drop=FALSE] %*% OmegaM1[[iC]] %*% d2Omega[[iP1]][[iP2]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[[iC]]) * residuals[iC,iIndex,drop=FALSE])
                }else if(grid.varparam[iG,"deriv21"]){
                    term1b <- diag(OmegaM1[[iC]] %*% d2Omega[[iP2]][[iP1]][iIndex,iIndex,drop=FALSE])
                    term3 <- 1/2 * sum((residuals[iC,iIndex,drop=FALSE] %*% OmegaM1[[iC]] %*% d2Omega[[iP2]][[iP1]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[[iC]]) * residuals[iC,iIndex,drop=FALSE])
                }else{
                    term1b <- 0
                    term3 <- 0
                }
                hessian[iP1,iP2,iC] <- hessian[iP1,iP2,iC] - 1/2 * sum( (1-leverage[iC,iIndex,drop=FALSE]) * (term1a + term1b) ) + term2 + term3
                hessian[iP2,iP1,iC] <- hessian[iP1,iP2,iC]
            }

            ## *** second derivative relative to the mean and variance parameters
            for(iG in index.hybridparam){ # iG <- 1
                iP1 <- grid.hybridparam[iG,1]
                iP2 <- grid.hybridparam[iG,2]

                term1 <- - sum((dmu[[iP1]][iC,iIndex,drop=FALSE] %*% OmegaM1[[iC]] %*% dOmega[[iP2]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[[iC]]) * residuals[iC,iIndex,drop=FALSE])
                if(!is.null(dmu[[iP2]]) && !is.null(dOmega[[iP1]])){
                    term2 <- - sum((dmu[[iP2]][iC,iIndex,drop=FALSE] %*% OmegaM1[[iC]] %*% dOmega[[iP1]][iIndex,iIndex,drop=FALSE] %*% OmegaM1[[iC]]) * residuals[iC,iIndex,drop=FALSE])
                }else{
                    term2 <- 0
                }

                hessian[iP1,iP2,iC] <- hessian[iP1,iP2,iC] + term1 + term2
                hessian[iP2,iP1,iC] <- hessian[iP1,iP2,iC]
            }
        }        
    }


    ## ** export
    return(hessian)
}

## * .dInformation2
#' @title Compute the First Derivative of the Expected Information Matrix
#' @description Compute the first derivative of the expected information matrix.
#' @name dInformation2-internal
#' 
#' @details \code{.dInformation2} will perform the computation individually when the
#' argument \code{index.Omega} is not null.
#' 
#' @keywords internal
.dInformation2 <- function(dmu, d2mu, dOmega, d2Omega,
                           Omega, OmegaM1, n.corrected,
                           index.Omega, leverage, n.cluster,
                           name.param, name.3deriv){
    n.param <- length(name.param)
    index.deriv <- match(name.3deriv, name.param)

### ** prepare
    test.global <- is.null(index.Omega)
    if(!test.global){
        M.template <- Omega
        M.template[] <- 0
    }
    
    dInfo <-  array(0,
                    dim = c(n.param, n.param, length(name.3deriv)),
                    dimnames = list(name.param, name.param, name.3deriv))
    
### ** loop
    for(iDeriv in index.deriv){ # iDeriv <- 4
        for(iP1 in 1:n.param){ # iP1 <- 1
            for(iP2 in iP1:n.param){ # iP2 <- 1
                
                iNameD <- name.param[iDeriv]
                iName1 <- name.param[iP1]
                iName2 <- name.param[iP2]
                 ## cat(iNameD," ",iName1,"",iName2,"\n")
                
                ## *** identify relevant terms
                test.Omega1 <- !is.null(dOmega[[iNameD]]) && !is.null(dOmega[[iName1]]) && !is.null(dOmega[[iName2]])
                test.Omega2a <- !is.null(d2Omega[[iNameD]][[iName1]]) && !is.null(dOmega[[iName2]])
                test.Omega2b <- !is.null(d2Omega[[iName1]][[iNameD]]) && !is.null(dOmega[[iName2]])
                test.Omega3a <- !is.null(d2Omega[[iNameD]][[iName2]]) && !is.null(dOmega[[iName1]])
                test.Omega3b <- !is.null(d2Omega[[iName2]][[iNameD]]) && !is.null(dOmega[[iName1]])
                
                test.mu1a <- !is.null(d2mu[[iNameD]][[iName1]]) && !is.null(dmu[[iName2]])
                test.mu1b <- !is.null(d2mu[[iName1]][[iNameD]]) && !is.null(dmu[[iName2]])
                test.mu2a <- !is.null(d2mu[[iNameD]][[iName2]]) && !is.null(dmu[[iName1]])
                test.mu2b <- !is.null(d2mu[[iName2]][[iNameD]]) && !is.null(dmu[[iName1]])
                test.mu3 <- !is.null(dOmega[[iNameD]]) && !is.null(dmu[[iName1]]) && !is.null(dmu[[iName2]])

                if((test.Omega1 + test.Omega2a + test.Omega2b + test.Omega3a + test.Omega3b + test.mu1a + test.mu1b + test.mu2a + test.mu2b + test.mu3) == 0){
                    next
                }

                ## *** extract quantities for computations 
                if(test.mu1a){
                    d2mu.D1 <- d2mu[[iNameD]][[iName1]]
                }else if(test.mu1b){
                    d2mu.D1 <- d2mu[[iName1]][[iNameD]]
                }
                if(test.mu2a){
                    d2mu.D2 <- d2mu[[iNameD]][[iName2]]
                }else if(test.mu2b){
                    d2mu.D2 <- d2mu[[iName2]][[iNameD]]
                }
                if(test.Omega2a){
                    d2Omega.D1 <- d2Omega[[iNameD]][[iName1]]
                }else if(test.Omega2b){
                    d2Omega.D1 <- d2Omega[[iName1]][[iNameD]]
                }
                if(test.Omega3a){
                    d2Omega.D2 <- d2Omega[[iNameD]][[iName2]]
                }else{
                    d2Omega.D2 <- d2Omega[[iName2]][[iNameD]]
                }
                
                if(test.global){
                    ## *** Global: extract quantities for computations
                    if(!is.null(dOmega[[iNameD]])){
                        OmegaM1.dOmega.D <- OmegaM1 %*% dOmega[[iNameD]]
                    }
                    if(!is.null(dOmega[[iName1]])){
                        OmegaM1.dOmega.1 <- OmegaM1 %*% dOmega[[iName1]]
                    }
                    if(!is.null(dOmega[[iName2]])){
                        OmegaM1.dOmega.2 <- OmegaM1 %*% dOmega[[iName2]]
                    }                    

                    ## *** Global: compute
                    if(test.Omega1){
                        iDiag1 <- diag(OmegaM1.dOmega.D %*% OmegaM1.dOmega.1 %*% OmegaM1.dOmega.2)
                        iDiag2 <- diag(OmegaM1.dOmega.1 %*% OmegaM1.dOmega.D %*% OmegaM1.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum(iDiag1 * n.corrected + iDiag2 * n.corrected)
                    }

                    if(test.Omega2a || test.Omega2b){
                        iDiag <- diag(OmegaM1 %*% d2Omega.D1 %*% OmegaM1.dOmega.2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * n.corrected)
                    }

                    if(test.Omega3a || test.Omega3b){
                        iDiag <- diag(OmegaM1.dOmega.1 %*% OmegaM1 %*% d2Omega.D2)
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * n.corrected)
                    }

                    if(test.mu1a || test.mu1b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(d2mu.D1 %*% OmegaM1 * dmu[[iName2]])
                    }

                    if(test.mu2a || test.mu2b){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(dmu[[iName1]] %*% OmegaM1 * d2mu.D2)
                    }
                    if(test.mu3){
                        dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - sum(dmu[[iName1]] %*% OmegaM1.dOmega.D %*% OmegaM1 * dmu[[iName2]])
                    }
                    
                }else{
                    for(iC in 1:n.cluster){ # iC <- 1
                        iIndex <- index.Omega[[iC]]
                        
                        if(!is.null(dOmega[[iNameD]])){
                            OmegaM1.dOmega.D <- OmegaM1[[iC]] %*% dOmega[[iNameD]][iIndex,iIndex]
                        }
                        if(!is.null(dOmega[[iName1]])){
                            OmegaM1.dOmega.1 <- OmegaM1[[iC]] %*% dOmega[[iName1]][iIndex,iIndex]
                        }
                        if(!is.null(dOmega[[iName2]])){
                            OmegaM1.dOmega.2 <- OmegaM1[[iC]] %*% dOmega[[iName2]][iIndex,iIndex]
                        }

                        if(test.Omega1){
                            iDiag1 <- diag(OmegaM1.dOmega.D %*% OmegaM1.dOmega.1 %*% OmegaM1.dOmega.2)
                            iDiag2 <- diag(OmegaM1.dOmega.1 %*% OmegaM1.dOmega.D %*% OmegaM1.dOmega.2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - 1/2 * sum((iDiag1+iDiag2) * (1 - leverage[iC,iIndex]))
                        }
                        if(test.Omega2a || test.Omega2b){
                            iDiag <- diag(OmegaM1[[iC]] %*% d2Omega.D1[iIndex,iIndex] %*% OmegaM1.dOmega.2)
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * (1 - leverage[iC,iIndex]))
                        }

                        if(test.Omega3a || test.Omega3b){
                            iDiag <- diag(OmegaM1.dOmega.1 %*% OmegaM1[[iC]] %*% d2Omega.D2[iIndex,iIndex])
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + 1/2 * sum(iDiag * (1 - leverage[iC,iIndex]))
                        }

                        if(test.mu1a || test.mu1b){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(d2mu.D1[iC,iIndex] %*% OmegaM1[[iC]] * dmu[[iName2]][iC,iIndex])
                        }
                        
                        if(test.mu2a || test.mu2b){
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] + sum(dmu[[iName1]][iC,iIndex] %*% OmegaM1[[iC]] * d2mu.D2[iC,iIndex])
                        }
                        
                        if(test.mu3){                            
                            dInfo[iName1,iName2,iNameD] <- dInfo[iName1,iName2,iNameD] - sum(dmu[[iName1]][iC,iIndex] %*% OmegaM1.dOmega.D %*% OmegaM1[[iC]] * dmu[[iName2]][iC,iIndex])
                        }
                    }
                }
            }
            
        }
        ## *** Symmetrize
        dInfo[,,iNameD] <- symmetrize(dInfo[,,iNameD], update.upper = NULL)
    }

    ### ** export
    return(dInfo)
}


