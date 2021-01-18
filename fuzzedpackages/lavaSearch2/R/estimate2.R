### estimate2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 16 2018 (16:38) 
## Version: 
## Last-Updated: feb 15 2019 (14:08) 
##           By: Brice Ozenne
##     Update #: 864
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * estimate2
#' @title Compute Bias Corrected Quantities.
#' @description Compute bias corrected residuals variance covariance matrix
#' and information matrix.
#' Also provides the leverage values and corrected sample size when adjust.n is set to TRUE.
#' @name estimate2
#' 
#' @keywords internal
.estimate2 <- function(object, epsilon, n.cluster,
                       name.param, name.endogenous, name.meanparam, name.varparam,
                       index.Omega,
                       adjust.Omega, adjust.n, tol, n.iter, trace){

    ## ** Prepare
    Omega <- object$conditionalMoment$Omega
    dmu <- object$conditionalMoment$dmu
    dOmega <- object$conditionalMoment$dOmega
    
    name.hybridparam <- intersect(name.meanparam, name.varparam)

    n.param <- length(name.param)
    n.meanparam <- length(name.meanparam)
    n.varparam <- length(name.varparam)
    n.hybridparam <- length(name.hybridparam)

    n.endogenous <- length(name.endogenous)
    grid.meanparam <- .combination(name.meanparam, name.meanparam)    
    n.grid.meanparam <- NROW(grid.meanparam)
    grid.varparam <- .combination(name.varparam, name.varparam)
    n.grid.varparam <- NROW(grid.varparam)

    ## check low diagonal
    name2num <- setNames(1:n.param,name.param)
    if(!all(name2num[grid.meanparam[,1]]-name2num[grid.meanparam[,2]]>=0)){
        stop("Incorrect allocation of the computation of the information matrix (mean parameter) \n")
    }
    name2num <- setNames(1:n.param,name.param)
    if(!all(name2num[grid.varparam[,1]]-name2num[grid.varparam[,2]]>=0)){
        stop("Incorrect allocation of the computation of the information matrix (variance parameter) \n")
    }
    ##
    
    leverage <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                       dimnames = list(NULL, name.endogenous))
    ls.dmu <- vector(mode = "list", length = n.cluster)
    for(iC in 1:n.cluster){ # iC <- 1
        if(is.null(index.Omega)){            
            leverage[iC,] <- 0
            ls.dmu[[iC]] <- matrix(0, nrow = n.param, ncol = n.endogenous,
                                   dimnames = list(name.param, name.endogenous))
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu[name.meanparam],function(x){x[iC,]}))
        }else{
            leverage[iC,index.Omega[[iC]]] <- 0
            ls.dmu[[iC]] <- matrix(0, nrow = n.param, ncol = length(index.Omega[[iC]]),
                                   dimnames = list(name.param, name.endogenous[index.Omega[[iC]]]))
            ls.dmu[[iC]][name.meanparam,] <- do.call(rbind, lapply(dmu[name.meanparam],function(x){x[iC,index.Omega[[iC]]]}))
        }        
    }
    
    ## ** Initialisation (i.e. first iteration without correction)
    if(any(eigen(Omega)$value<=0)){
        stop("the residual variance-covariance matrix is not positive definite \n")
    }

    if(is.null(index.Omega)){
        n.corrected <- rep(n.cluster, n.endogenous)
    }else{
        n.corrected <- NULL
    }
    ls.Psi <- vector(mode = "list", length = n.cluster)

    Omega.adj <- Omega
    if(!adjust.n){
       epsilon.adj <- epsilon
    }

    if(trace>0){
        cat("* Reconstruct estimated information matrix ")
    }

    iInfo <- .information2(dmu = dmu,
                           dOmega = dOmega,
                           Omega = Omega,
                           n.corrected = n.corrected,
                           leverage = leverage, index.Omega = index.Omega, n.cluster = n.cluster,
                           grid.meanparam = grid.meanparam,
                           n.grid.meanparam = n.grid.meanparam,
                           grid.varparam = grid.varparam,
                           n.grid.varparam = n.grid.varparam,
                           name.param = name.param,
                           n.param = n.param)
    iVcov.param <- try(chol2inv(chol(iInfo)), silent = TRUE)
    if(inherits(iVcov.param, "try-error")){
        iVcov.param <- solve(iInfo)
    }
    if(trace>0){
        cat("- done \n")
    }
    
    ## ** Loop    
    if(adjust.Omega || adjust.n){
        if(trace>0){
            cat("* iterative small sample correction: ")
        }
        iIter <- 0
        iTol <- Inf
        Omega_save <- Omega
        iOmega.adj <- Omega.adj
    }else{
        iIter <- Inf
        iTol <- -Inf        
    }
    
    while(iIter < n.iter & iTol > tol){
        if(trace>0){
            cat("*")
        }

        ## *** Step (i-ii): compute individual bias, expected bias
        Psi <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                      dimnames = list(name.endogenous, name.endogenous))
        M.countCluster <- matrix(0, nrow = n.endogenous, ncol = n.endogenous,
                                 dimnames = list(name.endogenous, name.endogenous))
        for(iC in 1:n.cluster){
            ## individual bias
            ls.Psi[[iC]] <- t(ls.dmu[[iC]])  %*% iVcov.param %*% ls.dmu[[iC]]
            ## cumulated bias            
            if(is.null(index.Omega)){
                Psi <- Psi + ls.Psi[[iC]]
                M.countCluster <- M.countCluster + 1
            }else{
                Psi[index.Omega[[iC]],index.Omega[[iC]]] <- Psi[index.Omega[[iC]],index.Omega[[iC]]] + ls.Psi[[iC]]
                M.countCluster[index.Omega[[iC]],index.Omega[[iC]]] <- M.countCluster[index.Omega[[iC]],index.Omega[[iC]]] + 1
            }
        }

        ## update
        for(iPsi in 1:length(Psi)){
            if(M.countCluster[iPsi]>0){
                Psi[iPsi] <- Psi[iPsi]/M.countCluster[iPsi]
            }
        }
        
        ## *** Step (iii): compute leverage
        if(adjust.n){
            epsilon.adj <- .adjustResiduals(Omega = Omega.adj,
                                            Psi = Psi,
                                            epsilon = epsilon,
                                            index.Omega = index.Omega,
                                            name.endogenous = name.endogenous,
                                            n.endogenous = n.endogenous,
                                            n.cluster = n.cluster)

            leverage <- .adjustLeverage(Omega = Omega.adj,
                                        epsilon = epsilon.adj,
                                        ls.dmu = ls.dmu,
                                        dOmega = dOmega,
                                        vcov.param = iVcov.param,
                                        index.Omega = index.Omega,
                                        name.endogenous = name.endogenous,
                                        n.endogenous = n.endogenous,
                                        name.varparam = name.varparam,
                                        n.varparam = n.varparam,
                                        n.cluster = n.cluster)

            n.corrected <- rep(n.cluster, n.endogenous) - colSums(leverage, na.rm = TRUE)
        }
        
        ## *** Step (v): correct residual covariance matrix, estimates, and derivatives
        if(adjust.Omega){
            ## corrected residual covariance variance
            Omega.adj <- Omega + Psi
            
            ## correct estimates
            object$conditionalMoment <- .adjustMoment(object, Omega = Omega.adj)
            dOmega <- object$conditionalMoment$dOmega
            ## conditionalMoment.adj$param - coef(object)
           
        }

        ## *** Step (vii): expected information matrix
        iInfo <- .information2(dmu = dmu,
                               dOmega = dOmega,
                               Omega = Omega.adj,
                               n.corrected = n.corrected,
                               leverage = leverage,
                               index.Omega = index.Omega,
                               n.cluster = n.cluster,
                               grid.meanparam = grid.meanparam,
                               n.grid.meanparam = n.grid.meanparam,
                               grid.varparam = grid.varparam,
                               n.grid.varparam = n.grid.varparam,
                               name.param = name.param,
                               n.param = n.param)
        iVcov.param <- try(chol2inv(chol(iInfo)), silent = TRUE)
        if(inherits(iVcov.param, "try-error")){
            iVcov.param <- solve(iInfo)
        }
        
        ## *** Update cv
        iIter <- iIter + 1
        iTol <- norm(Omega.adj-Omega_save, type = "F")
        Omega_save <- Omega.adj
        ## cat("Omega.adj: ",Omega.adj," | n:",n.corrected," | iTol:",iTol,"\n")
    }
    
    ## ** Post processing
    if(!is.infinite(iIter)){

        if(iTol > tol){
            warning("small sample correction did not reach convergence after ",iIter," iterations \n")

            if(trace>0){
                cat(" - incomplete \n")
            }
        }else{
            if(trace>0){
                cat(" - done \n")
            }
        }
        
    }

    vcov.param <- try(chol2inv(chol(iInfo)), silent = TRUE)
    if("try-error" %in% class(vcov.param)){
        errorMessage <- vcov.param
        vcov.param <- solve(iInfo)
        attr(vcov.param, "warning") <- errorMessage
    }
    dimnames(vcov.param) <- dimnames(iInfo)

    ## update object
    object$conditionalMoment$Omega <- Omega.adj
    object$dVcov <- list(param = object$conditionalMoment$param,
                         score = NULL,
                         vcov.param = vcov.param,
                         dVcov.param = NULL,
                         Omega = Omega.adj,
                         residuals = epsilon.adj,
                         leverage = leverage,
                         n.corrected = rep(n.cluster, n.endogenous) - colSums(leverage, na.rm = TRUE),
                         opt = list(objective = iTol, iterations = iIter, convergence = (iTol <= tol), grid.meanparam = grid.meanparam, grid.varparam = grid.varparam))

    ## ** Export
    return(object)
}

## * .adjustResiduals
.adjustResiduals <- function(Omega, Psi, epsilon,
                             index.Omega,
                             name.endogenous, n.endogenous, n.cluster){

    if(is.null(index.Omega)){ ## no missing values
        
        Omega.chol <- matrixPower(Omega, symmetric = TRUE, power = 1/2)
        H <- Omega %*% Omega - Omega.chol %*% Psi %*% Omega.chol
        HM1 <- tryCatch(matrixPower(H, symmetric = TRUE, power = -1/2), warning = function(w){w})
        if(inherits(HM1,"warning")){
            stop("Cannot compute the adjusted residuals \n",
                 "Estimated bias too large compared to the estimated variance-covariance matrix \n",
                 "Consider setting argument \'adjust.n\' to FALSE when calling sCorrect \n")
        }
        epsilon.adj <- epsilon %*% Omega.chol %*% HM1 %*% Omega.chol
        
    }else{ ## missing values
        
        epsilon.adj <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                              dimnames = list(NULL, name.endogenous))

        for(iC in 1:n.cluster){
            iIndex <- index.Omega[[iC]]
            iOmega <- Omega[iIndex,iIndex,drop=FALSE]
            iOmega.chol <- matrixPower(iOmega, symmetric = TRUE, power = 1/2)
            iH <- iOmega %*% iOmega - iOmega.chol %*% Psi[iIndex,iIndex,drop=FALSE] %*% iOmega.chol
            iHM1 <- tryCatch(matrixPower(iH, symmetric = TRUE, power = -1/2), warning = function(w){w})
            if(inherits(iHM1,"warning")){
                stop("Cannot compute the adjusted residuals \n",
                     "Estimated bias too large compared to the estimated variance-covariance matrix \n",
                     "Consider setting argument \'adjust.n\' to FALSE when calling sCorrect \n")
            }
            epsilon.adj[iC,iIndex] <- epsilon[iC,iIndex] %*% iOmega.chol %*% iHM1 %*% iOmega.chol
        }
        
    }
    dimnames(epsilon.adj) <- list(NULL,name.endogenous)
    return(epsilon.adj)
}

## * .adjustLeverage
.adjustLeverage <- function(Omega, epsilon, ls.dmu, dOmega, vcov.param,
                            index.Omega,
                            name.endogenous, n.endogenous, name.varparam, n.varparam, n.cluster){

    ## ** prepare
    leverage <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                       dimnames = list(NULL, name.endogenous))

    if(is.null(index.Omega)){
        iIndex <- 1:n.endogenous
        iOmegaM1 <- chol2inv(chol(Omega)) ## solve(Omega)
        iOmegaM1.dOmega.OmegaM1 <- lapply(dOmega, function(x){iOmegaM1 %*% x %*% iOmegaM1})
    }

    ## ** compute
    for(iC in 1:n.cluster){                 # iC <- 1
        if(!is.null(index.Omega)){
            iIndex <- index.Omega[[iC]]
            iOmegaM1 <- chol2inv(chol(Omega[iIndex,iIndex,drop=FALSE]))
            iOmegaM1.dOmega.OmegaM1 <- lapply(dOmega, function(x){iOmegaM1 %*% x[iIndex,iIndex] %*% iOmegaM1})
        }
        ## derivative of the score regarding Y
        scoreY <- ls.dmu[[iC]] %*% iOmegaM1

        for(iP in 1:n.varparam){ ## iP <- 1
            scoreY[name.varparam[iP],] <- scoreY[name.varparam[iP],] + 2 * epsilon[iC,iIndex] %*% iOmegaM1.dOmega.OmegaM1[[name.varparam[iP]]]
        }
        ## leverage
        leverage[iC,iIndex] <- colSums(vcov.param %*% ls.dmu[[iC]] * scoreY) ## NOTE: dimensions of ls.dmu and scoreY matches even when there are missing values
                                        # same as
                                        # diag(t(ls.dmu[[iC]])  %*% iVcov.param %*% scoreY)
    }

    return(leverage)            
}

## * .adjustMoment
`.adjustMoment` <-
    function(object, ...) UseMethod(".adjustMoment")

## * .adjustMoment.lm
.adjustMoment.lm <- function(object, Omega){

    object$conditionalMoment$param["sigma2"] <- as.double(Omega)
    return(object$conditionalMoment)
    
}

## * .adjustMoment.gls
.adjustMoment.gls <- function(object, Omega, ...){

    ## ** extract information
    class.cor <- object$conditionalMoment$skeleton$class.cor
    class.var <- object$conditionalMoment$skeleton$class.var
    name.corcoef <- object$conditionalMoment$skeleton$name.corcoef
    name.otherVar <- object$conditionalMoment$skeleton$name.otherVar
    name.varcoef <- object$conditionalMoment$skeleton$name.varcoef
    ref.group <- object$conditionalMoment$skeleton$ref.group
    M.corcoef <- object$conditionalMoment$skeleton$M.corcoef
    name.endogenous <- object$conditionalMoment$skeleton$name.endogenous
    n.endogenous <- object$conditionalMoment$skeleton$n.endogenous
    
    ## ** identify parameters

    if(identical(class.var, "NULL")){
        object$conditionalMoment$param["sigma2"] <- mean(diag(Omega))
    }else{            
        index.Sigma2 <- which(ref.group %in% name.otherVar == FALSE)
        object$conditionalMoment$param["sigma2"] <- mean(diag(Omega)[index.Sigma2])

        vec.k <- tapply(diag(Omega)/Omega[index.Sigma2,index.Sigma2], ref.group, mean)            
        object$conditionalMoment$param[name.otherVar] <- vec.k[name.otherVar]
    }

    if(identical(class.cor, "NULL")){
        ## do nothing
    }else if("corCompSymm" %in% class.cor){
        object$conditionalMoment$param[name.corcoef] <- mean(stats::cov2cor(Omega)[lower.tri(Omega)])
    }else if("corSymm" %in% class.cor){
        vec.cor <- tapply(stats::cov2cor(Omega)[lower.tri(Omega)],
                          M.corcoef[lower.tri(Omega)],
                          mean)            
        object$conditionalMoment$param[name.corcoef] <- vec.cor[name.corcoef]
    } 

    ## ** update conditional moments
    object$conditionalMoment$Omega <- .getVarCov2(object,
                                                  param = object$conditionalMoment$param,
                                                  attr.param = attributes(object$conditionalMoment$param),
                                                  name.endogenous = name.endogenous,
                                                  n.endogenous = n.endogenous,
                                                  ref.group = ref.group)
    
    ## ** update first derivative of the conditional variance
    object$conditionalMoment$dOmega <- skeletonDtheta(object, class.cor = class.cor, class.var = class.var, 
                                                      sigma2.base0 = object$conditionalMoment$skeleton$sigma2.base0,
                                                      Msigma2.base0 = object$conditionalMoment$skeleton$Msigma2.base0,
                                                      M.corcoef = M.corcoef, ref.group = ref.group,
                                                      index.lower.tri = object$conditionalMoment$skeleton$index.lower.tri,
                                                      indexArr.lower.tri = object$conditionalMoment$skeleton$indexArr.lower.tri,
                                                      name.endogenous =  name.endogenous, n.endogenous = n.endogenous,
                                                      cluster = object$conditionalMoment$skeleton$cluster,
                                                      n.cluster = object$conditionalMoment$skeleton$n.cluster,
                                                      var.coef = object$conditionalMoment$param[name.varcoef],
                                                      name.varcoef = name.varcoef, name.otherVar = name.otherVar,
                                                      n.varcoef = object$conditionalMoment$skeleton$n.varcoef,
                                                      cor.coef = object$conditionalMoment$param[name.corcoef],
                                                      name.corcoef = name.corcoef,
                                                      n.corcoef = object$conditionalMoment$skeleton$n.corcoef,
                                                      update.mean = FALSE, update.variance = TRUE, ...)$dOmega

    ## ** export
    ## names(object$conditionalMoment)
    ## object$conditionalMoment$param["sigma2"] <- as.double(Omega)
    return(object$conditionalMoment)
    
}

## * .adjustMoment.lme
.adjustMoment.lme <- function(object, Omega){

    name.rancoef <- attr(object$conditionalMoment$param,"ran.coef")
    
    ## ** Identify random effect
    if(!identical(object$conditionalMoment$skeleton$class.cor,"NULL")){
        stop("Does not know how to identify the correlation coefficients when corStruct is not NULL \n")
    }
    object$conditionalMoment$param[name.rancoef] <- mean(Omega[lower.tri(Omega)])

    ## ** save derivative regarding random effect
    save <- object$conditionalMoment$dOmega$ranCoef1

    ## ** compute moments 
    conditionalMoment <- .adjustMoment.gls(object, Omega = Omega - object$conditionalMoment$param["ranCoef1"],
                                           name.rancoef = name.rancoef)

    ## ** restaure derivative regarding random effect
    conditionalMoment$dOmega$ranCoef1 <- save
    return(conditionalMoment)
}

## * .adjustMoment.lvmfit
.adjustMoment.lvmfit <- function(object, Omega){

    ## ** extract info
    n.endogenous <- NROW(Omega)
    df.param <- object$conditionalMoment$df.param
    
    index.matrix <- object$conditionalMoment$adjustMoment$index.matrix
    index.Psi <- object$conditionalMoment$adjustMoment$index.Psi
    A <- object$conditionalMoment$adjustMoment$A
    name.var <- object$conditionalMoment$adjustMoment$name.var
    n.rhs <- object$conditionalMoment$adjustMoment$n.rhs
    index.LambdaB <- object$conditionalMoment$adjustMoment$index.LambdaB
    name.endogenous <- object$conditionalMoment$adjustMoment$name.endogenous
    name.latent <- object$conditionalMoment$adjustMoment$name.latent
    
    skeleton <- object$conditionalMoment$skeleton

    param <- object$conditionalMoment$param
    Lambda <- object$conditionalMoment$value$Lambda
    iIB <- object$conditionalMoment$value$iIB
    iIB.Lambda <- object$conditionalMoment$value$iIB.Lambda
    dLambda <- object$conditionalMoment$dMoment.init$dLambda
    dB <- object$conditionalMoment$dMoment.init$dB

    ## ** right hand side of the equation
    eq.rhs <- Omega[index.matrix$index]

    ## ** left hand side of the equation
    if(NROW(index.Psi)>0){
        n.index.Psi <- NROW(index.Psi)
        n.latent <- NROW(skeleton$Psi)        
        Z <- iIB %*% Lambda

        ## A = t(Z) Psi Z + Sigma
        ## (t(Z) Psi Z)_{ij} = \sum_{k,l} Z_{k,i} Psi_{k,l} Z_{l,j}
        for(iIndex in 1:n.rhs){ # iIndex <- 1
            iRow <- index.matrix[iIndex,"row"]
            iCol <- index.matrix[iIndex,"col"]
            for(iPsi in 1:n.index.Psi){
                iRowPsi <- index.Psi[iPsi,"row"]
                iColPsi <- index.Psi[iPsi,"col"]
                A[iIndex,skeleton$Psi[iRowPsi,iColPsi]] <- A[iIndex,skeleton$Psi[iRowPsi,iColPsi]] + Z[iRowPsi,iRow]*Z[iColPsi,iCol]
            }
        }
    }

    ## ** solve equation
    ## microbenchmark::microbenchmark(svd = {asvd <- svd(A) ; asvd$v %*% diag(1/asvd$d) %*% t(asvd$u) %*% eq.rhs;},
    ## qr = qr.coef(qr(A), eq.rhs),
    ## Rcpp = OLS_cpp(A, eq.rhs),
    ## RcppTry = try(OLS_cpp(A, eq.rhs)[,1], silent = TRUE),
    ## Rcpp2 = OLS2_cpp(A, eq.rhs),
    ## OLS1 = solve(crossprod(A), crossprod(A, eq.rhs)),
    ## OLS2 = solve(t(A) %*% A) %*% t(A) %*% eq.rhs,
    ## OLS_stats = stats::lsfit(x = A, y = eq.rhs),
    ## OLS_LINPACK = .Call(stats:::C_Cdqrls, x = A, y = eq.rhs, tolerance = 1e-7, FALSE)$coefficients, times = 500)
    if(lava.options()$method.estimate2=="svd"){
        asvd <- svd(A)
        iSolution <- try((asvd$v %*% diag(1/asvd$d) %*% t(asvd$u) %*% eq.rhs)[,1], silent = TRUE)
    }else if(lava.options()$method.estimate2=="ols"){
        iSolution <- try(OLS_cpp(A, eq.rhs)[,1], silent = TRUE)
    }else{
        stop("unknown OLS methods \n")
    }
    
    if(inherits(iSolution, "try-error")){
        if(abs(det(t(A) %*% A)) <  1e-10){            
            stop("Singular matrix: cannot update the estimates \n")
        }else{
            stop(iSolution)
        }
    }

    ## ** update parameters in conditional moments
    object$conditionalMoment$param[name.var] <- setNames(iSolution, name.var)

    ## ** update conditional moments
    object$conditionalMoment$skeleton$toUpdate <- object$conditionalMoment$adjustMoment$toUpdate
    object$conditionalMoment$value <- skeleton.lvmfit(object,
                                                      param = param,
                                                      data = NULL,
                                                      name.endogenous = name.endogenous,
                                                      name.latent = name.latent)
    object$conditionalMoment$Omega <- Omega
    

    ## ** update first derivative of the conditional variance
    if(length(index.LambdaB)>0){
        object$conditionalMoment$dMoment.init$toUpdate[] <- FALSE
        object$conditionalMoment$dMoment.init$toUpdate[index.LambdaB] <- TRUE

        object$conditionalMoment$dOmega <- skeletonDtheta.lvmfit(object,
                                                                name.endogenous = name.endogenous,
                                                                name.latent = name.latent)$dOmega
    }
    ## ** export
    return(object$conditionalMoment)
}

##----------------------------------------------------------------------
### estimate2.R ends here
