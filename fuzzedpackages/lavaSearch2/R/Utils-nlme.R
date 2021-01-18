### utils-nlme.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 15 2017 (17:29) 
## Version: 
## Last-Updated: dec 10 2018 (23:28) 
##           By: Brice Ozenne
##     Update #: 671
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .coef2
#' @title Export Mean and Variance Coefficients
#' @description Export mean and variance coefficients
#' from a \code{lm}, \code{gls}, or \code{lme} object.
#' @name coef2-internal
#'
#' @param object a \code{lm}, \code{gls} or \code{lme} object.
#' @param name.Y [character] the name of the endogenous variable. Used to name certain variance parameters.
#' 
#' @details The variance coefficients that are exported are the residual variance of each outcome. 
#' This is \eqn{\sigma^2} for the first one and \eqn{k^2 \sigma^2} for the remaining ones.
#'
#' @return A numeric vector named with the names of the coefficient with three attributes:
#' \itemize{
#' \item mean.coef: the name of the mean coefficients.
#' \item var.coef: the name of the variance coefficients.
#' \item cor.coef:  the name of the correlation coefficients.
#' }
#'
#' @concept extractor
#' @keywords internal
`.coef2` <-
    function(object) UseMethod(".coef2")

## * .coef2.lm
#' @rdname coef2-internal
.coef2.lm <- function(object){
    coef.object <- coef(object)
    p <- c(coef.object,sigma2=sigma(object)^2)
    attr(p, "mean.coef") <- names(coef.object)
    attr(p, "var.coef") <- "sigma2"
    attr(p, "cor.coef") <- NULL
    return(p)
}

## * .coef2.gls
#' @rdname coef2-internal
.coef2.gls <- function(object){

    ## *** mean coefficients
    mean.coef <- stats::coef(object)

    ## *** variance coefficients
    var.coef <- c(sigma2 = stats::sigma(object)^2)
    if(!is.null(object$modelStruct$varStruct)){
        var.coef <- c(var.coef,
                      stats::coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE)^2)          
    }

    ## *** covariance coefficients
    if(!is.null(object$modelStruct$corStruct)){
        cor.coef <- stats::coef(object$modelStruct$corStruct, unconstrained = FALSE)

        n.var <- length(var.coef)
        n.cor <- length(cor.coef)

        ## check unstructured covariance matrix
        if(!is.null(object$modelStruct$varStruct) && ((n.var*(n.var-1))/2 == n.cor)){

            vecgroup <- attr(unclass(object$modelStruct$corStruct), "group")
            veccov.cor <- unname(unlist(attr(object$modelStruct$corStruct, "covariate")))
            veccov.var <- attr(object$modelStruct$varStruct, "groups")

            table.covvar <- table(veccov.cor,veccov.var)
            newlevels.cor <- colnames(table.covvar)[apply(table.covvar, 1, which.max)]
            veccov.cor2 <- factor(veccov.cor, levels = 0:max(veccov.cor), labels = newlevels.cor)
            
            if(identical(as.character(veccov.cor2),as.character(veccov.var))){

                cor.coefName <- apply(utils::combn(newlevels.cor, m = 2), MARGIN = 2, FUN = paste, collapse = "")
                names(cor.coef) <- paste0("corCoef",cor.coefName)

            }else{
                names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
            }
            
            
        }else{
            names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
        }
        
        
    }else{
        cor.coef <- NULL
    }

    p <- c(mean.coef, cor.coef, var.coef)
    attr(p, "mean.coef") <- names(mean.coef)
    attr(p, "var.coef") <- names(var.coef)
    attr(p, "cor.coef") <- names(cor.coef)
    return(p)
}




## * .coef2.lme
#' @rdname coef2-internal
.coef2.lme <- function(object){

     ## *** mean coefficients
    mean.coef <- nlme::fixef(object)

    ## *** variance coefficients
    var.coef <- c(sigma2 = stats::sigma(object)^2)
    if(!is.null(object$modelStruct$varStruct)){
       var.coef <- c(var.coef,
                      stats::coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE)^2)   
    }

    ## *** random effect coefficients
    random.coef <- as.double(nlme::getVarCov(object))    
    names(random.coef) <- paste0("ranCoef",1:length(random.coef))

     ## *** correlation coefficients
    if(!is.null(object$modelStruct$corStruct)){
        cor.coef <- stats::coef(object$modelStruct$corStruct, unconstrained = FALSE)
        names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
    }else{
        cor.coef <- NULL
    }
    
    p <- c(mean.coef, cor.coef, var.coef, random.coef)
    attr(p, "mean.coef") <- names(mean.coef)
    attr(p, "var.coef") <- names(var.coef)
    attr(p, "cor.coef") <- names(cor.coef)
    attr(p, "ran.coef") <- names(random.coef)
    return(p)
}

## * .getFormula2
`.getFormula2` <-
    function(object) UseMethod(".getFormula2")

## * .getFormula2.gls
.getFormula2.gls <- function(object){
    return(evalInParentEnv(object$call$model))
}

## * .getFormula2.lme
.getFormula2.lme <- function(object){
    return(evalInParentEnv(object$call$fixed))
}


## * .getCluster2
#' @title Reconstruct the Cluster Variable from a nlme Model
#' @description Reconstruct the cluster variable from a nlme model.
#' @name getCluster2-internal
#'
#' @param object a \code{gls} or \code{lme} object.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param data [data.frame] the data set.
#' @param ... [internal] Only used by the generic method.
#'  
#' @return A list containing:
#' \itemize{
#' \item cluster: the cluster index for each observation.
#' \item n.cluster: the number of clusters.
#' }
#'
#' @concept extractor
#' @keywords internal
`.getCluster2` <-
    function(object, ...) UseMethod(".getCluster2")

## * .getCluster2.gls
#' @rdname getCluster2-internal
.getCluster2.gls <- function(object, cluster, data, ...){

### ** get cluster
    if(is.null(object$modelStruct$corStruct)){
        if(missing(cluster)){
            stop("cluster must be specified for gls object with no correlation structure \n")
        }        
        if(length(cluster) == 1 && is.character(cluster)){
            if(cluster %in% names(data) == FALSE){
                stop("Variable \"",cluster,"\" not in data \n")
            }
            ## cluster <- as.numeric(factor(data[[cluster]], levels = unique(data[[cluster]])))
            cluster <- as.numeric(as.factor(data[[cluster]]))
            
        }else if(length(cluster)==NROW(data)){
            cluster <- as.numeric(as.factor(cluster))
        }else{
            stop("length of cluster and data do not match \n")
        }
    }else{
        cluster <- as.numeric(nlme::getGroups(object))
    }
    n.cluster <- length(unique(cluster))

### ** reorder cluster according to the data ordering
    levels.cluster <- unique(cluster)
    cluster <- as.numeric(factor(cluster, levels = levels.cluster))

### ** export
    return(list(cluster = cluster,
                levels.cluster = levels.cluster,
                n.cluster = n.cluster))
}

## * .getCluster2.lme
#' @name getCluster2-internal
.getCluster2.lme <- function(object, ...){

### ** get cluster
    if(NCOL(object$groups)!=1){
        stop("cannot only handle one random effect \n")
    }
    cluster <- as.numeric(nlme::getGroups(object))
    n.cluster <- length(unique(cluster))

### ** reorder cluster according to the data ordering
    levels.cluster <- unique(cluster)
    cluster <- as.numeric(as.factor(cluster))
    ## cluster <- as.numeric(factor(cluster, levels = levels.cluster))

### ** export
    return(list(cluster = cluster,
                levels.cluster = levels.cluster,
                n.cluster = n.cluster))
}

## * .getIndexOmega2
#' @title Extract the name of the endogenous variables
#' @description Extract the name of the endogenous variables from a nlme model.
#' @name getIndexOmega2-internal
#'
#' @param object a \code{gls} or \code{lme} object.
#' @param param [numeric vector] the mean and variance coefficients.
#' @param attr.param [character vector] the type of each coefficients (mean or variance).
#' @param name.Y [character] name of the endogenous variable.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param data [data.frame] the data set.
#' @param ... [internal] Only used by the generic method.
#'  
#' @return A list containing:
#' \itemize{
#' \item index.Omega: [list of integer vector] For each cluster of observations,
#' the index of the endogenous variable relative to each observation.
#' \item n.endogenous: [integer] the number of endogenous variables.
#' \item name.endogenous: [character vector] the name of the endogenous variables.
#' \item ref.group: [character vector] the levels of the variable defining the variance component in a generic covariance matrix.
#' }
#'
#' @concept extractor
#' @keywords internal
`.getIndexOmega2` <-
    function(object, ...) UseMethod(".getIndexOmega2")

## * .getIndexOmega2.gls
#' @rdname getIndexOmega2-internal
.getIndexOmega2.gls <- function(object, param, attr.param, 
                                name.Y, cluster, levels.cluster, data){

    class.cor <- class(object$modelStruct$corStruct)
    class.re <- class(object$modelStruct$reStruct)
    if("NULL" %in% class.cor == FALSE){
        formula.cor <- attr(object$modelStruct$corStruct,"formula")
        varIndex.cor <- all.vars(nlme::getCovariateFormula(formula.cor))
    }else{
        varIndex.cor <- NULL
    }

    class.var <- class(object$modelStruct$varStruct)
    if("NULL" %in% class.var == FALSE){
        formula.var <- attr(object$modelStruct$varStruct,"formula")
        varIndex.var <- all.vars(nlme::getCovariateFormula(formula.var))
        groupValue.var <- attr(object$modelStruct$varStruct,"groups")
        if("NULL" %in% class.cor == FALSE || "NULL" %in% class.re == FALSE){ ## nlme automaticly sort data when corStruct or reSruct
            groupValue.var <- groupValue.var[order(order(cluster))] ## undo automatic sort
        } 
    }else{
        varIndex.var <- NULL
    }

    ## ** Check admissible var-cor structure
    validClass.cor <- c("NULL","corCompSymm","corSymm","corStruct")
    if(any(class.cor %in% validClass.cor == FALSE)){
        stop("Can only handle corStruct of class \"corCompSymm\" or \"corSymm\"\n")
    }

    validClass.var <- c("NULL","varIdent","varFunc")
    if(any(class.var %in% validClass.var == FALSE)){
        stop("Can only handle varStruct of class \"varIdent\"\n")
    }

    ## ** Check compatible ordering between var and cor
    if(length(varIndex.cor) > 0 && length(varIndex.cor) > 0 && !identical(varIndex.cor,varIndex.cor)){
        stop("Inconsistency between the left hand side of the formula in corStruct and the left hand side of the formula in varStruct. \n",
             "it should be something like: correlation = corStruct(form = ~index|groupA) \n",
             "                             weight = varStruct(form = ~index|groupB) \n")
    }
    
    ## ** Identify the index and name of the endogenous variables
    if("NULL" %in% class.var && "NULL" %in% class.cor){ ## basic lme models or lm-ish models
        ## order of the variables does not matter
        index.Omega <- tapply(cluster,cluster,function(iC){list(1:length(iC))})
        norm <- FALSE
    }else if("NULL" %in% class.var && "corCompSymm" %in% class.cor){
        ## order of the variables does not matter
        index.Omega <- attr(object$modelStruct$corStruct, "covariate")[levels.cluster]
        norm <- TRUE
    }else if(length(varIndex.cor)!=0){ 
        ## order of the variables matters: use index variable in corStruct
        index.Omega <- attr(object$modelStruct$corStruct, "covariate")[levels.cluster]
        norm <- TRUE
    }else if(length(varIndex.var)!=0){
        ## order of the variables matters: using index variable in varStruct
        index.tempo <- data[[varIndex.var]]
        if(!is.numeric(index.tempo)){
            stop("The variable in the left hand side of the formula in varStruct must be numeric \n")
        }
        if(!is.null(object$na.action)){
            index.tempo <- index.tempo[-object$na.action]
        }
        index.Omega <- tapply(index.tempo, cluster, function(iC){list(iC)})
        norm <- TRUE
    }else{
        ## order of the variables matters: check missing values
        if("NULL" %in% class.var == FALSE){
            test.duplicated <- unique(unlist(tapply(groupValue.var, cluster, duplicated)))
        }else{
            groupValue.cor <- attr(object$modelStruct$corStruct, "covariate")
            ref.tempo <- unname(sort(groupValue.cor[[1]]))
            test.identical <- unique(unlist(lapply(groupValue.cor,function(x){ # x <- groupValue.cor[[1]]
                identical(unname(sort(x)),ref.tempo)
            })))
            test.duplicated <- TRUE
        }
        ## recover order
        if(("NULL" %in% class.var == FALSE) && all(test.duplicated==FALSE)){
            ## from varIdent when no missing values
            tempo <- as.numeric(factor(groupValue.var, levels = attr(object$modelStruct$varStruct,"groupNames")))
            index.Omega <- tapply(tempo, cluster, function(iC){iC})
            norm <- FALSE
        }else{
            ## from the order of the data
            index.Omega <- tapply(cluster,cluster,function(iC){list(1:length(iC))})
            norm <- FALSE
            if("NULL" %in% class.cor == FALSE){
                if(any(test.identical == FALSE)){
                    warning("The residuals covariance matrice is subset based on the ordering of the data\n",
                            "It is safer to define the ordering within group by adding a variable in the left hand side of the formula in corStruct \n",
                            "e.g. correlation = corStruct(form = ~index|group) \n")
                }
            }else{
                warning("The attribution of the repetition number is based on the ordering of the data\n",
                        "It is safer to define the ordering within group by adding a variable in the left hand side of the formula in varStruct \n",
                        "e.g. correlation = varStruct(form = ~index|group) \n")
            }
        }
    }
    attr(object$modelStruct$varStruct,"groupNames")
    attr(object$modelStruct$varStruct,"groups")

    ## ** Normalize index.Omega
    if(norm){
        level.index <- unique(unlist(index.Omega))
        convertion <- setNames(order(level.index), level.index)
        index.Omega <- lapply(index.Omega, function(x){as.double(convertion[as.character(x)])})
    }

    ## ** Define the name and number of endogenous variables
    vecIndex.Omega <- unlist(index.Omega)
    n.endogenous <- max(vecIndex.Omega)
    name.endogenous <- paste0(name.Y,".",1:n.endogenous)

    if("corSymm" %in% class.cor){
        ## second order polynomial equation
        ## m(m-1)/2 = n.cor
        ## i.e. m = 1/2 + sqrt(1+8 n.cor)/2
        cor.coef <- param[attr.param$cor.coef]
        if(n.endogenous != ( 1 + sqrt(1 + 8 * length(cor.coef)) ) / 2){
            stop("The values of ",varIndex.cor," does not match the number of correlation coefficients \n",
                 "i.e. the maximum value of ",varIndex.cor," should equal ( 1 + sqrt(1 + 8 * n.cor.coef)) ) / 2 \n"
                 )
        }
    }

    if("NULL" %in% class.var == FALSE){

        groupValue.var.ordered <- groupValue.var[order(cluster)] ## reorder by cluster
        table.unique <- tapply(1:length(vecIndex.Omega),vecIndex.Omega,function(x){
            length(unique(groupValue.var.ordered[x]))
        })
        if(any(table.unique!=1)){
            stop("The residual covariance matrix should not differ between clusters \n")
        }                
        ref.group <- groupValue.var.ordered[!duplicated(vecIndex.Omega)]
    }else{
        ref.group <- NULL
    }
        
    ## ** Export
    return(list(index.Omega = index.Omega,
                n.endogenous = n.endogenous,
                name.endogenous = name.endogenous,
                ref.group = ref.group))
}

## * .getIndexOmega2.lme
#' @rdname getIndexOmega2-internal
.getIndexOmega2.lme <- .getIndexOmega2.gls

## * .getVarCov2
#' @title Reconstruct the Marginal Variance Covariance Matrix from a nlme Model
#' @description Reconstruct the marginal variance covariance matrix from a nlme model.
#' @name getVarCov2-internal
#'
#' @param object a \code{gls} or \code{lme} object
#' @param param [numeric vector] the mean and variance coefficients.
#' @param attr.param [character vector] the type of each coefficients (mean or variance).
#' @param name.endogenous [character vector] name of each repetition of the endogenous variable. 
#' @param n.endogenous [integer >0] number of repetitions of the endogenous variable.
#' @param ref.group [character vector] the levels of the variable defining the variance component in a generic covariance matrix.
#' @param ... [internal] Only used by the generic method.
#'  
#' @return [matrix] the marginal variance covariance matrix for a full sample.
#'
#' @concept extractor
#' @keywords internal
`.getVarCov2` <-
    function(object, ...) UseMethod(".getVarCov2")

## * .getVarCov2.gls
#' @rdname getVarCov2-internal
.getVarCov2.gls <- function(object, param, attr.param,
                            name.endogenous, n.endogenous, ref.group, ...){

    ## ** Extract information
    var.coef <- param[attr.param$var.coef]
    cor.coef <- param[attr.param$cor.coef]

    ## ** Diagonal terms
    if(length(ref.group)>0){        
        factor.varcoef <- setNames(c(1,var.coef[-1]),
                                   attr(object$modelStruct$varStruct,"groupNames"))
        sigma2.base <- var.coef["sigma2"] * factor.varcoef[ref.group]
    }else{
        sigma2.base <- rep(var.coef["sigma2"], n.endogenous)
    }
    ## re-order according to the order of the correlation coefficients (if possible)
    if(length(cor.coef)>1 & length(var.coef)>1){
        cor.level <- gsub("corCoef","",names(cor.coef))
        var.level <- names(sigma2.base)
        var.permlevel <- .allPermutations(var.level)

        M.try <- apply(var.permlevel, MARGIN = 1, function(iLevel){
            all(apply(utils::combn(iLevel, m = 2), MARGIN = 2, FUN = paste, collapse = "") == cor.level)
        })
        if(any(M.try)){
            sigma2.base <- sigma2.base[var.permlevel[which(M.try),]]
        }
    }
    
    Omega <- diag(as.double(sigma2.base),
                  nrow = n.endogenous, ncol = n.endogenous)

    ## ** Extra-diagonal terms
    if(length(cor.coef)>0){
        index.lower <- which(lower.tri(Omega))
        index.lower.arr <- which(lower.tri(Omega),arr.ind = TRUE)
        vec.sigma.tempo <- apply(index.lower.arr,1,function(x){prod(sqrt(sigma2.base[x]))})        
        Omega[index.lower] <- cor.coef*vec.sigma.tempo
        Omega <- symmetrize(Omega, update.upper = TRUE)
    }    

    ## ** names
    if(all(!duplicated(names(sigma2.base)))){
        dimnames(Omega) <- list(names(sigma2.base), names(sigma2.base))
    }else{
        dimnames(Omega) <- list(name.endogenous, name.endogenous)
    }

    
    ## ** export
    return(Omega)
}

## * .getVarCov2.lme
#' @rdname getVarCov2-internal
.getVarCov2.lme <- function(object, param, attr.param, ...){

    ## ** prepare with gls
    out <- .getVarCov2.gls(object, param = param, attr.param = attr.param, ...)

    ## ** add contribution of the random effect
    ran.coef <- param[attr.param$ran.coef]
    out <- out + ran.coef

    ## ** export
    return(out)    
}

## * .allPermutations
## .allPermutations(1:3)
## .allPermutations(2:3)
.allPermutations <- function(vec){
    X <- lapply(vec, function(x){
        cbind(x, .allPermutations(setdiff(vec, x)))
    })
    return(unname(do.call(rbind,X)))
}

##----------------------------------------------------------------------
### utils-nlme.R ends here
