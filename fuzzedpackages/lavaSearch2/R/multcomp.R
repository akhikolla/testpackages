
### mlf2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (12:56) 
## Version: 
## Last-Updated: jun 14 2019 (13:44) 
##           By: Brice Ozenne
##     Update #: 556
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## User 
### Code:


## * estfun
#' @title Extract Empirical Estimating Functions (lvmfit Object)
#' @description Extract the empirical estimating functions of a lvmfit object.
#' This function is for internal use but need to be public to enable its use by \code{multcomp::glht}.
#' @name estfun
#' 
#' @param x an \code{lvmfit} object.
#' @param ... arguments passed to methods.
#'
#' @details This function enables to use the \code{glht} function with lvmfit object.
#' Otherwise when calling \code{multcomp:::vcov.mmm} then \code{sandwich::sandwich} and then \code{sandwich::meat}, \code{sandwich::meat} will complain that \code{estfun} is not defined for \code{lvmfit} objects.
#'
#' @examples
#' library(multcomp)
#' 
#' #### generative model ####
#' mSim <- lvm(X ~ Age + 0.5*Treatment,
#'             Y ~ Gender + 0.25*Treatment,
#'             c(Z1,Z2,Z3) ~ eta, eta ~ 0.75*treatment,
#'             Age[40:5]~1)
#' latent(mSim) <- ~eta
#' categorical(mSim, labels = c("placebo","SSRI")) <- ~Treatment
#' categorical(mSim, labels = c("male","female")) <- ~Gender
#'
#' #### simulate data ####
#' n <- 5e1
#' set.seed(10)
#' df.data <- lava::sim(mSim, n = n, latent = FALSE)
#'
#' #### fit separate models ####
#' lmX <- lm(X ~ Age + Treatment, data = df.data)
#' lvmY <- estimate(lvm(Y ~ Gender + Treatment), data = df.data)
#' lvmZ <- estimate(lvm(c(Z1,Z2,Z3) ~ eta, eta ~ Treatment), 
#'                  data = df.data)
#'
#' #### create mmm object #### 
#' e.mmm <- mmm(X = lmX, Y = lvmY, Z = lvmZ)
#'
#' #### create contrast matrix ####
#' resC <- createContrast(e.mmm, var.test = "Treatment", add.variance = FALSE)
#'
#' #### adjust for multiple comparisons ####
#' e.glht <- glht(e.mmm, linfct = resC$mlf)
#' summary(e.glht)
#' @concept multiple comparison

## * estfun.lvmfit
#' @rdname estfun
#' @method estfun lvmfit
#' @export
estfun.lvmfit <- function(x, ...){
    U <- lava::score(x, indiv = TRUE)
    return(U)
}

## * estfun.gls
#' @rdname estfun
#' @method estfun gls
#' @export
estfun.gls <- function(x, ...){
    if(inherits(x,"gls2")){
        U <- score2(x)
    }else{
        U <- score2(x, bias.correct = FALSE)
    }
    return(U)
}

## * estfun.lme
#' @rdname estfun
#' @method estfun lme
#' @export
estfun.lme <- function(x, ...){
    if(inherits(x,"lme2")){
        U <- score2(x)
    }else{
        U <- score2(x, bias.correct = FALSE)
    }
    return(U)
}

## * Documentation - glht2
#' @title General Linear Hypothesis
#' @description Test general linear hypotheses and across latent variable models with small sample corrections.
#' @name glht2
#' 
#' @param model a \code{lvmfit} or \code{mmm} object.
#' The \code{mmm} object can only contain lm/gls/lme/lvmfit objects.
#' @param linfct [matrix or vector of character] the linear hypotheses to be tested. Same as the argument \code{par} of \code{\link{createContrast}}.
#' @param rhs [vector] the right hand side of the linear hypotheses to be tested.
#' @param bias.correct [logical] should the standard errors of the coefficients be corrected for small sample bias?
#' @param df [logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
#' @param robust [logical] should robust standard error be used? 
#' Otherwise rescale the influence function with the standard error obtained from the information matrix.
#' @param cluster  [integer vector] the grouping variable relative to which the observations are iid.
#'
#' @details
#' Whenever the argument linfct is not a matrix, it is passed to the function \code{createContrast} to generate the contrast matrix and, if not specified, rhs. \cr \cr
#'
#' Since only one degree of freedom can be specify in a glht object and it must be an integer, the degree of freedom of the denominator of an F test simultaneously testing all hypotheses is retained, after rounding. \cr \cr
#'
#' Argument rhs and null are equivalent.
#' This redondance enable compatibility between \code{lava::compare}, \code{compare2}, \code{multcomp::glht}, and \code{glht2}.
#' @return A \code{glht} object.
#' 
#' @seealso
#' \code{\link{createContrast}} to create contrast matrices. \cr
#' \code{\link{sCorrect}} to pre-compute quantities for the small sample correction.
#' 
#' @concept multiple comparisons
#'
#' @examples
#' library(multcomp)
#' 
#' ## Simulate data
#' mSim <- lvm(c(Y1,Y2,Y3)~ beta * eta, Z1 ~ E, Z2 ~ E, Age[40:5]~1)
#' latent(mSim) <- "eta"
#' set.seed(10)
#' n <- 1e2
#'
#' df.data <- lava::sim(mSim, n, latent = FALSE, p = c(beta = 1))
#'
#' #### Inference on a single model ####
#' e.lvm <- estimate(lvm(Y1~E), data = df.data)
#' summary(glht2(e.lvm, linfct = c("Y1~E + Y1","Y1")))
#' 
#' #### Inference on separate models ####
#' ## fit separate models
#' lmX <- lm(Z1 ~ E, data = df.data)
#' lvmY <- estimate(lvm(Z2 ~ E + Age), data = df.data)
#' lvmZ <- estimate(lvm(c(Y1,Y2,Y3) ~ eta, eta ~ E), 
#'                  data = df.data)
#'
#' #### create mmm object #### 
#' e.mmm <- mmm(X = lmX, Y = lvmY, Z = lvmZ)
#'
#' #### create contrast matrix ####
#' resC <- createContrast(e.mmm, var.test = "E", add.variance = TRUE)
#'
#' #### adjust for multiple comparisons ####
#' e.glht2 <- glht2(e.mmm, linfct = resC$contrast, df = FALSE)
#' summary(e.glht2)
#'
#' @concept multiple comparison
#' @export
`glht2` <-
    function(model, linfct, rhs,
             bias.correct, df, robust, cluster) UseMethod("glht2")


## * glht2.lvmfit
#' @rdname glht2
#' @export
glht2.lvmfit <- function(model, linfct, rhs = 0,
                         bias.correct = TRUE, df = TRUE, robust = FALSE, cluster = NULL){

    if(robust==FALSE && !is.null(cluster)){
        stop("Argument \'cluster\' must be NULL when argument \'robust\' is FALSE \n")
    }
    
    ### ** define contrast matrix
    if(!is.matrix(linfct)){
        resC <- createContrast(model, par = linfct)
        linfct <- resC$contrast
        if("rhs" %in% names(match.call()) == FALSE){
            rhs <- resC$null
        }
    }

### ** pre-compute quantities for the small sample correction
    if(!inherits(model,"lvmfit2")){
        sCorrect(model, df = df) <- bias.correct
    }

### ** Wald test with small sample correction
    name.param <- colnames(linfct)
    n.param <- NCOL(linfct)
    n.hypo <- NROW(linfct)

    resWald <- compare2(model, contrast = linfct, null = rhs, as.lava = FALSE)
    ## update name according to multcomp, i.e. without second member
    rownames(linfct) <- .contrast2name(linfct, null = NULL) 

### ** Global degree of freedom
    if(df){
        df.global <- round(resWald["global","df"], digits = 0)
    }else{
        df.global <- 0
    }
    
### ** compute variance-covariance matrix
    if(robust){
        vcov.model <- crossprod(iid2(model, cluster = cluster))
    }else{
        vcov.model <- vcov2(model)
    }

### ** convert to the appropriate format
    out <- list(model = model,
                linfct = linfct,
                rhs = unname(rhs),
                coef = coef(model),
                vcov = vcov.model,
                df = df.global,
                alternative = "two.sided",
                type = NULL,
                robust = robust,
                bias.correct = bias.correct)
    class(out) <- c("glht2","glht")
        
    ### ** export
    return(out)
}


## * glht2.mmm
#' @rdname glht2
#' @export
glht2.mmm <- function (model, linfct, rhs = 0,
                       bias.correct = TRUE, df = TRUE, robust = FALSE, cluster = NULL){


    if(robust==FALSE && !is.null(cluster)){
        stop("Argument \'cluster\' must be NULL when argument \'robust\' is FALSE \n")
    }
    
    ### ** check the class of each model
    n.model <- length(model)
    name.model <- names(model)    
    if(is.null(name.model)){
        stop("Argument \'model\' must be named list. \n")
    }
    
    test.lm <- sapply(model, inherits, what = "lm")
    test.gls <- sapply(model, inherits, what = "gls")
    test.lme <- sapply(model, inherits, what = "lme")
    test.lvmfit <- sapply(model, inherits, what = "lvmfit")
    if(any(test.lm + test.gls + test.lme + test.lvmfit == 0)){
        index.wrong <- which(test.lm + test.gls + test.lme + test.lvmfit == 0)
        stop("Argument \'model\' must be a list of objects that inherits from lm/gls/lme/lvmfit. \n",
             "Incorrect element(s): ",paste(index.wrong, collapse = " "),".\n")
    }


    ## ** define the contrast matrix
    out <- list()
    if (is.character(linfct)){
        resC <- createContrast(model, par = linfct, add.variance = TRUE)
        contrast <- resC$contrast
        ls.contrast <- resC$ls.contrast
        if("rhs" %in% names(match.call()) == FALSE){
            rhs <- resC$null
        }
    }else if(is.matrix(linfct)){

        ls.contrast <- lapply(name.model, function(x){ ## x <- name.model[2]
            iColnames <- grep(paste0("^",x,": "), colnames(linfct), value = FALSE, fixed = FALSE)
            iRownames <- rowSums(linfct[,iColnames]!=0)>0
            linfct[iRownames, iColnames,drop=FALSE]            
        })
        names(ls.contrast) <- name.model
        contrast <- linfct
        if("rhs" %in% names(match.call()) == FALSE){ ## left rhs to default value
            rhs <- rep(0, NROW(contrast))
        }else if(length(rhs)!=NROW(contrast)){
            stop("mismatch between the dimensions of argument \'rhs\' and argument \'contrast\' \n")
        }
    }else{
        stop("Argument \'linfct\' must be a matrix or a vector of characters. \n",
             "Consider using  out <- createContrast(...) and pass out$contrast to linfct. \n")
    }

    ## ** check whether it is possible to compute df
    if(identical(df, TRUE)){
           
        ## does each model has the same df?
        ## test.df <- try(lapply(model, df.residual), silent = TRUE)
        ## if(inherits(test.df, "try-error")){
        ##     stop("Cannot check the degrees of freedom for each model - no \'df.residual\' method available \n",
        ##          "Consider setting the argument \'df\' to FALSE \n")
        ## }
        
        ## if(any(sapply(test.df,is.null))){
        ##     stop("Cannot compute residual degrees of freedom for all models \n",
        ##          "Consider setting the argument \'df\' to FALSE \n")
        ## }

        ## if(length(unique(unlist(test.df)))>1){
        ##     stop("Residual degrees of freedom differ across models \n",
        ##          "Consider setting the argument \'df\' to FALSE \n")
        ## }

        ## are linear hypothesis model specific?
        ls.testPerModel <- lapply(ls.contrast, function(iModel){
            rowSums(contrast[,colnames(iModel)]!=0)>0
        })
        
        if(any(Reduce("+",ls.testPerModel)>1)){
            stop("Cannot compute the degrees of freedom for tests performed across several models \n",
                 "Consider setting the argument \'df\' to FALSE \n")
        }    
    }
    
    ## ** Extract influence functions from all models    
    ls.res <- lapply(1:n.model, function(iM){ ## iM <- 1

### *** Pre-compute quantities
        if(!inherits(model[[iM]],"lm2") && !inherits(model[[iM]],"gls2") && !inherits(model[[iM]],"lme2") && !inherits(model[[iM]],"lvmfit2")){
            sCorrect(model[[iM]], df = df) <- bias.correct
        }
        out$param <- model[[iM]]$sCorrect$param
        name.param <- names(out$param)
        name.object.param <- paste0(name.model[iM],": ",name.param)
        out$param <- setNames(out$param, name.object.param)
        
### *** Compute df for each test
        if(df){
            ## here null does not matter since we only extract the degrees of freedom
            iContrast <- ls.contrast[[iM]]
            colnames(iContrast) <- name.param
        
            iWald <- compare2(model[[iM]], contrast = iContrast, as.lava = FALSE)
            out$df <- iWald[1:(NROW(iWald)-1),"df"]
        }else{
            out$df <- Inf
        }
### *** get iid decomposition
        index.missing <- model[[iM]]$na.action
        n.obs <- stats::nobs(model[[iM]]) + length(index.missing)

        out$iid <- matrix(NA, nrow = n.obs, ncol = length(name.param),
                          dimnames = list(NULL, name.param))
        
        ## dim(iid2(model[[iM]], robust = robust, cluster = cluster))
        ## dim(out$iid)
        ## names(model[[iM]])
        out$iid[setdiff(1:n.obs,index.missing),] <- iid2(model[[iM]], robust = robust, cluster = cluster)
        
        colnames(out$iid) <- name.object.param

        return(out)
        
    })
    seq.df <- unlist(lapply(ls.res,"[[","df"))
    seq.param <- unlist(lapply(ls.res,"[[","param"))

    if(df){
        if(length(unique(seq.df))>1){
            warning("Unequal degrees of freedom for the Wald statistics \n",
                    "The median of the degrees of freedom is used.")
        }
        df.global <- round(stats::median(seq.df), digits = 0)
    }else{
        df.global <- 0
    }
        
    ls.iid <- lapply(ls.res,"[[","iid")
    n.obs <- unique(unlist(lapply(ls.iid, NROW)))
    if(length(n.obs)>1){
        stop("Mismatch between the number of observations in the iid \n",
                "Likely to be due to the presence of missing values \n")
        
    }
    M.iid <- do.call(cbind,ls.iid)
    if(any(is.na(M.iid))){
        M.iid[is.na(M.iid)] <- 0
    }
    vcov.model <- crossprod(M.iid)
    
    
### ** convert to the appropriate format
    out <- list(model = model,
                linfct = linfct,
                rhs = unname(rhs),
                coef = seq.param,
                vcov = vcov.model,
                df = df.global,
                alternative = "two.sided",
                type = NULL,
                robust = robust,
                bias.correct = bias.correct)
    class(out) <- c("glht2","glht")
        
    ### ** export
    return(out)    
}


## * .calcClosure
.calcClosure <- function(name, estimate, covariance, type, df){

    n.hypo <- length(name)
    correlation <- stats::cov2cor(covariance)

    ## ** create all possible hypotheses
    ls.closure <- lapply(n.hypo:1, function(iNtest){ ## iNtest <- 1  
        iList <- list(M = utils::combn(name, m = iNtest))
        iList$vec <- apply(iList$M, 2, paste, collapse = ",")
        return(iList)
    })

    ## ** compute all p.values
    for(iLevel in 1:length(ls.closure)){ ## iLevel <- 1
        ls.closure[[iLevel]]$test <- t(apply(ls.closure[[iLevel]]$M, 2, function(iHypo){
            index <- which(name %in% iHypo)
            if(type == "chisq"){
                return(.ChisqTest(estimate[index], covariance = covariance[index,index,drop=FALSE], df = df))
            }else if(type == "max"){
                return(.tTest(estimate[index],
                              covariance = covariance[index,index,drop=FALSE],
                              correlation = correlation[index,index,drop=FALSE], df = df))
            }
        }))
        rownames(ls.closure[[iLevel]]$test) <- ls.closure[[iLevel]]$vec
    }
    
    ## ** find all hypotheses in the closure related to an individual hypothesis
    ls.hypo <- vector(mode = "list", length = n.hypo)
    for(iHypo in 1:n.hypo){ ## iHypo <- 1
        ls.hypo[[iHypo]] <- do.call(rbind,lapply(ls.closure, function(iClosure){ ## iClosure <- 1
            iIndex <- which(colSums(iClosure$M==name[iHypo])>0)
            data.frame(hypothesis = iClosure$vec[iIndex],
                       statistic = as.double(iClosure$test[iIndex,"statistic"]),
                       p.value = as.double(iClosure$test[iIndex,"p.value"]))
        }))
    }
    names(ls.hypo) <- name
        
    ## ** adjusted p.values
    vec.p.value <- unlist(lapply(ls.hypo, function(x){max(x$p.value)}))
    return(list(closure = ls.closure,
                test = ls.hypo,
                p.value = vec.p.value))
    
}

## * .tTest
.tTest <- function(estimate, covariance, correlation, df, ...){
    df1 <- length(estimate)
    statistic <- max(abs(estimate/sqrt(diag(covariance))))
    if(is.null(df)){
        distribution <-  "gaussian"
    }else{
        distribution <- "student"
    }
    p.value <- .calcPmaxIntegration(statistic, p = df1, Sigma = correlation, df = df,
                                    distribution = distribution)
    return(c("statistic" = statistic,
             "p.value" = p.value))
}

## * .ChisqTest
.ChisqTest <- function(estimate, covariance, df, ...){
    df1 <- length(estimate)
    ## q * statistic ~ chisq or fisher
    statistic <- as.double(matrix(estimate, nrow = 1) %*% solve(covariance) %*% matrix(estimate, ncol = 1)) / df1
    if(!is.null(df)){
        return(c("statistic" = statistic,
                 "p.value" = 1-stats::pf(statistic, df1 = df1, df2 = df)))
    }else{
        return(c("statistic" = statistic,
                 "p.value" = 1-stats::pchisq(statistic, df = df1)))
        
    }
}


 
