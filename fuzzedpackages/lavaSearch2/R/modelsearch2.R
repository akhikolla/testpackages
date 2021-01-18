## * modelsearch2 (documentation)
#' @title Data-driven Extension of a Latent Variable Model
#' @description Procedure adding relationship between variables that are supported by the data.
#' @name modelsearch2
#' 
#' @param object a \code{lvmfit} object.
#' @param link [character, optional for \code{lvmfit} objects] the name of the additional relationships to consider when expanding the model. Should be a vector containing strings like "Y~X". See the details section.
#' @param data [data.frame, optional] the dataset used to identify the model
#' @param method.p.adjust [character] the method used to adjust the p.values for multiple comparisons.
#' Can be any method that is valid for the \code{stats::p.adjust} function (e.g. \code{"fdr"}).
#' Can also be \code{"max"}, \code{"fastmax"}, or \code{"gof"}.
#' @param method.maxdist [character] the method used to estimate the distribution of the max statistic.
#' \code{"resampling"} resample the score under the null to estimate the null distribution.
#' \code{"bootstrap"} performs a wild bootstrap of the iid decomposition of the score to estimate the null distribution.
#' \code{"approximate"} attemps to identify the latent gaussian variable corresponding to each score statistic (that is chi-2 distributed).
#' It approximates the correlation matrix between these latent gaussian variables and uses numerical integration to compute the distribution of the max.
#' @param n.sample [integer, >0] number of samples used in the resampling approach.
#' @param na.omit should tests leading to NA for the test statistic be ignored. Otherwise this will stop the selection process.
#' @param alpha [numeric 0-1] the significance cutoff for the p-values.
#' When the p-value is below, the corresponding link will be added to the model
#' and the search will continue. Otherwise the search will stop.
#' @param nStep the maximum number of links that can be added to the model.
#' @param trace [logical] should the execution of the function be traced?
#' @param cpus the number of cpus that can be used for the computations.
#'
#' @details
#' method.p.adjust = \code{"max"} computes the p-values based on the distribution of the max statistic.
#' This max statistic is the max of the square root of the score statistic.
#' The p-value are computed integrating the multivariate normal distribution.
#' 
#' method.p.adjust = \code{"fastmax"} only compute the p-value for the largest statistic.
#' It is faster than \code{"max"} and lead to identical results.
#' 
#' method.p.adjust = \code{"gof"} keep adding links until the chi-squared test (of correct specification of the covariance matrix) is no longer significant.
#' @return A list containing:
#' \itemize{
#' \item sequenceTest: the sequence of test that has been performed.
#' \item sequenceModel: the sequence of models that has been obtained.
#' \item sequenceQuantile: the sequence of rejection threshold. Optional. 
#' \item sequenceIID: the influence functions relative to each test. Optional. 
#' \item sequenceSigma: the covariance matrix relative to each test. Optional. 
#' \item initialModel: the model before the sequential search.
#' \item statistic: the argument \code{statistic}.
#' \item method.p.adjust: the argument \code{method.p.adjust}.
#' \item alpha: [numeric 0-1] the significance cutoff for the p-values.
#' \item cv: whether the procedure has converged.
#' } 
#'
#' @concept modelsearch
#' @export
`modelsearch2` <-
    function(object, link, data,
             method.p.adjust, method.maxdist, n.sample, na.omit, 
             alpha,  nStep, trace, cpus) UseMethod("modelsearch2")


## * modelsearch2 (example)
#' @rdname modelsearch2
#' @examples
#'
#' ## simulate data 
#' mSim <- lvm()
#' regression(mSim) <- c(y1,y2,y3,y4)~u
#' regression(mSim) <- u~x1+x2
#' categorical(mSim,labels=c("A","B","C")) <- "x2"
#' latent(mSim) <- ~u
#' covariance(mSim) <- y1~y2
#' transform(mSim, Id~u) <- function(x){1:NROW(x)}
#'
#' set.seed(10)
#' df.data <- lava::sim(mSim, n = 1e2, latent = FALSE)
#' 
#' ## only identifiable extensions
#' m <- lvm(c(y1,y2,y3,y4)~u)
#' latent(m) <- ~u
#' addvar(m) <- ~x1+x2
#' 
#' e <- estimate(m, df.data)
#'
#' \dontrun{
#' resSearch <- modelsearch(e)
#' resSearch
#'
#' resSearch2 <- modelsearch2(e, nStep = 2)
#' resSearch2
#' }
#' \dontshow{
#' search.link <- c("u~x1","u~x2","y1~x1","y1~x2","y1~~y2","y1~~y3")
#' resSearch2 <- modelsearch2(e, nStep = 2, link = search.link)
#' resSearch2
#' }
#'
#' ## some extensions are not identifiable
#' m <- lvm(c(y1,y2,y3)~u)
#' latent(m) <- ~u
#' addvar(m) <- ~x1+x2 
#'
#' e <- estimate(m, df.data)
#'
#' \dontrun{
#' resSearch <- modelsearch(e)
#' resSearch
#' resSearch2 <- modelsearch2(e)
#' resSearch2
#' }
#'
#' ## for instance
#' mNI <- lvm(c(y1,y2,y3)~u)
#' latent(mNI) <- ~u
#' covariance(mNI) <- y1~y2
#' ## estimate(mNI, data = df.data)
#' ## does not converge
#'
#' 
#' 

## * modelsearch2.lvmfit (code)
#' @rdname modelsearch2
#' @export
modelsearch2.lvmfit <- function(object, link = NULL, data = NULL, 
                                method.p.adjust = "fastmax", method.maxdist = "approximate", n.sample = 1e5, na.omit = TRUE, 
                                alpha = 0.05, nStep = NULL, 
                                trace = TRUE, cpus = 1){

    ## ** check arguments
    ## object
    if(any(is.na(model.frame(object))) && method.p.adjust %in% c("max","fastmax")){
        warning("Missing values - the iid decomposition of the test statistics will only be computed on complete data \n")
    }
    
    ## methods
    method.p.adjust <- match.arg(method.p.adjust, c("fastmax", "max", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none","gof"))    
    if(method.p.adjust == "gof" ){
        method.p.adjust <- "none"
        stop.gof <- TRUE
    }else{
        stop.gof <- FALSE
    }

    if(n.sample<0 || (n.sample %% 1 != 0) ){
        stop("Argument \'n.sample\' must be a positive integer \n")
    }
    method.maxdist <- match.arg(method.maxdist, c("approximate","resampling","bootstrap"))    

    ## cpus
    if(is.null(cpus)){ cpus <- parallel::detectCores()}

        if(is.null(cpus) || cpus > 1){
        test.package <- try(requireNamespace("foreach"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'foreach\' \n",
                 "This package is necessary when argument \'cpus\' is greater than 1 \n")
        }
    }
    
    if(!is.null(cpus) && cpus>1){
        if(cpus > parallel::detectCores()){
            stop("Argument \'cpus\' is greater than the number of available CPU cores \n",
                 "available CPU cores: ",parallel::detectCores(),"\n")
        }
    }

    ## ## extra arguments 
    ## dots <- list(...)
    ## if(length(dots)>0){
    ##     stop("modelsearch2 does not take any extra arguments \n",
    ##          "name of the extra arguments: \"",paste(names(dots), collapse = "\" \""),"\" \n")
    ## }

    ## ** prepare
    ## *** data
    if(is.null(data)){
        data <- as.data.frame(stats::model.frame(object, all = TRUE))
    }
    
    ## *** normalize the links
    if(is.null(link)){
        res.find <- do.call(findNewLink,
                            args = c(list(object$model,
                                          data = data,
                                          output = "names")))
        directive <- res.find$directional
        restricted <- res.find$M.links
        link <- res.find$link
        if(is.null(link)){
            stop("Automatic search has not found any possible additional link \n",
                 "Consider specifying manually the argument \'link\' \n")
        }
    }else{
        resLink <- .initializeLinks(object, data = data, link = link)
        object <- resLink$object
        link <- resLink$link
        directive <- resLink$directive
        restricted <- resLink$restricted
    }    
    
    ## ** initialization
    if(is.null(nStep)){
        nStep <- NROW(restricted)
    }
    iStep <- 1
    iRestricted <- restricted
    iDirective <- directive
    iLink <- link
    iObject <- object

    ## update of the model
    add.args <- setdiff(names(object$call), c("","object","data","control"))
    ls.call <- lapply(add.args, function(arg){object$call[[arg]]})
    names(ls.call) <- add.args

    ls.call$data <- data
    if(!is.null(data)){
        index.cols <- which(names(data)%in%names(ls.call$data)==FALSE)
        if(length(index.cols)>0){
            ls.call$data <- cbind(ls.call$data,
                                  subset(as.data.frame(data), select = index.cols))
        }
    }
    if(is.null(object$control)){
        ls.call$control <- list()
    }else{
        ls.call$control <- object$control
    }
    ls.call$control$trace <- FALSE
    
    ## output
    ls.seqTests <- list()
    ls.seqModels <- list()
    ls.seqIID <- list() # only for method.p.adjust = "max"
    ls.seqSigma <- list() # only for method.p.adjust = "max"
    vec.seqQuantile <- NULL # only for method.p.adjust = "max"
    
    ## criterion    
    cv <- FALSE

    ## define cluster
    if(cpus>1){
        if(trace>0){
            cl <- parallel::makeCluster(cpus, outfile = "")
        }else{
            cl <- parallel::makeCluster(cpus)
        }
        doParallel::registerDoParallel(cl)
        
        vec.packages <- c("lavaSearch2","lava")
        parallel::clusterCall(cl, fun = function(x){
            sapply(vec.packages, function(iP){
                suppressPackageStartupMessages(attachNamespace(iP)) ## requireNamespace did not worked
            })
        })
        
    }else{
        cl <- NULL
    }

    ## ** display a summary of the call
    if(trace>0){

        cat("\n",
            "** Sequential variable selection using the score statistic ** \n",
            " Number of possible additional links         : ",length(link)," \n",
            " Maximum number of steps                     : ",nStep,"\n",        
            " Adjustment method for multiple comparisons  : ",method.p.adjust,"\n",
            " Confidence level                            : ",1-alpha,"\n",
            " Number of cpus                              : ",cpus,"\n\n",
            sep="")
    }

    ## ** Forward search
    if(stop.gof){
        if(trace>0){
            cat("p.Chi-squared test = ",gof(iObject)$fit$p.value,"\n", sep = "")
        }
        if(gof(iObject)$fit$p.value >= alpha){
            cv <- TRUE
        }
    }
        
    while(iStep <= nStep && NROW(iRestricted)>0 && cv==FALSE){
        if(trace >= 1){cat("Step ",iStep,":\n",sep="")}

        
        
        resStep <- .oneStep_scoresearch(iObject, data = data,
                                        restricted = iRestricted, link = iLink, directive = iDirective,
                                        method.p.adjust = method.p.adjust, method.maxdist = method.maxdist, n.sample = n.sample,
                                        cl = cl, trace = trace)

        ## ** update according the most significant p.value
        ## *** check convergence
        if(stop.gof){
            cv <- FALSE
            test.na <- FALSE
        }else if(na.omit || method.p.adjust == "fastmax"){
            cv <- all(stats::na.omit(resStep$test$adjusted.p.value) > alpha)
            test.na <- FALSE
        }else{
            cv <- all(resStep$test$adjusted.p.value > alpha)
            if(is.na(cv)){                    
                cv <- TRUE
                test.na <- TRUE
            }else{
                test.na <- FALSE
            }
        }

        ## *** identify most promising test
        index.maxTest <- which.max(abs(resStep$test$statistic))[1]
        resStep$test$selected <- FALSE
        resStep$test[index.maxTest,"selected"] <- (resStep$test[index.maxTest,"adjusted.p.value"] <= alpha)
        resStep$test$nTests <- NROW(resStep$test)
        resStep$test <- resStep$test[order(resStep$test$statistic),]

        ## *** update the model
        if(cv==FALSE){
            ls.call$x <- addLink(iObject$model, var1 = iRestricted[index.maxTest,1], var2 = iRestricted[index.maxTest,2],
                                 covariance = 1-iDirective[index.maxTest])

            ## first attempt
            ls.call$control$start <- stats::coef(iObject)
            suppressWarnings(
                iObject <- tryCatch(do.call(lava::estimate, args = ls.call),
                                    error = function(x){x},
                                    finally = function(x){x})
            )

            ## second attempt
            if(inherits(iObject,"error") || iObject$opt$convergence>0){
                ls.call$control$start <- NULL
                suppressWarnings(
                    iObject <- do.call(lava::estimate, args = ls.call)
                )
                if(inherits(iObject,"error") || iObject$opt$convergence>0){
                    stop("Estimation of the extended latent variable model did not converge \n")
                }
            }

            ## update links
            iLink <- iLink[-index.maxTest]
            iRestricted <- iRestricted[-index.maxTest,,drop=FALSE]
            iDirective <- iDirective[-index.maxTest]
        }
        


        ## *** update the output
        ls.seqTests[[iStep]] <- resStep$test
        ls.seqModels[[iStep]] <- iObject
        if(method.p.adjust %in% c("max","fastmax")){
            ls.seqIID[[iStep]] <- resStep$iid
            ls.seqSigma[[iStep]] <- resStep$Sigma
            vec.seqQuantile <- c(vec.seqQuantile,resStep$test$quantile[1])
        }


        ## *** display results
        if(trace > 0){
            rowSelected <- NROW(resStep$test)
            if(cv==FALSE){
                cat("add ",as.character(resStep$test[rowSelected, "link"]),
                    " (statistic = ",resStep$test[rowSelected, "statistic"],
                    ", adjusted.p.value = ",resStep$test[rowSelected, "adjusted.p.value"],
                    ")\n",sep="")
            }else{
                if(test.na){
                    cat("NA among the test statistics \n")
                }else{
                    cat("no variable to add",
                        " (statistic = ",resStep$test[rowSelected, "statistic"],
                        ", adjusted.p.value = ",resStep$test[rowSelected, "adjusted.p.value"],
                        ")\n",sep="")
                }
            }
        }

        ## *** check convergence gof
        if(stop.gof){
            if(trace>0){
                cat("p.Chi-squared test = ",gof(iObject)$fit$p.value,"\n", sep = "")
            }
            if(gof(iObject)$fit$p.value >= alpha){
                cv <- TRUE
            }
        }
    

        
        iStep <- iStep + 1
    }

    if(cpus>1){
        parallel::stopCluster(cl)
    }

    ## ** export
    if(length(ls.seqIID)==0){ls.seqIID <- NULL}
    if(length(ls.seqSigma)==0){ls.seqSigma <- NULL}
    output <- list(sequenceTest = ls.seqTests,
                   sequenceModel = ls.seqModels,
                   sequenceQuantile = vec.seqQuantile,
                   sequenceIID = ls.seqIID,
                   sequenceSigma = ls.seqSigma,
                   initialModel = object,
                   method.p.adjust = method.p.adjust,
                   alpha = alpha,
                   cv = cv)
    class(output) <- "modelsearch2"
    return(output)
}

## * .initializeLinks
.initializeLinks <- function(object, data, link){
    restricted <- do.call(cbind,initVarLinks(link))
    directive <- rep(TRUE, length(link))

    ## ** identify covariance link
    index.Ndir <- grep(lava.options()$symbols[2],link,fixed=TRUE)
    if(length(index.Ndir)>0){
        directive[index.Ndir] <- FALSE
    }
   
    ## ** get all vars
    if(is.null(data)){            
        data <- evalInParentEnv(object$call$data)             
        if(is.null(data)){
            data <- lava::sim(object,1)
        }            
    }
    ## ** take care of categorical variables
    iData <- try(eval(object$call$d), silent = TRUE)
    if(!inherits(iData, "data.frame")){
        iData <- evalInParentEnv(object$call$data)
        if(!inherits(iData, "data.frame")){
            stop("Could not identify argument data in object$call \n")
        }
    }

    M.linkvar <- do.call(rbind,lapply(1:NROW(restricted), function(row){ ## row <- 2
        data.frame(Y = unname(restricted[row,1]),
                   X = unname(var2dummy(object$model,
                                        data = iData,
                                        var = restricted[row,2])),
                   dir = directive[row],
                   stringsAsFactors = FALSE)
    }))

    ## check no missing covariance links
    index.regression <- which(M.linkvar[,"dir"]==TRUE)
    if(length(index.regression)>0){
        index.covariance <- which(M.linkvar[index.regression,"X"] %in% lava::endogenous(object) + M.linkvar[index.regression,"X"] %in% lava::endogenous(object) ==2)
        if(length(index.covariance)>0){
            stop("Covariance links must be indicated with the symbol \"",lava.options()$symbols[2],"\" \n",
                 "Possible covariance links: ",paste0(link[index.regression][index.covariance], collapse = " "),"\n")
        }
    }
    restricted2 <- as.matrix(M.linkvar[,1:2,drop=FALSE])
    directive <- M.linkvar[,3]
    link <- paste0(restricted2[,1],lava.options()$symbols[2-directive],restricted2[,2])

    ## ** check links
    allVars.link <- setdiff(unique(as.vector(restricted2)), lava::latent(object$model))
    allVars.model <- lava::vars(object$model)
    allVars.data <- names(data)

    if(any(allVars.link %in% allVars.model == FALSE)){
        missing.var <- allVars.link[allVars.link %in% allVars.model == FALSE]
        if(any(allVars.link %in% allVars.data == FALSE)){
            missing.var <- allVars.link[allVars.link %in% allVars.data == FALSE]
            stop("Some links contains variables that are not in the latent variable model \n",
                 "variables(s) : \"",paste(missing.var,collapse ="\" \""),"\"\n")
        }

    }
    
    ## ** check covariance links
    if(any(directive==FALSE)){
        if(any(restricted2[directive==FALSE,1] %in% lava::exogenous(object)) || any(restricted2[directive==FALSE,2] %in% lava::exogenous(object))){
            wrong <- union(which(restricted2[directive==FALSE,1] %in% exogenous(object)),
                           which(restricted2[directive==FALSE,2] %in% exogenous(object)))
            stop("Covariance links can only relate endogenous variables \n",
                 "Covariance link(s) involving exogenous variables: ,", paste(link[wrong], collapse = " ; "),"\n")
        }
    }
    
     return(list(object = object,
                link = link,
                directive = directive,
                restricted = restricted2))
}
## * .oneStep_scoresearch
.oneStep_scoresearch  <- function(object, data,
                                  restricted, link, directive,
                                  method.p.adjust, alpha, method.maxdist, n.sample,
                                  cl, trace){

    ## ** initialization
    n.link <- NROW(restricted)
    coef.object <- coef(object)
    namecoef.object <- names(coef.object)
    ncoef.object <- length(coef.object)
    type.information <- lava.options()$search.type.information
    type.statistic <- lava.options()$search.sample.stat
    
    ## ** warper
    warper <- function(iterI){ # iterI <- 1

        out <- list(table = data.frame(statistic = as.numeric(NA),
                                       p.value = as.numeric(NA),
                                       adjusted.p.value = as.numeric(NA),
                                       dp.Info = as.numeric(NA),
                                       stringsAsFactors = FALSE),
                    iid = NULL)

        ## *** define extended model
        newModel <- addLink(object$model, var1 = restricted[iterI,1], var2 = restricted[iterI,2],
                            covariance = 1-directive[iterI])

        ## *** remove useless variables
        Mlink <- newModel$M + (newModel$cov - diag(1, NROW(newModel$cov), NCOL(newModel$cov)))
        noLink.var <- names(which((rowSums(Mlink)==0)+(colSums(Mlink)==0)==2))
        if(length(noLink.var)>0){
            rmvar(newModel) <- noLink.var
        }

        ## *** compute sufficient statistics
        ## necessary otherwise information can have a weird behavior, e.g.
        ## library(lava)
        ##
        ## mSim <- lvm(Y ~ X1, X2 ~ eta, X3 ~ eta, X4 ~ eta)
        ## latent(mSim) <- ~eta
        ## d <- sim(mSim, 100)
        ##
        ## m <- lvm(Y~X1+X2) 
        ## e <- estimate(m, d)
        ## information(e, data = d, p = coef(e)) ## gold standard
        ## information(m, data = d, p = coef(e)) ## issue
        ##
        ## fix
        ## ss <- lava:::procdata.lvm(m, data = d, missing = FALSE) 
        ## mm <- lava::fixsome(m, measurement.fix=TRUE, S=ss$S, mu=ss$mu, n = ss$n, debug=FALSE)
        ## information(mm, data = d, p = coef(e)) ## ok
        suffStat <- lava_procdata.lvm(newModel, data = data, missing = inherits(object,"lvm.missing")) 
        newModel <- lava::fixsome(newModel, measurement.fix=TRUE, S=suffStat$S, mu=suffStat$mu, n = suffStat$n, debug=FALSE)

        ## *** define constrained coefficients
        coef0.new <- setNames(rep(0, ncoef.object+1), coef(newModel))
        coef0.new[namecoef.object] <- coef.object

        ## *** compute the iid decomposition and statistic
        namecoef.newobject <- names(coef0.new)
        Info <- lava::information(newModel, p = coef0.new, n = NROW(data), type = type.information, data = data)
        dimnames(Info) <- list(namecoef.newobject,namecoef.newobject)
        
        if(method.p.adjust %in% c("max","fastmax")){
            iid.score <- lava::score(newModel, p = coef0.new, data = data, indiv = TRUE)
            ## rm na
            iid.score <- iid.score[rowSums(is.na(iid.score))==0,]
            if(method.maxdist == "approximate"){
                ## compute decomposition
                out$iid <-  iid.score %*% solve(Info) %*% cbind(colSums(iid.score))
                colnames(out$iid) <- link[iterI]
                ## out$iid <- out$iid/sqrt(sum(out$iid^2))
                out$table$statistic <- sum(out$iid)
                out$table$dp.Info <- TRUE
                
            }else if(method.maxdist %in% c("resampling","bootstrap")){
                n.sample <- NROW(iid.score)

                InfoM12 <- matrixPower(Info, power  = -1/2, symmetric = TRUE, tol = 1e-15, print.warning = FALSE)
                out$table$dp.Info <- !("warning" %in% names(attributes(InfoM12)))
                dimnames(InfoM12) <- list(namecoef.newobject,namecoef.newobject)
                
                ## initial version
                ## linComb <- cbind(1, -solve(Info[link[iterI],link[iterI],drop=FALSE]) %*% Info[link[iterI],namecoef.object,drop=FALSE]) %*% Info[c(link[iterI],namecoef.object),c(link[iterI],namecoef.object)]
                if(FALSE){
                    InfoM1 <- crossprod(InfoM12)
                    dimnames(InfoM1) <- list(namecoef.newobject,namecoef.newobject)
                    linComb <- cbind(1, -Info[link[iterI],namecoef.object,drop=FALSE] %*% solve(Info[namecoef.object,namecoef.object,drop=FALSE])) %*% Info[c(link[iterI],namecoef.object),c(link[iterI],namecoef.object)]
                    iid.theta <- iid.score %*% InfoM1
                    colnames(iid.theta) <- namecoef.newobject
                    out$iid <- iid.theta[,link[iterI]] %*% linComb[,namecoef.newobject,drop=FALSE] %*% InfoM12
                    colnames(out$iid) <- paste0(link[iterI],":",namecoef.newobject)
                }
                ## short version
                out$iid <- (iid.score[,link[iterI],drop=FALSE] - iid.score[,namecoef.object,drop=FALSE] %*% solve(Info[namecoef.object,namecoef.object,drop=FALSE]) %*% Info[namecoef.object,link[iterI],drop=FALSE]) %*% InfoM12[link[iterI],,drop=FALSE]
                colnames(out$iid) <- paste0(link[iterI],":",namecoef.newobject)

                out$table$statistic <- as.double(crossprod(colSums(out$iid))) ## first order approximation (almost identical to exact value)
            }
        }else{
            ## ee.lvm <- estimate(newModel, data = data)
            ## SS <- score(ee.lvm, p = coef0.new)
            ## II <- information(ee.lvm, p = coef0.new)
            ## SS %*% solve(I) %*% t(SS)
            score <- lava::score(newModel, p = coef0.new, indiv = FALSE, data = data)
            
            out$table$statistic <- as.double(score %*% solve(Info) %*% t(score))
            ## range(Info - II)
            ## range(score - SS)
        }
        return(out)
    }

    ## ** compute score tests
    if(trace>0){
        cat(" - compute score test for all possible additional links \n")
    }
    
    if(!is.null(cl)){
    
        if(trace > 0){
            pb <- utils::txtProgressBar(max = n.link, style = 3) 
        }

        
        ## get influence function
        i <- NULL # [:for CRAN check] foreach
        res <- foreach::`%dopar%`(
                            foreach::foreach(i = 1:n.link,
                                             .export = c("lava_procdata.lvm"),
                                             .combine = function(res1,res2){
                                                 res <- list(table = rbind(res1$table,res2$table),
                                                             iid = cbind(res1$iid,res2$iid))
                                                 return(res)
                                             }), {
                                                 if(trace>0){utils::setTxtProgressBar(pb, i)}
                                                 return(warper(i))
                                             })

        if(trace>0){close(pb)}
        
    }else{
        
        if(trace>0){
            resApply <- pbapply::pblapply(1:n.link, warper)            
        }else{
            resApply <- lapply(1:n.link, warper)
        }
        res <- list(table = do.call(rbind, lapply(resApply,"[[","table")),
                    iid = do.call(cbind,lapply(resApply,"[[","iid")))
        
    }
    ## index.iid <- unlist(lapply(1:n.link, function(iL){ ## iL <- 1
    ## rep(iL,times = NCOL(res$iid[[iL]]))
    ## }))
    table.test <- data.frame(link = link, res$table, error = NA, stringsAsFactors = FALSE)
    iid.link <- res$iid

    ## ** p.value
    statistic <- as.numeric(table.test[,"statistic"])
    if(any(statistic<0)){
        stop("Negative score statistic \n")
    }
    ## univariate rejection area
    table.test[,"p.value"] <- 1-stats::pchisq(statistic, df = 1)

    ## ** adjusted p.value
    if(method.p.adjust %in% c("fastmax","max")){
        if(method.maxdist == "approximate"){
            outDistMax <- .approxMaxDistChi2(table = table.test, iid = iid.link, statistic = statistic, method.p.adjust = method.p.adjust,
                                            link = link, n.link = n.link,
                                            search.calc.quantile.int = lava.options()$search.calc.quantile.int, alpha = alpha,
                                            cl = cl, trace = trace)
        }else if(method.maxdist %in% c("resampling","bootstrap")){
            outDistMax <- .sampleMaxDistChi2(table = table.test, iid = iid.link, statistic = statistic, method.p.adjust = method.p.adjust,
                                        link = link, n.link = n.link, n.sample = n.sample, method = method.maxdist,
                                        cl = cl, trace = trace)
        }

        table.test <- outDistMax$table
        Sigma <- outDistMax$Sigma

    }else{        
        table.test[, "adjusted.p.value"] <- stats::p.adjust(table.test$p.value, method = method.p.adjust)
        table.test[, "quantile"] <- as.numeric(NA)
        Sigma <- NULL        
    }


    return(list(test = table.test,
                Sigma = Sigma,
                iid = iid.link))
}


## * approxMaxDistChi2
.approxMaxDistChi2 <- function(table, iid, statistic, method.p.adjust,
                               link, n.link,
                               search.calc.quantile.int, alpha,
                               cl, trace){
        Sigma <- stats::cor(iid)
        dimnames(Sigma) <- list(link,link)
        if(method.p.adjust == "fastmax"){
            index.maxstat <- which.max(statistic)
        
            resInt <- .calcPmaxIntegration(statistic = sqrt(statistic[index.maxstat]), p = n.link,
                                           Sigma = Sigma, df = NULL, distribution = "gaussian")
        
            table[index.maxstat, "adjusted.p.value"] <- as.double(resInt)
            table[index.maxstat, "error"] <- attr(resInt,"error")
        }else if(method.p.adjust == "max"){

            if(is.null(cl)){
                if(trace>0){      
                    ls.resInt <- pbapply::pblapply(1:n.link, function(i){
                        .calcPmaxIntegration(statistic = sqrt(statistic[i]), p = n.link,
                                             Sigma = Sigma, df = NULL, distribution = "gaussian")
                    })            
                }else{
                    ls.resInt <- lapply(1:n.link, function(i){
                        .calcPmaxIntegration(statistic = sqrt(statistic[i]), p = n.link,
                                             Sigma = Sigma, df = NULL, distribution = "gaussian")
                    })            
                }
            }else{
            
                if(trace>0){
                    pb <- utils::txtProgressBar(max = n.link, style = 3)                   
                }

                ## export package
                parallel::clusterCall(cl, fun = function(x){
                    suppressPackageStartupMessages(requireNamespace("mvtnorm", quietly = TRUE))
                })
        
                value <- NULL # [:for CRAN check] foreach
                ls.resInt <- foreach::`%dopar%`(
                                          foreach::foreach(value = 1:n.link,
                                                           .export = c(".calcPmaxIntegration")),
                                          {
                                              if(trace>0){utils::setTxtProgressBar(pb, value)}
                                              return(.calcPmaxIntegration(statistic = sqrt(statistic[value]), p = n.link,
                                                                          Sigma = Sigma, df = NULL, distribution = "gaussian"))
                                          })

                if(trace>0){close(pb)}
            }
            names(ls.resInt) <- link
        
            table[, "adjusted.p.value"] <- unlist(lapply(ls.resInt,as.double))
            table[, "error"] <- unlist(lapply(ls.resInt,attr,"error"))
        }
        
        if(lava.options()$search.calc.quantile.int){
            table[, "quantile"] <- .calcQmaxIntegration(alpha = alpha, p = n.link,
                                                             Sigma = Sigma,
                                                             df = NULL, distribution = "gaussian")
        }

        return(list(table = table,
                    Sigma = Sigma))
}
## * sampleMaxDistChi2
.sampleMaxDistChi2 <- function(table, iid, statistic, method.p.adjust,
                               link, n.link, n.sample, method,
                               cl, trace){
    p <- NCOL(iid)
    n <- NROW(iid)
    ls.name <- strsplit(colnames(iid), split = ":")
    vec.model <- unlist(lapply(ls.name,"[[",1))
    Umodel <- unique(vec.model)
    ls.indexModel <- tapply(1:length(vec.model),vec.model,list)

    ## ** sampling under H0

    ## *** resampling    
    if(method == "resampling"){
        Sigma <- crossprod(iid)
        sample2 <- mvtnorm::rmvnorm(n.sample, mean = rep(0,p), sigma = Sigma)^2
        M.scoreStat <- do.call(cbind,lapply(1:n.link, function(iModel){ ## iModel <- 1
            return(rowSums(sample2[,ls.indexModel[[iModel]],drop=FALSE]))
        }))
    }

    ## *** wild bootstrap
    if(method == "bootstrap"){
        Sigma <- NULL
        M.scoreStat <- wildBoot_cpp(iid = iid,
                                    lsIndexModel = lapply(ls.indexModel, function(x){x-1}),
                                    nSample = n.sample,
                                    nObs = n,
                                    nModel = n.link,
                                    p = p)
    }
    ## *** check
    ## apply(M.scoreStat,2, function(x){1-mean(x <= qchisq(0.95, df = 1))})
    ## hist(M.scoreStat[,1], freq = FALSE)
    ## points(seq(0,15,0.1), dchisq(seq(0,15,0.1), df = 1), col = "red", type = "l")
    
    
    ## ** p-value for each statistic
    p.value <- colMeans(sweep(M.scoreStat, MARGIN = 2, FUN = ">", STATS = statistic)) + 1/2 * colMeans(sweep(M.scoreStat, MARGIN = 2, FUN = "==", STATS = statistic))
    ## mean(M.scoreStat[,1] > statistic[1])

    ## ** p-value for the max statistic
    maxScoreStat <- apply(M.scoreStat,1,max)
    p.value.max <- sapply(statistic, function(iT){mean( maxScoreStat>iT + 0.5*(maxScoreStat==iT))})
    ## p.value.max / p.value
    table[, "p.value"] <- p.value
    table[, "adjusted.p.value"] <- p.value.max
    table[, "error"] <- NA
    return(list(table = table,
                Sigma = Sigma))

}



## * iidConstrainscore (obsolete)
## iidConstrainScore <- function(object, newobject){
##     if(all(sapply(newobject$mean,function(x){all(is.na(x))}))){
##         suffStat <- lava:::procdata.lvm(newobject, data = object$data$model.frame, missing = FALSE) 
##         newobject <- lava::fixsome(newobject, measurement.fix=TRUE, S=suffStat$S, mu=suffStat$mu, n = suffStat$n, debug=FALSE)
##     }
    
##     param0 <- coef(object)
##     name.param0 <- names(param0)
##     name.param <- coef(newobject)
##     n.param <- length(name.param)

##     newparam <- setNames(rep(0,n.param), name.param)
##     newparam[name.param0] <- param0
##     extraparam <- setdiff(name.param, name.param0)

##     S <- score(newobject, p = newparam, data = object$data$model.frame, indiv = TRUE)
##     I <- information(newobject, p = newparam, data = object$data$model.frame)
##     dimnames(I) <- list(name.param,name.param)
##     iInfo <- solve(I)
    
##     linComb <- cbind(1, solve(I[extraparam,extraparam,drop=FALSE]) %*% I[extraparam,name.param0,drop=FALSE]) %*% I[c(extraparam,name.param0),c(extraparam,name.param0)]
##     iid <- S %*% iInfo
##     out <- sweep(iid, MARGIN = 2, FUN = "*", STATS = as.double(linComb[,name.param]))
##     return(out/sqrt(object$data$n))
## }

