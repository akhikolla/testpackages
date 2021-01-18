### calibrateType1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (10:23) 
## Version: 
## Last-Updated: mar 13 2019 (11:55) 
##           By: Brice Ozenne
##     Update #: 813
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * Documentation - calibrateType1
##' @title Simulation Study Assessing Bias and Type 1 Error
##' @description Perform a simulation study over one or several sample size
##' to assess the bias of the estimate
##' and the type 1 error of the Wald test and robust Wald test
##' @name calibrateType1
##' 
##' @param object a \code{lvm} object defining the model to be fitted.
##' @param warmup [list of lvm] a list of \code{lvm} objects that will be sequentially fitted with for
##' starting values the parameter of the previous model in the list (if any). The parameters of the final
##' model of the list are used to initialize the fit of the model of interest (i.e. object).
##' @param param [character vector] names of the coefficient whose value will be tested. 
##' @param null [numeric vector] vector of null hypotheses, one for each model coefficient.
##' By default a vector of 0. 
##' @param n [integer vector, >0] sample size(s) considered in the simulation study.
##' @param n.rep [integer, >0] number of simulations per sample size.
##' @param cluster  [integer vector] the grouping variable relative to which the observations are iid.
##' Will be passed to \code{lava::estimate}.
##' @param correction [logical] should the type 1 error after correction be computed?
##' @param generative.object [lvm] object defining the statistical model generating the data.
##' @param generative.coef [name numeric vector] values for the parameters of the generative model.
##' Can also be \code{NULL}: in such a case the coefficients are set to default values decided by lava (usually 0 or 1).
##' @param true.coef [name numeric vector] expected values for the parameters of the fitted model.
##' @param n.true [integer, >0] sample size at which the estimated coefficients will be a reliable approximation of the true coefficients.
##' @param round.true [integer, >0] the number of decimal places to be used for the true value of the coefficients. No rounding is done if \code{NULL}.
##' @param bootstrap [logical] should bootstrap resampling be performed?
##' @param n.bootstrap [integer, >0] the number of bootstrap sample to be used for each bootstrap.
##' @param checkType1 [logical] returns an error if the coefficients associated to the null hypotheses do not equal 0.
##' @param checkType2 [logical] returns an error if the coefficients associated to the null hypotheses equal 0.
##' @param dir.save [character] path to the directory were the results should be exported.
##' Can also be \code{NULL}: in such a case the results are not exported.
##' @param F.test [logical] should a multivariate Wald test be perform testing simultaneously all the null hypotheses?
##' @param label.file [character] element to include in the file name.
##' @param seed [integer, >0] seed value that will be set at the beginning of the simulation to enable eproducibility of the results.
##' Can also be \code{NULL}: in such a case no seed is set.
##' @param cpus [integer >0] the number of processors to use.
##' If greater than 1, the simulations are performed in parallel. 
##' @param trace [integer] should the execution of the function be trace. Can be 0, 1 or 2.
##' @param ... [internal] Only used by the generic method.
##' 
##' @return An object of class \code{calibrateType1}.
##' @seealso \code{link{autoplot.calibrateType1}} for a graphical display of the bias or of the type 1 error.
##' 
##' @author Brice Ozenne
##'
##' @examples
##' \dontrun{
##' #### simulate data ####
##' m.Sim <- lvm(c(Y1[mu1:sigma]~1*eta,
##'                Y2[mu2:sigma]~1*eta,
##'                Y3[mu3:sigma]~1*eta,
##'                eta~beta1*Group+beta2*Gender))
##' latent(m.Sim) <- ~eta
##' categorical(m.Sim, labels = c("M","F")) <- ~Gender
##'
##' d <- lava::sim(m.Sim, 1e2)
##'
##' #### calibrate type 1 error on the estimated model ####
##' m <- lvm(Y1~eta,
##'          Y2~eta,
##'          Y3~eta,
##'          eta~Group+Gender)
##' e <- lava::estimate(m, data = d)
##' res <- calibrateType1(e, param = "eta~Group", n.rep = 100)
##' res <- calibrateType1(e, param = c("eta~Group","Y1~eta"), F.test = TRUE, n.rep = 100)
##' res <- calibrateType1(e, param = "eta~Group", n.rep = 100, cpus = 4)
##' summary(res)
##' }
##' 
##' @export
`calibrateType1` <-
    function(object, param, n.rep, ...) UseMethod("calibrateType1")

## * calibrateType1.lvm
##' @rdname calibrateType1
##' @export
calibrateType1.lvm <- function(object, param, n.rep, n, correction = TRUE, warmup = NULL, null = NULL,
                               F.test = FALSE, cluster = NULL,
                               generative.object = NULL, generative.coef = NULL, 
                               true.coef = NULL, n.true = 1e6, round.true = 2,              
                               bootstrap = FALSE, n.bootstrap = 1e3,
                               checkType1 = FALSE, checkType2 = FALSE,
                               dir.save = NULL, label.file = NULL,             
                               seed = NULL, cpus = 1, trace = 2, ...){

### ** test
    if(cpus>1){
        test.package <- try(requireNamespace("foreach"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'foreach\' \n",
                 "This package is necessary when argument \'cpus\' is greater than 1 \n")
        }

        if(cpus > parallel::detectCores()){
            stop("Argument \'cpus\' is greater than the number of available CPU cores \n",
                 "available CPU cores: ",parallel::detectCores(),"\n")
        }
    }
    
### ** prepare
    n.n <- length(n)

    ## *** generative model
    if(is.null(generative.object)){
        generative.object <- object
    }
    if(any(lava::manifest(object) %in% lava::manifest(generative.object) == FALSE)){
        missingVar <- lava::manifest(object)[lava::manifest(object) %in% lava::manifest(generative.object) == FALSE]
        stop("The object contains manifest variables that are not in the generative model \n",
             "missing manifest variables: \"",paste0(missingVar, collapse = "\" \""),"\" \n")
    }
    df.type_generative <- coefType(generative.object, as.lava = FALSE)
    name.param_generative <- df.type_generative[!is.na(df.type_generative$lava),"param"]
    n.param_generative <- length(name.param_generative)

    ## *** coef generative
    if(!is.null(generative.coef) && any(names(generative.coef) %in% name.param_generative == FALSE)){
        extraParam <- names(generative.coef)[names(generative.coef) %in% name.param_generative == FALSE]
        stop("Invalid argument \'generative.coef\': some of the coefficient names do not match those of the generative model \n",
             "extra coefficients: \"",paste0(extraParam, collapse = "\" \""),"\"\n")
    }

    ## *** coef of the fitted model
    if(is.null(true.coef)){
        if(trace>1){
            cat("  Estimate true coefficients using a sample size of n=",n.true," ", sep="")
        }
        e.true <- lava::estimate(object,
                                 data = lava::sim(generative.object, n = n.true, p = generative.coef, latent = FALSE))
        coef.true <- coef(e.true)
        if(!is.null(round.true)){
            coef.true <- round(coef.true, digits = round.true)
        }
        if(trace>1){
            cat("- done \n")
        }
    }else{
        if(trace>1){
            cat("  Check true coefficients ")
        }
        n.true <- n[n.n]
        e.true <- suppressWarnings(lava::estimate(object, cluster = cluster,
                                                  data = lava::sim(generative.object, n = n.true, p = generative.coef, latent = FALSE)))
        name.test <- names(coef(e.true))
        if(e.true$opt$convergence>0 && trace>1){
            cat("(incorrect convergence of the model) ")
        }
        
        if(!identical(sort(name.test),sort(names(true.coef)))){
            extraNames <- setdiff(names(true.coef),name.test)
            missingNames <- setdiff(name.test,names(true.coef))
            stop("Names of the coefficients in argument \'true.coef\' does not matches those of the estimated coefficients \n",
                 "missing names: \"",paste0(missingNames, collapse = "\" \""),"\" \n",
                 "extra names: \"",paste0(extraNames, collapse = "\" \""),"\" \n")
        }
        
        coef.true <- true.coef
        if(trace>1){
            cat("- done \n")
        }
        
    }
    name.coef <- names(coef.true)
    n.coef <- length(name.coef)
       
    ## *** type of the coef of the fitted model
    df.type <- coefType(e.true, as.lava = FALSE)
    df.type <- df.type[df.type$name %in% name.coef,]
    type.coef <- setNames(df.type$detail, df.type$name)

    ## *** null hypothesis
    n.param <- length(param)
    if(any(param %in% name.coef == FALSE)){
        incorrect.name <- param[param %in% name.coef == FALSE]
        possible.name <- setdiff(name.coef, param)
        ls.name <- lapply(incorrect.name, function(iN){
            dist.tempo <- utils::adist(x = iN, y = possible.name)
            return(possible.name[which.min(dist.tempo)])
        })
        ex.name <- unique(unlist(ls.name))

        stop("Invalid argument \'param\': some of the coefficient names does not match those of the estimate model \n",
             "incorrect names: \"",paste(incorrect.name, collapse = "\" \""),"\" \n",
             "example of valid names: \"",paste(ex.name, collapse = "\" \""),"\"\n")
    }

    if(checkType1 && any(coef.true[param]!=0)){
        txtCoef <- paste(param[coef.true[param]!=0], collapse = "\" \"")
        stop("Control type 1 error: coefficients \"",txtCoef,"\" are not 0 while their belong to the param hypothesis\n")
    }
    if(checkType2 && any(coef.true[param]==0)){
        txtCoef <- paste(param[coef.true[param]==0], collapse = "\" \"")
        stop("Control type 2 error: coefficients \"",txtCoef,"\" are 0 while their belong to the param hypothesis\n")
    }
    res.C <- createContrast(param, name.param = name.coef, add.rowname = TRUE, rowname.rhs = FALSE)
    contrast <- res.C$contrast
    if(is.null(null)){
        rhs <- res.C$null
    }else{
        if(length(null)!=n.param){
            stop("Argument \'null\' have the same length as argument \'param\' \n")
        }
        if(!is.null(names(null)) && any(names(null)!=param)){
            stop("If named, argument \'null\' must have for name argument \'param\' \n")
        }
        rhs <- null
    }

    ## *** warmup
    if(!is.null(warmup) && inherits(warmup,"lvm")){
        stop("Argument \'warmup\' must be a list of lvm objects \n")
    }
 
### ** display
    if(trace>1){
        cat("  Settings: \n")
        cat("  > simulation for n=",paste(n,collapse = " "),"\n",sep="")
        cat("  > model: \n")
        print(object)
        cat("  > expected coefficients: \n")
        print(coef.true)
        cat("  > null hypotheses: \n")
        print(rhs)
        if(bootstrap){
            cat("  > bootstrap: ",bootstrap,"\n")
        }
        if(!is.null(seed)){
            cat("  > seed: ",seed,"\n")
        }        
    }

    
### ** loop
    store.coef <- param
    if(F.test){
        store.coef <- c(store.coef, "global")
    }
    n.store <- length(store.coef)
    
    if(trace>1){
        cat("\n")
        cat(" Perform simulation: \n")
    }
    if(cpus>1){
        ## *** define cluster
        if(trace>0){
            cl <- suppressMessages(parallel::makeCluster(cpus, outfile = ""))
            pb <- utils::txtProgressBar(max = n.rep, style = 3)          
        }else{
            cl <- parallel::makeCluster(cpus)
        }
        ## *** link to foreach
        doParallel::registerDoParallel(cl)

        ## *** export package
        parallel::clusterCall(cl, fun = function(x){
            suppressPackageStartupMessages(requireNamespace("lava", quietly = TRUE))
            suppressPackageStartupMessages(requireNamespace("lavaSearch2", quietly = TRUE))
        })
        
        ## *** seed
        cpus.name <- unlist(parallel::clusterApply(cl = cl, 1:cpus, function(x){
            myName <- paste(Sys.info()[["nodename"]], Sys.getpid(), sep="-")
            return(myName)
        }))
        if(length(seed)==0){
            seed <- rep(as.numeric(NA), cpus)
        }else if(length(seed)==1){
            set.seed(seed)
            seed <- sample(1:1e5,size=cpus,replace=FALSE)                
        }else{
            if(length(seed)!=cpus){
                stop("Length of argument \'seed\' does not match argument \'cpus\' \n")
            }            
        }
        names(seed) <- cpus.name

        
        ## *** parallel computation
        toExport <- c(".warperType1", "cpus.name")

        iRep <- NULL # [:for CRAN check] foreach
        resSim <- foreach::`%dopar%`(
                               foreach::foreach(iRep = 1:n.rep,
                                                .export = toExport,
                                                .packages = c("lavaSearch2")),{ # iRep <- 1

                                                    if(trace>0){utils::setTxtProgressBar(pb, iRep)}

                                                    myName <- paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
                                                    iSeed <- as.double(seed[myName])
                                                    
                                                    ls.pvalue <- vector(mode = "list", length = n.n)
                                                    ls.estimate <- vector(mode = "list", length = n.n)
                                                    iIndex <- 1
                                                    
                                                    for(iN in n){ # iN <- n[1]
                                                        iRes <- .warperType1(iRep,
                                                                             n = iN,
                                                                             correction = correction,
                                                                             generative.object = generative.object,
                                                                             generative.coef = generative.coef,
                                                                             object = object,
                                                                             warmup = warmup,
                                                                             cluster = cluster,
                                                                             coef.true = coef.true,
                                                                             type.coef = type.coef,
                                                                             name.coef = name.coef,
                                                                             store.coef = store.coef,
                                                                             n.coef = n.coef,
                                                                             n.store = n.store,
                                                                             F.test = F.test,
                                                                             param = param,
                                                                             contrast = contrast,
                                                                             rhs = rhs,
                                                                             bootstrap = bootstrap,
                                                                             n.bootstrap = n.bootstrap,
                                                                             seed = iSeed,
                                                                             ...)

                                                        ls.pvalue[[iIndex]] <- iRes$pvalue
                                                        ls.estimate[[iIndex]] <- iRes$estimate
                                                        iIndex <- iIndex + 1
                                                    }
                                                    
                                                    return(list(pvalue = do.call("rbind",ls.pvalue),
                                                                estimate = do.call("rbind",ls.estimate)))
                                                })
        
        parallel::stopCluster(cl)
        if(trace>0){close(pb)}
        ## *** post process
        dt.pvalue <- do.call("rbind",lapply(resSim,"[[","pvalue"))
        dt.pvalue <- dt.pvalue[order(dt.pvalue$n,dt.pvalue$rep),,drop=FALSE]
        dt.estimate <- do.call("rbind",lapply(resSim,"[[","estimate"))
        dt.estimate <- dt.estimate[order(dt.estimate$n,dt.estimate$rep),,drop=FALSE]
        
    }else{

        ## *** filename
        if(is.null(label.file)){label.file <- seed}
        filename_tempo.pvalue <- paste0("type1error-S",label.file,"(tempo).rds")
        filename_tempo.estimate <- paste0("estimate-S",label.file,"(tempo).rds")
        filename.pvalue <- gsub("\\(tempo\\)","",filename_tempo.pvalue)
        filename.estimate <- gsub("\\(tempo\\)","",filename_tempo.estimate)

        if(!is.null(dir.save)){
            validPath(dir.save, type = "dir")
        }

        if(!is.null(dir.save)){
            cat("  > export results in ",dir.save,"\n")
        }

        if(trace>0){
            test.package <- try(requireNamespace("pbapply"), silent = TRUE)
            if(inherits(test.package,"try-error")){
                stop("There is no package \'pbapply\' \n",
                     "This package is necessary when argument \'trace\' is TRUE \n")
            }
            FCTapply <- pbapply::pblapply
        }else{
            FCTapply <- lapply
        }
        if(!is.null(seed)){
            set.seed(seed)
        }else{
            seed <- as.numeric(NA)
        }
        
        ## *** sequential simulation
        dt.pvalue <- NULL
        dt.estimate <- NULL

        for(iN in n){ ## iN <- n[1]

            if(trace>0){
                cat("  > sample size=",iN,"\n", sep = "")                    
            }

            resSim <- do.call(FCTapply, args = list(X = 1:n.rep, FUN = function(iRep){
                .warperType1(iRep,
                             n = iN,
                             correction = correction,
                             generative.object = generative.object,
                             generative.coef = generative.coef,
                             object = object, warmup = warmup, cluster = cluster,                             
                             coef.true = coef.true, type.coef = type.coef, name.coef = name.coef,
                             store.coef = store.coef, n.coef = n.coef, n.store = n.store,
                             F.test = F.test, param = param, contrast = contrast, rhs = rhs,
                             bootstrap = bootstrap,
                             n.bootstrap = n.bootstrap,
                             seed = seed,
                             ...)
            }))
            dt.pvalue <- rbind(dt.pvalue,
                               do.call("rbind",lapply(resSim,"[[","pvalue"))
                               )
            dt.estimate <- rbind(dt.estimate,
                                 do.call("rbind",lapply(resSim,"[[","estimate")))
                             
            
            ## export (tempo)
            if(!is.null(dir.save)){
                saveRDS(dt.pvalue, file = file.path(dir.save,filename_tempo.pvalue))
                saveRDS(dt.estimate, file = file.path(dir.save,filename_tempo.estimate))
            }
            
        }   
    }

    ## ** export
    if(!is.null(dir.save)){
        saveRDS(dt.pvalue, file = file.path(dir.save,filename.pvalue))
        saveRDS(dt.estimate, file = file.path(dir.save,filename.estimate))
    }
    out <- list(estimate = dt.estimate,
                p.value = dt.pvalue,
                e.true = e.true,
                param = param)
    class(out) <- append("calibrateType1",class(out))
    return(out)


}

## * calibrateType1.lvmfit
##' @rdname calibrateType1
##' @export
calibrateType1.lvmfit <- function(object, param, n.rep, correction = TRUE, F.test = FALSE,
                                  bootstrap = FALSE, n.bootstrap = 1e3,
                                  seed = NULL, trace = 2, cpus = 1, ...){

    ## ** Prepare
    ## *** model
    object.model <- object$model

    ## *** coef
    coef.true <- coef(object)
    name.coef <- names(coef.true)
    if(any(param %in% name.coef == FALSE)){
        txt <- param[param %in% name.coef == FALSE]
        txt2 <- setdiff(name.coef, param)
        stop("Argument \'param\' does not match the names of the model coefficients \n",
             "Incorrect param: \"",paste(txt, collapse = "\" \""),"\" \n",
             "Possible param: \"",paste(txt2, collapse = "\" \""),"\" \n")
    }
    coef.true[param] <- 0

    ## *** data
    n <- object$data$n

    ## ** Run
    out <- calibrateType1(object.model,
                          param = param,
                          n.rep = n.rep,
                          n = n,
                          correction = correction,
                          F.test = F.test,
                          generative.object = object.model,
                          generative.coef = coef.true, 
                          true.coef = coef.true,              
                          bootstrap = bootstrap,
                          n.bootstrap = n.bootstrap,
                          checkType1 = FALSE,
                          checkType2 = FALSE,
                          dir.save = NULL,
                          label.file = NULL,             
                          seed = seed,
                          cpus = cpus,
                          trace = trace,
                          ...)


    ## ** Export
    return(out)
}

## * .warperType1
.warperType1 <- function(iRep, n, correction, generative.object, generative.coef,
                         object, warmup, cluster,
                         coef.true, type.coef, name.coef, store.coef, n.coef, n.store,
                         F.test, param, contrast, rhs,
                         bootstrap, n.bootstrap,
                         seed,
                         ...){

    ls.iE <- list(estimate.truth = coef.true)  ## temporary
    ls.iP <- list()  ## temporary
    out <- list()
    dots <- list(...)
    if("control" %in%  names(dots) == FALSE){
        dots$control <- list()
    }
    ## ** simulation
    dt.sim <- lava::sim(generative.object, n = n, p = generative.coef, latent = FALSE)

    ## ** initialisation
    if(!is.null(warmup)){
        n.warmup <- length(warmup)
        for(iW in 1:n.warmup){ ## iW <- 1
            iE.lvm <- suppressWarnings(do.call(lava::estimate, args = c(list(warmup[[iW]], data = dt.sim, cluster = cluster), dots)))
            if(("convergence" %in% names(iE.lvm$opt)) && (iE.lvm$opt$convergence==1)){return(list(pvalue=NULL,estimate=NULL))} ## exclude lvm that has not converged
            dots$control$start <- coef(iE.lvm)
        }
    }
    
    ## ** model fit
    e.lvm <- suppressWarnings(do.call(lava::estimate, args = c(list(object, data = dt.sim, cluster = cluster), dots)))
    eS.lvm <- suppressWarnings(try(summary(e.lvm)$coef, silent = TRUE))
    
    ## check correct convergence of the latent variable model
    if(("convergence" %in% names(e.lvm$opt)) && (e.lvm$opt$convergence==1)){return(list(pvalue=NULL,estimate=NULL))} ## exclude lvm that has not converged
    if(any(eigen(getVarCov2(e.lvm))$values<=0)){return(list(pvalue=NULL,estimate=NULL))} ## exclude lvm where the residual covariance matrix is not semipositive definite
    ratio_sd_beta <- sqrt(diag(vcov(e.lvm)))/(abs(coef(e.lvm))+1)
    if(max(na.omit(ratio_sd_beta))>1e3){return(list(pvalue=NULL,estimate=NULL))} ## exclude if standard error much larger than coefficient
    if(inherits(eS.lvm, "try-error")){return(list(pvalue=NULL,estimate=NULL))} ## exclude lvm where we cannot compute the summary

    ## ** corrections
    if(correction){
        e.lvm.Satt <- e.lvm    
        testError.Satt <- try(sCorrect(e.lvm.Satt) <- FALSE, silent = TRUE)
        e.lvm.KR <- e.lvm
        testError.KR <- try(suppressWarnings(sCorrect(e.lvm.KR, safeMode = TRUE) <- TRUE), silent = TRUE)
    }else{
        e.lvm.Satt <- 1
        class(e.lvm.Satt) <- "try-error"
        e.lvm.KR <- 1
        class(e.lvm.KR) <- "try-error"
    }
    
    ## ** extract p.values
    eS.ML <- summary2(e.lvm, robust = FALSE, df = FALSE, bias.correct = FALSE)$coef
    F.ML <- compare2(e.lvm, robust = FALSE, df = FALSE, bias.correct = FALSE,
                     contrast = contrast, null = rhs, F.test = F.test, as.lava = FALSE)

    eS.robustML <- summary2(e.lvm, robust = TRUE, df = FALSE, bias.correct = FALSE)$coef
    F.robustML <- compare2(e.lvm, robust = TRUE, df = FALSE, bias.correct = FALSE,
                           contrast = contrast, null = rhs, F.test = F.test, as.lava = FALSE)
    

    if(!inherits(e.lvm.Satt,"try-error")){
        
        eS.Satt <- summary2(e.lvm.Satt, robust = FALSE)$coef
        F.Satt <- compare2(e.lvm.Satt, robust = FALSE,
                           contrast = contrast, null = rhs, F.test = F.test,
                           as.lava = FALSE)
            
        eS.robustSatt <- summary2(e.lvm.Satt, robust = TRUE)$coef
        F.robustSatt <- compare2(e.lvm.Satt, robust = TRUE,
                                  contrast = contrast, null = rhs, F.test = F.test,
                                  as.lava = FALSE)
    }
    
    if(!inherits(e.lvm.KR,"try-error")){

        ## eS.SSC <- summary2(e.lvm.KR, robust = FALSE, df = FALSE)$coef
        F.SSC <- compare2(e.lvm.KR, robust = FALSE, df = FALSE,
                           contrast = contrast, null = rhs, F.test = F.test,
                          as.lava = FALSE)
        
        ## eS.robustSSC <- summary2(e.lvm.KR, robust = TRUE, df = FALSE)$coef
        F.robustSSC <- compare2(e.lvm.KR, robust = TRUE, df = FALSE,
                                 contrast = contrast, null = rhs, F.test = F.test,
                                 as.lava = FALSE)

        eS.KR <- summary2(e.lvm.KR, robust = FALSE)$coef
        F.KR <- compare2(e.lvm.KR, robust = FALSE,
                          contrast = contrast, null = rhs, F.test = F.test,
                          as.lava = FALSE)
        
        eS.robustKR <- summary2(e.lvm.KR, robust = TRUE)$coef
        F.robustKR <- compare2(e.lvm.KR, robust = TRUE,
                                contrast = contrast, null = rhs, F.test = F.test,
                                as.lava = FALSE)
    }
    
    ## ** store

    ## *** estimates
    ls.iE$estimate.ML <- eS.ML[name.coef,"Estimate"]
    if(!inherits(e.lvm.KR,"try-error")){
        ls.iE$estimate.MLcorrected <- eS.KR[name.coef,"Estimate"]
    }else{
        ls.iE$estimate.MLcorrected <- rep(as.numeric(NA), n.coef)
    }

    ## *** standard errors
    if(is.null(cluster)){
        ls.iE$se.ML <- eS.ML[name.coef,"Std. Error"]
        if(!inherits(e.lvm.KR,"try-error")){
            ls.iE$se.MLcorrected <- eS.KR[name.coef,"Std. Error"]
        }else{
            ls.iE$se.MLcorrected <- rep(as.numeric(NA), n.coef)
        }

        ls.iE$se.robustML <- eS.robustML[name.coef,"robust SE"]
        if(!inherits(e.lvm.KR,"try-error")){
            ls.iE$se.robustMLcorrected <- eS.robustKR[name.coef,"robust SE"]
        }else{
            ls.iE$se.robustMLcorrected <- rep(as.numeric(NA), n.coef)
        }
    }

    ## *** degree of freedom
    if(!inherits(e.lvm.Satt,"try-error")){
        if(is.null(cluster)){
            ls.iE$df.ML <- eS.Satt[name.coef,"df"]
        }
        ls.iE$df.robustML <- eS.robustSatt[name.coef,"df"]
    }
    
    if(!inherits(e.lvm.KR,"try-error")){
        if(is.null(cluster)){
            ls.iE$df.MLcorrected <- eS.KR[name.coef,"df"]
        }
        ls.iE$df.robustMLcorrected <- eS.robustKR[name.coef,"df"]
    }
    
    ## *** p-value
    if(is.null(cluster)){
        ls.iP$p.Ztest <- F.ML[store.coef,"p-value"]
        if(!inherits(e.lvm.Satt,"try-error")){
            ls.iP$p.Satt <- F.Satt[store.coef,"p-value"]
        }
        if(!inherits(e.lvm.KR,"try-error")){
            ls.iP$p.SSC <- F.SSC[store.coef,"p-value"]
            ls.iP$p.KR <- F.KR[store.coef,"p-value"]
        }
    }
    ls.iP$p.robustZtest <- F.robustML[store.coef,"p-value"]
    if(!inherits(e.lvm.Satt,"try-error")){
        ls.iP$p.robustSatt <- F.robustSatt[store.coef,"p-value"]
    }
    if(!inherits(e.lvm.KR,"try-error")){
        ls.iP$p.robustSSC <- F.robustSSC[store.coef,"p-value"]
        ls.iP$p.robustKR <- F.robustKR[store.coef,"p-value"]
    }

    ## *** niter.correct / warning
    if(!inherits(e.lvm.KR,"try-error")){
        test.warning <- inherits(attr(e.lvm.KR$sCorrect,"warning"),"try-error")
        niter.correct <- as.double(e.lvm.KR$sCorrect$opt$iterations)
    }else{
        test.warning <- NA
        niter.correct <- as.double(NA)
    }
    
    ## ** bootstrap
    if(bootstrap>0){
        e.boot <- eval(parse(text = "butils::bootReg(e.lvm, type = \"coef\", n.boot = n.bootstrap"))

        index.coef.boot <- match(param, name.coef)
        boot.perc <- summary(e.boot, p.value = TRUE, type = "perc", print = FALSE, index = index.coef.boot)
        boot.stud <- summary(e.boot, p.value = TRUE, type = "stud", print = FALSE, index = index.coef.boot)
        boot.bca <- summary(e.boot, p.value = TRUE, type = "bca", print = FALSE, index = index.coef.boot)

        ls.iP$p.bootPerc <- boot.perc[store.coef,"p.value"]
        ls.iP$p.bootStud <- boot.stud[store.coef,"p.value"]
        ls.iP$p.bootBca <- boot.bca[store.coef,"p.value"]
    }else{
        n.bootstrap <- as.numeric(NA)
    }

    ## ** collect results
    ## estimates
    ## use rep to avoid warning
    ## return(list(n = n,
                ## rep = iRep,
                ## seed = seed,
                ## ninter = niter.correct,
                ## warning = test.warning,
                ## name = names(coef.true),
                ## type = type.coef[name.coef]))
    df1 <- data.frame(n = n,
                      rep = iRep, n.coef,
                      seed = seed, n.coef,
                      niter = niter.correct,
                      warning = test.warning,
                      name = names(coef.true),
                      type = unname(type.coef[name.coef]),
                      stringsAsFactors = FALSE)
    rownames(df1) <- NULL

    df2 <- do.call(cbind,ls.iE)
    rownames(df2) <- NULL

    out$estimate <- cbind(df1, df2)

    ## p.value
    df1 <- data.frame(n = n,
                      rep = iRep,
                      seed = seed,
                      nboot = n.bootstrap,
                      niter = niter.correct,
                      warning = test.warning,
                      link = unname(store.coef),
                      stringsAsFactors = FALSE)
    rownames(df1) <- NULL

    df2 <- do.call(cbind,ls.iP)
    rownames(df2) <- NULL

    out$pvalue <- cbind(df1, df2)

    ## ** export
    return(out)
}
    
######################################################################
### calibrateType1.R ends here

## * .mergeCalibrateType1
.mergeCalibrateType1 <- function(list){

    
    if(any(sapply(list,function(iL){inherits(iL,"calibrateType1")})==FALSE)){
        stop("Argument \'list\' must only contain \"calibrateType1\" objects \n")
    }

    newObject <- list[[1]]
    newObject$estimate <- do.call(rbind,lapply(list,"[[","estimate"))
    newObject$p.value <- do.call(rbind,lapply(list,"[[","p.value"))

    return(newObject)
}
