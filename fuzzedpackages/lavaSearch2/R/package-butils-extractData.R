## * Documentation
#' @title Extract Data From a Model
#' 
#' @description Extract data from a model using \code{nlme::getData}, \code{riskRegression::coxDesign} or \code{model.frame}.. 
#' If it fails it will try to extract it by its name according to \code{model$call$data}.
#' 
#' @param object the fitted model.
#' @param design.matrix [logical] should the data be extracted after transformation (e.g. conversion of categorical variables to dummy variables)?
#' Otherwise the original data will be returned.
#' @param as.data.frame [logical] should the output be converted into a \code{data.frame} object?
#' @param envir [environment] the environment from which to search the data.
#'
#' @return a dataset.
#' 
#' @examples
#' set.seed(10)
#' n <- 101
#'
#' #### linear regression ####
#' Y1 <- rnorm(n, mean = 0)
#' Y2 <- rnorm(n, mean = 0.3)
#' Id <- findInterval(runif(n), seq(0.1,1,0.1))
#' data.df <- rbind(data.frame(Y=Y1,G="1",Id = Id),
#'            data.frame(Y=Y2,G="2",Id = Id)
#'            )
#' m.lm <- lm(Y ~ G, data = data.df)
#' a <- extractData(m.lm, design.matrix = TRUE)
#' b <- extractData(m.lm, design.matrix = FALSE)
#' 
#' library(nlme)
#' m.gls <- gls(Y ~ G, weights = varIdent(form = ~ 1|Id), data = data.df)
#' c <- extractData(m.gls)
#' m.lme <- lme(Y ~ G, random = ~ 1|Id, data = data.df)
#' d <- extractData(m.lme)
#' 
#' library(lava)
#' e.lvm <- estimate(lvm(Y ~ G), data = data.df)
#' e <- extractData(e.lvm)
#' e <- extractData(e.lvm, design.matrix = TRUE)
#' 
#' #### survival #### 
#' library(survival)
#'
#' \dontrun{
#'   library(riskRegression) ## needs version >=1.4.3
#'   dt.surv <- sampleData(n, outcome = "survival")
#'   m.cox <- coxph(Surv(time, event) ~ X1 + X2, data = dt.surv, x = TRUE, y = TRUE)
#'   f <- extractData(m.cox, design.matrix = FALSE)
#'   f <- extractData(m.cox, design.matrix = TRUE)
#'   m.cox <- coxph(Surv(time, event) ~ strata(X1) + X2, data = dt.surv, x = TRUE, y = TRUE)
#'   f <- extractData(m.cox, design.matrix = TRUE)
#' }
#' 
#' #### nested fuuctions ####
#' fct1 <- function(m){
#'    fct2(m)
#' }
#' fct2 <- function(m){ 
#'    extractData(m)
#' }
#' g <- fct1(m.gls)
#' @concept extractor
#' @export
`extractData` <-
    function(object, design.matrix, as.data.frame, envir){
        UseMethod("extractData", object)
    }

## * method extractData.lm
#' @rdname extractData
#' @export
extractData.lm <- function(object, design.matrix = FALSE, as.data.frame = TRUE,
                           envir = environment()){
    ## ** check arguments
    validLogical(design.matrix, valid.length = 1)
    validLogical(as.data.frame, valid.length = 1)

    ## ** extract data
    if(design.matrix){
        data <- model.matrix(object)
    }else{
        ## cannot use model.frame because it only returned the part of the dataset relevant for fitting the model
        ## this is not enougth for modelsearch2
        ## data <- try(model.frame(object), silent = TRUE)
        ## data <- object$model
        data <- evalInParentEnv(object$call$data)
        if("function" %in% class(data)){
            stop("data has the same name as a function \n",
                 "consider renaming data before generating object \n")
        }
        if(!inherits(data, "data.frame")){
            stop("Could not extract the data from the model \n")
        } 
    }

    ## ** normalize data
    if(as.data.frame){
        data <- as.data.frame(data)        
    }

    ## ** export
    return(data)
}

## * method extractData.coxph
#' @rdname extractData
#' @export
extractData.coxph <- function(object, design.matrix = FALSE, as.data.frame = TRUE,
                              envir = environment()){
    ## ** check arguments
    validLogical(design.matrix, valid.length = 1)
    validLogical(as.data.frame, valid.length = 1)

    ## ** extract data
    if(design.matrix){
         tryPkg <- requireNamespace("riskRegression")
         if("try-error" %in% class(tryPkg)){
            stop(tryPkg)
        }else if(utils::packageVersion("riskRegression")<="1.4.3"){
            stop("riskRegression version must be > 1.4.3 \n",
                 "latest version available on Github at tagteam/riskRegression \n")
        }else{
            #### [:toUpdate]
            ##  data <- try(riskRegression::coxDesign(object), silent = TRUE)
            ##  strataVar <- riskRegression::coxVariableName(object)$stratavars.original

            ## this is a temporary modification waiting for the update of riskRegression on CRAN
            coxDesign.rr <- get("coxDesign", envir = asNamespace("riskRegression"), inherits = FALSE)
            coxVariableName.rr <- get("coxVariableName", envir = asNamespace("riskRegression"), inherits = FALSE)
            data <- try(coxDesign.rr(object), silent = TRUE)
            strataVar <- coxVariableName.rr(object)$stratavars.original
        } 
      
        if(length(strataVar)>0){         
            data2 <- evalInParentEnv(object$call$data)
            data <- cbind(as.data.frame(data),
                          as.data.frame(data2)[,strataVar,drop=FALSE])        
        }
    }else{
        data <- evalInParentEnv(object$call$data)
        
        if("function" %in% class(data)){
            stop("data has the same name as a function \n",
                 "consider renaming data before generating object \n")
        }
        if(!inherits(data, "data.frame")){
            stop("Could not extract the data from the model \n")
        } 
    }

    ## ** normalize data
    if(as.data.frame){
        data <- as.data.frame(data)        
    }

    ## ** export
    return(data)
    
}

## * method extractData.cph
#' @rdname extractData
#' @export
extractData.cph <- extractData.coxph

## * method extractData.lvmfit
#' @rdname extractData
#' @export
extractData.lvmfit <- function(object, design.matrix = FALSE, as.data.frame = TRUE,
                               envir = environment()){
    ## ** check arguments
    validLogical(design.matrix, valid.length = 1)
    validLogical(as.data.frame, valid.length = 1)

    ## ** extract data
    if(design.matrix){
        data <- object$data$model.frame
        keep.cols <- intersect(c("(Intercept)",lava::vars(object)), names(data))
        data <- data[,keep.cols,drop=FALSE]
    }else{
        data <- evalInParentEnv(object$call$data)
        
        if("function" %in% class(data)){
            stop("data has the same name as a function \n",
                 "consider renaming data before generating object \n")
        }
        
        if(!inherits(data, "data.frame")){
            data <- model.frame(object)
        }
    }

    ## ** normalize data
    if(as.data.frame){
        data <- as.data.frame(data)        
    }

    ## ** export
    return(data)
    
}

## * method extractData.gls
#' @rdname extractData
#' @export
extractData.gls <- function(object, design.matrix = FALSE, as.data.frame = TRUE,
                            envir = environment()){
    ## ** check arguments
    validLogical(design.matrix, valid.length = 1)
    validLogical(as.data.frame, valid.length = 1)

    ## ** extract data
    if(design.matrix){

        # assign the dataset to the object if not in the current environment
        name.data <- as.character(object$call$data)
        if((length(name.data) == 1) && (name.data %in% ls() == FALSE)){
            object$data <- evalInParentEnv(object$call$data)
        }
      
        data <- try(nlme::getData(object), silent = TRUE)

    }else{
        data <- evalInParentEnv(object$call$data)
        
        if("function" %in% class(data)){
            stop("data has the same name as a function \n",
                 "consider renaming data before generating object \n")
        }
        if(!inherits(data, "data.frame")){
            stop("Could not extract the data from the model \n")
        } 
    }

    ## ** normalize data
    if(as.data.frame){
        data <- as.data.frame(data)        
    }

    ## ** export
    return(data)
    
}

## * method extractData.lme
#' @rdname extractData
#' @export
extractData.lme <- extractData.gls
