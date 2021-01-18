### methods-modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 22 2017 (16:43) 
## Version: 
## last-updated: feb 18 2019 (09:33) 
##           By: Brice Ozenne
##     Update #: 216
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * nStep
## ** documentation - nStep
#' @title Find the Number of Steps Performed During the Sequential Testing
#' @description Find the number of steps performed during the sequential testing.
#' @name nStep
#' 
#' @param object a \code{modelsearch2} object.
#' 
#' @return an integer.
#' 
#' @examples
#' \dontrun{
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, method.p.adjust = "holm")
#' nStep(res)
#' }
#' 
#' @concept modelsearch
#' @concept extractor
#' @export
`nStep` <-
  function(object) UseMethod("nStep")

## ** function - nStep
#' @rdname nStep
#' @export
nStep.modelsearch2 <- function(object){

    return(length(object$sequenceTest))
    
}

## * getStep
## ** documentation - getStep
#' @title Extract one Step From the Sequential Procedure
#' @description Extract one step from the sequential procedure.
#' @name getStep
#' 
#' @param object a \code{modelsearch2} object
#' @param step [integer >0] which test should be extracted?
#' @param slot [character] the element from the modelsearch2 object that should be extracted.
#'
#' @examples
#' \dontrun{
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, method.p.adjust = "holm")
#'
#' getStep(res)
#' getStep(res, slot = "sequenceTest")
#' getStep(res, step = 1)
#' }
#' @concept modelsearch
#' @concept extractor
#' @export
`getStep` <-
  function(object, step, slot) UseMethod("getStep")

## ** function - getStep
#' @rdname getStep
#' @export
getStep.modelsearch2 <- function(object, step = nStep(object), slot = NULL){
    
    ## ** normalize arguments
    lastStep <- nStep(object)
    if(any(step %in% 1:lastStep == FALSE)){
        stop("step must be an integer between 1 and ",lastStep,"\n")
    }
    if(!is.null(slot) && slot %in% names(object) == FALSE){
        stop("argument \'slot\' must be one of \"",paste(names(object),collapse="\" \""),"\"\n")
    }

    ## ** subset
    new.object <- object
    new.object$sequenceTest <- object$sequenceTest[step]
    new.object$sequenceModel <- object$sequenceModel[step]
    class(new.object) <- NULL

    ## ** export
    if(is.null(slot)){
        return(new.object)
    }else{
        if(is.list(new.object[[slot]])){
            return(new.object[[slot]][[1]])
        }else{
            return(new.object[[slot]])
        }
    }
}

## * getNewLink
## ** documentation - getNewLink
#' @title Extract the Links that Have Been Found by the modelsearch2.
#' @description Extract the links that have been found relevant by modelsearch2.
#' @name getNewLink
#' 
#' @param object a \code{modelsearch2} object.
#' @param step [logical] which test should be extracted?
#'
#' @return A character vector.
#' 
#' @examples
#' \dontrun{
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#'
#' set.seed(10)
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, method.p.adjust = "holm")
#' getNewLink(res)
#' }
#' @concept modelsearch
#' @concept extractor
#' 
#' @export
`getNewLink` <-
  function(object, step) UseMethod("getNewLink")

## ** function - getNewLink
#' @rdname getNewLink
#' @export
getNewLink.modelsearch2 <- function(object, step = 1:nStep(object)){

    selected <- link <- NULL
    
    ## ** normalize arguments
    lastStep <- nStep(object)
    if(any(step %in% 1:lastStep == FALSE)){
        stop("step must be an integer between 1 and ",lastStep,"\n")
    }

    ## ** extract
    ls.link <- lapply(step, function(x){
        iStep <- getStep(object,step=x,slot="sequenceTest")
        return(subset(iStep, subset = selected == TRUE, select = "link", drop = TRUE))
    })

    return(unlist(ls.link))    
}

## * getNewModel
## ** documentation - getNewModel
#' @title Extract the Model that Has Been Retains by the modelsearch2.
#' @description Extract the model that has been retained by modelsearch2.
#' @name getNewModel
#' 
#' @param object a \code{modelsearch2} object.
#' @param step [integer >=0] the step at which the model should be extracted.
#' 0 returns the initial model, i.e. before adding any links.
#'
#' @return A \code{lvmfit} object.
#' 
#' @examples
#' \dontrun{
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#'
#' set.seed(10)
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, method.p.adjust = "holm")
#' getNewModel(res)
#' }
#' @concept modelsearch
#' @concept extractor
#' 
#' @export
`getNewModel` <-
  function(object, step) UseMethod("getNewModel")

## ** function - getNewLink
#' @rdname getNewModel
#' @export
getNewModel.modelsearch2 <- function(object, step = nStep(object)){

    ## ** normalize arguments
    lastStep <- nStep(object)
    if(any(step %in% 0:lastStep == FALSE)){
        stop("Argument \'step\' must be an integer between 0 and ",lastStep,"\n")
    }
    if(length(step)!=1){
        stop("Argument \'step\' must have length one \n")
    }
    
    ## ** extract
    if(step>0){
        return(object$sequenceModel[[step]])
    }else{
        return(object$initialModel)
        }
}

## * autoplot
## ** documentation - autoplot
#' @title Display the Value of a Coefficient across the Steps.
#' @description Display the value of a coefficient across the steps.
#' @name autplot-modelsearch2
#' 
#' @param object a \code{modelsearch2} object.
#' @param param [character vector] the name of the coefficient(s) to be displayed.
#' @param ci [logical] should the confidence intervals of the coefficient(s) be displayed.
#' @param conf.level [numeric, 0-1] confidence level of the interval.
#' @param step [integer >0] the steps at which the coefficient value should be displayed.
#' @param plot [logical] should the graph be displayed?
#' @param add.0 [logical] should an horizontal line representing no effect be displayed?
#' @param ...  [internal] only used by the generic method.
#'
#' @return A list containing \itemize{
#' \item \code{plot}: a ggplot object.
#' \item \code{data}: the data used to generate the ggplot object.
#' }
#' 
#' @examples
#' \dontrun{
#' mSim <- lvm(Y~G+X1+X2+X3+X4+X5)
#' addvar(mSim) <- ~Z1+Z2
#'
#' set.seed(10)
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+X3+X4+X5+Z1+Z2
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, method.p.adjust = "holm", alpha = 0.05)
#' autoplot(res, param = "Y~G")
#' autoplot(res, param = c("Y","Y~G"))
#' }
#' @concept modelsearch
#' @concept extractor
#' 

## ** code - autoplot
#' @rdname autoplot.modelsearch2
#' @method autoplot modelsearch2
#' @export
autoplot.modelsearch2 <- function(object, param, ci = TRUE, step = 0:nStep(object), conf.level = 0.95,
                                  plot = TRUE, add.0 = TRUE, ...){

    ## ** check arguments
    lastStep <- nStep(object)
    if(any(step %in% 0:lastStep == FALSE)){
        stop("Argument \'step\' must be an integer between 0 and ",lastStep,"\n")
    }
    if(any(diff(step)<=0)){
        stop("Argument \'step\' must be a vector of strictly increasing integers \n")
    }
    
    name.coef <- names(coef(object$initialModel))
    if(any(param %in% name.coef == FALSE)){
        stop("Argument \'param\' incorrect: some parameters do not belong to the initial model \n",
             "invalid parameters: \"",paste0(param[param %in% name.coef == FALSE], collapse ="\" \""),"\"\n")
    }

    ## ** get data
    n.points <- length(step)
    df.plot <- NULL

    for(iStep in 1:length(step)){
        iModel <- getNewModel(object, step = step[iStep])
        iCoef <- unname(cbind(coef(iModel)[param],stats::confint(iModel, level = conf.level)[param,,drop=FALSE]))
        colnames(iCoef) <- c("estimate", "ci.inf", "ci.sup")
        df.plot  <- rbind(df.plot,
                          data.frame(step = step[iStep], param  = param, iCoef)
                          )
    }

    ## ** graphical display
    gg <- ggplot2::ggplot(df.plot, aes_string(x = "step", y = "estimate"))
    if(add.0){
        gg <- gg + ggplot2::geom_hline(yintercept = 0, color = "red")
    }
    gg <- gg + ggplot2::geom_point() + ggplot2::geom_line()
    gg <- gg + ggplot2::facet_wrap(~param)
    if(ci){
        gg <- gg + ggplot2::geom_errorbar(aes_string(ymin = "ci.inf", ymax = "ci.sup"))        
        gg <- gg + ggplot2::ylab(paste0("estimate [",round(conf.level*100,2),"% confidence interval]"))
    }
    if(plot){
        print(gg)
    }

    ## ** export
    return(invisible(list(plot = gg,
                          data = df.plot)))

    

}
