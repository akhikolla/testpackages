#' @title Test deletion of edge from an interaction model
#' 
#' @description Tests if an edge can be deleted from an interaction
#'     model.
#' 
#' @details If the model is decomposable and the edge is contained in
#'     one clique only then the test is made in the marginal model
#'     given by that clique. In that case, if the model is a
#'     log-linear model then degrees of freedom are adjusted for
#'     sparsity
#' 
#' @aliases testdelete testdelete.iModel print.testdelete testdelete.mModel
#' 
#' @param object A model; an object of class \code{iModel}.
#' @param edge An edge in the model; either as a right-hand sided
#'     formula or as a vector
#' @param k Penalty parameter used when calculating change in AIC
#' @param details The amount of details to be printed; 0 surpresses
#'     all information
#' @param \dots Further arguments to be passed on to the underlying
#'     functions for testing.
#' @return A list.
#' @details If model is decomposable and edge is in one clique only,
#'     then degrees of freedom are adjusted for sparsity
#' 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{testadd}}
#' @keywords models htest
#' @examples
#' 
#' ## Discrete models
#' data(reinis)
#'
#' ## A decomposable model
#' mf <- ~smoke:phys:mental + smoke:systol:mental
#' object <- dmod(mf, data=reinis)
#' testdelete(object, c("phys", "mental"))
#' testdelete(object, c("smoke", "mental"))
#' 
#' ## A non-decomposable model
#' mf <- ~smoke:phys + phys:mental + smoke:systol + systol:mental
#' object <- dmod(mf, data=reinis)
#' 
#' testdelete(object, c("phys", "mental"))
#' 
#' ## Continuous models
#' data(math)
#'
#' ## A decomposable model
#' mf <- ~me:ve:al + me:al:an
#' object <- cmod(mf, data=math)
#' testdelete(object, c("ve", "al"))
#' testdelete(object, c("me", "al"))
#' 
#' ## A non-decomposable model
#' mf <- ~me:ve + ve:al + al:an + an:me
#' object <- cmod(mf, data=math)
#' testdelete(object, c("me", "ve"))


#' @export testdelete 
testdelete <- function(object, edge, k=2, details=1, ...)
  UseMethod("testdelete")

#' @export 
testdelete.iModel <- function(object, edge, k=2, details=1, ...){

    model.type <- class(object)[1]
    edge <- rhsFormula2list(edge)[[1]]
##    cat(sprintf("testdelete.iModel model.type=%s edge=%s\n", model.type, toString(edge)))

    .is_valid_edge(edge, object)

    if (is.null((amat <- list(...)$amat)))
        amat <- .as_amat(getmi(object, "glist"))
        
    if (amat[edge[1], edge[2]] != 1)
        stop(cat("edge:", edge, "not in model\n"))

    ##model.type="mModel"
    if (model.type == "mModel")
        ans <- .test_delete_edge(object, edge, details)
    else {
        ## Assume glist is minimal, i.e. no redundant elements
        ## Is model graphical? Is model decomposable?
        cliq     <- maxCliqueMAT(amat)$maxCliques
        isgraph  <- length(cliq) == length(getmi(object, "glist")) 
        isdecomp <- length(mcsMAT(amat)) > 0
        
        ## Is edge only in one clique in decomposable model?
        onlyinone <- FALSE
        if (isdecomp){
            ##idx   <- isin (cliq, edge, index=TRUE)
            idx   <- is_inset (edge, cliq, index=TRUE)
            onlyinone <- sum(idx) == 1
        }
        
        if (isdecomp && onlyinone && model.type %in% c("cModel", "dModel")){
            ## cat("If edge is in one clique only, do test in marginal table\n")
            hostcq <- cliq[idx == 1][[1]]
            ans <- .test_in_one_clique(object, edge, hostcq, details)
        } else {
            ## Make usual LR-test
            ans <- .test_delete_edge(object, edge, details)
        }
    }

    ret <- .finalize_test(ans, k)
    class(ret) <- "testdelete"
    ##cat("testdelete - exit\n"); str(ret)
    ret
}

#' @title Test addition of edge to graphical model
#' 
#' @description Performs a test of addition of an edge to a graphical
#'     model (an \code{iModel} object).
#' 
#' @aliases testadd testadd.iModel print.testadd testadd.mModel
#' 
#' @details Let M0 be the model and e={u,v} be an edge and let M1 be
#'     the model obtained by adding e to M0. If M1 is decomposable AND
#'     e is contained in one clique C only of M1 then the test is
#'     carried out in the C-marginal model. In this case, and if the
#'     model is a log-linear model then the degrees of freedom is
#'     adjusted for sparsity.
#' 
#' @param object A model; an object of class \code{iModel}.
#' @param edge An edge; either as a vector or as a right hand sided
#'     formula.
#' @param k Penalty parameter used when calculating change in AIC
#' @param details The amount of details to be printed; 0 surpresses
#'     all information
#' @param \dots Further arguments to be passed on to the underlying
#'     functions for testing.
#' @return A list
#'
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{testdelete}}
#' @keywords models htest
#' @examples
#' 
#' ## Discrete models
#' data(reinis)
#' 
#' ## A decomposable model
#' mf <- ~smoke:phys:mental + smoke:systol:mental
#' object <- dmod(mf, data=reinis)
#' testadd(object, c("systol", "phys"))
#' 
#' ## A non-decomposable model
#' mf <- ~smoke:phys + phys:mental + smoke:systol + systol:mental
#' object <- dmod(mf, data=reinis)
#' testadd(object, c("phys", "systol"))
#' 
#' ## Continuous models
#' data(math)
#' 
#' ## A decomposable model
#' mf <- ~me:ve:al + al:an
#' object <- cmod(mf, data=math)
#' testadd(object, c("me", "an"))
#' 
#' ## A non-decomposable model
#' mf <- ~me:ve + ve:al + al:an + an:me
#' object <- cmod(mf, data=math)
#' testadd(object, c("me", "al"))


#' @export testadd
testadd <- function(object, edge, k=2, details=1, ...)
  UseMethod("testadd")

#' @export 
testadd.iModel <- function(object, edge, k=2, details=1, ...){

    model.type <- class(object)[1]
    edge <- rhsFormula2list(edge)[[1]]
    ##cat(sprintf("testadd.iModel model.type=%s\n", model.type))

    .is_valid_edge(edge, object)
        
    if (is.null((amat <- list(...)$amat)))
        amat <- .glist2adjMAT(.glist(object))  ## FIXME use better function
        
    if (amat[edge[1], edge[2]] != 0)
        stop(cat("edge:", edge, "already in model\n"))

    ##model.type="mModel"
    if (identical(model.type, "mModel"))
        ans <- .test_add_edge(object, edge, details)
    else {   
        ## Add edge to model FIXME: Fails if amat is sparse!
        amat[edge[1], edge[2]] <- amat[edge[2], edge[1]] <- 1L
        
        ## Is model graphical? Is model decomposable?
        cliq <- maxCliqueMAT(amat)$maxCliques
        isgraph <- length(cliq) == length(.glist(object))
        isdecomp <- length(mcsMAT(amat)) > 0
        
        ## Is edge only in one clique in decomposable model?
        onlyinone <- FALSE
        if (isdecomp){
            ##idx   <- isin (cliq, edge, index=TRUE)
            idx   <- is_inset (edge, cliq, index=TRUE)            
            onlyinone <- sum(idx) == 1
        }
        
        if (isdecomp && onlyinone && model.type %in% c("cModel","dModel")){
            ## If edge is in one clique only, do test in marginal table
            hostcq <- cliq[idx==1][[1]]
            ans <- .test_in_one_clique(object, edge, hostcq, details)
        } else {
            ## Make usual LR-test
            ans <- .test_add_edge(object, edge, details)
        }
    }

    ret <- .finalize_test(ans, k)    
    class(ret) <- "testadd"
    ret
}



#' @method print testdelete
#' @export 
print.testdelete <- function(x,  ...){

  cat(sprintf("dev: %8.3f df:%3i p.value: %7.5f AIC(k=%3.1f): %6.1f edge: %s \n",
              x$statistic, x$df, x$p.value, x$k, x$aic, .toString(x$edge,':')))
  if (x$conmethod=="data.based"){
    if (x$details > 0){
      cat("host: ", x$hostcq, "\n")
      cat("Notice: Test performed in saturated marginal model\n")
    }
  } else {
    if (x$details > 0){
      cat("Notice: Test perfomed by comparing likelihood ratios\n")
    }
  }
  invisible(x)
}

#' @method print testadd
#' @export 
print.testadd <- function(x,  ...){

  cat(sprintf("dev: %8.3f df:%3i p.value: %7.5f AIC(k=%3.1f): %6.1f edge: %s \n",
              x$statistic, x$df, x$p.value, x$k, x$aic, .toString(x$edge,':')))
  if (x$conmethod=="data.based"){
    if (x$details > 0){
      cat("host: ", x$hostcq, "\n")
      cat("Notice: Test performed in saturated marginal model\n")
    }
  } else {
    if (x$details > 0){
      cat("Notice: Test perfomed by comparing likelihood ratios\n")
    }
  }
  invisible(x)
}


## dot-functions below here

.is_valid_edge <- function(edge, object){
    if (length(edge) != 2)
        stop(paste("Not a valid edge: ", paste(edge, collapse=":"), " \n"))
    if (!subsetof(edge, getmi(object, "varNames")))
        stop(cat("variables:", edge, "not in model\n"))
}


.finalize_test <- function(ans, k=2){
    extra2 <- list(aic=ans$statistic - k * ans$df, k=k)
    c(ans, extra2)
}

.test_delete_edge <- function(object, edge, details=1){
    sml   <- update(object, list(drop.edge=edge))
    ans   <- .comparemodels(object, sml)
    extra <- list(edge=edge, hostcq=NULL, details=details, conmethod='model.based')
    c(ans, extra)
}

.test_add_edge <- function(object, edge, details=1){
    lrg   <- update(object, list(add.edge=edge))
    ans   <- .comparemodels(lrg, object)
    extra <- list(edge=edge, hostcq=NULL, details=details, conmethod='model.based')
    c(ans, extra)
}

## FIXME: 31/12/17: For loglin er dev allerede "ganget med 2"... Check for de andre modeller...
.comparemodels <- function(lrg, sml) {
    ##devdiff <- 2 * (getmi(sml, "dev") - getmi(lrg, "dev"))
    devdiff <- (getmi(sml, "dev") - getmi(lrg, "dev"))
    dfdiff  <- getmi(sml, "dimension")['df'] - getmi(lrg, "dimension")['df']
    ##str(list(devdiff=devdiff, dfdiff=dfdiff))
    list('statistic'=devdiff, 'df'=dfdiff, 'p.value'=1 - pchisq(devdiff, dfdiff))
}

.test_in_one_clique <- function(object, edge, hostcq, details=1, ...){
    model.type <- class(object)[1]
    set <- c(edge, setdiffPrim(hostcq, edge))
    
    ans <- switch(model.type,
                  "cModel"={ 
                      ciTest_mvn(list(cov=getmi(object, "S"),
                                      n.obs=getmi(object, "n")), set=set, ...)
                  },
                  "dModel"={
                      ciTest_table(getmi(object, "data"),
                                   set=set, slice.info=FALSE, ...)
                  })
            
    extra <- list(edge=edge, hostcq=hostcq, details=details, conmethod='data.based')
    c(ans, extra)
}


























































############################################################################







## #############################################################################


## testadd.iModel <- function(object, edge, k=2, details=1, ...){

##     model.type <- class(object)[1]
##     ##cat(sprintf("testadd.iModel model.type=%s\n", model.type))
    
##     edge <- rhsFormula2list(edge)[[1]]
##     if (length(edge)!=2)
##         stop(paste("Not a valid edge: ", paste(edge, collapse=":"), " \n"))
    
##     ## ----- START USING amat
##     if (is.null((amat <- list(...)$amat)))
##         amat <- .glist2adjMAT(object$glist)

##     ## Is edge is in model? stop if not
##     if (!subsetof(edge, colnames(amat)))
##         stop(cat("variables:", edge, "not in model\n"))

##     if (amat[edge[1],edge[2]]!=0)
##         stop(cat("edge:", edge, "already in model\n"))

##     ## Add edge to model FIXME: Fails if amat is sparse!
##     amat[edge[1], edge[2]] <- amat[edge[2], edge[1]] <- 1L

##     ## Is model graphical?
##     cliq <- maxCliqueMAT(amat)$maxCliques
##     isgraph <- length(cliq) == length(object$glist)
    
##     ## Is model decomposable?
##     isdecomp <- length(mcsMAT(amat)) > 0
##     ## ----- STOP USING amat
    
##     ## Is edge only in one clique in decomposable model?
##     onlyinone <- FALSE
##     if (isdecomp){
##         idx   <- isin (cliq, edge, index=TRUE)
##         onlyinone <- sum(idx) == 1
##     }

##   if (isdecomp && onlyinone && model.type %in% c("cModel","dModel")){
##       ## If edge is in one clique only, do test in marginal table
##       ##
##       hostcq <- cliq[idx==1][[1]]
      
##       set <- c(edge, setdiffPrim(hostcq, edge))
      
##       ans <- switch(model.type,
##                     "cModel"={ 
##                         ciTest_mvn(list(cov=getmi(object, "S"),
##                                         n.obs=getmi(object, "n")),
##                                    set=set, ...)
##                     },
##                     "dModel"={
##                         ciTest_table(getmi(object, "data"),
##                                      set=set, ...)
##                     })
            
##       extra <- list(edge=edge, hostcq=hostcq, details=details, conmethod='data.based')
##   } else {
##       ## Make usual LR-test
##       ##
##       ob2   <- update(object, list(add.edge=edge))
##       ans   <- .comparemodels(ob2,object)
##       extra <- list(edge=edge, hostcq=NULL, details=details, conmethod='model.based')
##   }
##     extra2 <- list(aic=-(ans$statistic - k * ans$df), k=k)
##     ret <- c(ans, extra, extra2)
##     class(ret) <- "testadd"
##     return(ret)
## }







##     ## cat("ans--------------\n"); print(ans)         
##     ## stop()
##     ## cat("ans2-------------\n");
##     ## ans2 <- .test_delete_edge(object, edge, details)
##     ## print(ans2)

##     ## ret <- .finalize_test(ans, k)
##     ## str(ret)
    
##     ## ans <- .test_delete_edge(object, edge, details)
##     ## oo <<- object
##     ## ee <<- edge
##     ## aa <<- ans



## .comparemodels <- function(m1, m2) {
##   devdiff <- 2 * (getmi(m2, "dev") - getmi(m1, "dev"))
##   dfdiff  <- getmi(m2, "dimension")['df'] - getmi(m1, "dimension")['df']
##   str(list(devdiff=devdiff, dfdiff=dfdiff))
##   list('statistic'=devdiff, 'df'=dfdiff, 'p.value'=1 - pchisq(devdiff, dfdiff))
## }




## testdelete.iModel <- function(object, edge, k=2, details=1, ...){

##     model.type <- class(object)[1]
##     ##cat(sprintf("testdelete.iModel model.type=%s\n", model.type))    

##     edge <- rhsFormula2list(edge)[[1]]
##     if (length(edge) !=2 )
##         stop(paste("Not a valid edge: ", paste(edge, collapse=":"), " \n"))

##     if (is.null((amat <- list(...)$amat)))
##         amat <- .as_amat(getmi(object, "glist"))
        
##     ## Is edge is in model? stop if not
##     if (!subsetof(edge, colnames(amat)))
##         stop(cat("variables:", edge, "not in model\n"))

##     if (amat[edge[1], edge[2]] != 1)
##         stop(cat("edge:", edge, "not in model\n"))
    
##     ## Is model graphical?     ## FIXME: fails if model contains redundant elements..
##     cliq    <- maxCliqueMAT(amat)$maxCliques
##     isgraph <- length(cliq) == length(getmi(object, "glist"))

##     ## Is model decomposable?
##     isdecomp <- length(mcsMAT(amat)) > 0
    
##     ## Is edge only in one clique in decomposable model?
##     onlyinone <- FALSE
##     if (isdecomp){
##         idx   <- isin (cliq, edge, index=TRUE)
##         onlyinone <- sum(idx) == 1
##     }
    
##     if (isdecomp && onlyinone && model.type %in% c("cModel", "dModel")){
##         ## If edge is in one clique only, do test in marginal table
##         hostcq <- cliq[idx == 1][[1]]
##         set    <- c(edge, setdiff(hostcq, edge))
##         ##cat(sprintf("CHK: edge: {%15s} hostcq: {%s}\n", toString(edge), toString(hostcq)))        
##         ans <- switch(model.type,
##                       "cModel"={ 
##                           ciTest_mvn(list(cov=getmi(object, "S"),
##                                           n.obs=getmi(object, "n")),
##                                      set=set, ...)
##                       },
##                       "dModel"={
##                           ciTest_table(getmi(object, "data"),
##                                        set=set, ...)
##                       })
        
##         extra <- list(edge=edge, hostcq=hostcq, details=details, conmethod='data.based')
##     } else {
##         cat(sprintf("CHK: Make usual LR-test - edge = {%s}\n", toString(edge)))
##         ob2   <- update(object, list(drop.edge=edge))
##         ans   <- .comparemodels(object, ob2)
##         extra <- list(edge=edge, hostcq=NULL, details=details, conmethod='model.based')
##     }
##     extra2 <- list(aic=ans$statistic - k * ans$df, k=k)
##     ret    <- c(ans, extra, extra2)
##     class(ret) <- "testdelete"
##     ret
## }
