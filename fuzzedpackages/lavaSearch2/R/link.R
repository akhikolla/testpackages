## * findNewLink
## ** doc findNewLink
#' @title Find all New Links Between Variables
#' @description Find all new links between variables (adapted from lava::modelsearch).
#' 
#' @name findNewLink
#' 
#' @param object a \code{lvm} object.
#' @param data [optional] a dataset used to identify the categorical variables when not specified in the \code{lvm} object.
#' @param exclude.var [character vector] all links related to these variables will be ignore.
#' @param type [character vector] the type of links to be considered: \code{"regression"}, \code{"covariance"}, or \code{"both"}, .
#' @param rm.latent_latent [logical] should the links relating two latent variables be ignored?
#' @param rm.endo_endo [logical] should the links relating two endogenous variables be ignored?
#' @param rm.latent_endo [logical] should the links relating one endogenous variable and one latent variable be ignored?
#' @param output [character] Specify \code{"names"} to return the names of the variables to link
#' or specify \code{"index"} to return their position.
#' @param ... [internal] only used by the generic method.
#'
#' @return A list containing:
#' \itemize{
#' \item M.links: a matrix with two columns indicating (by name or position) the exogenous and endogenous variable corresponding to each link.
#' \item links: the name of the additional possible links
#' \item directional: a logical vector indicating for each link whether the link is unidirectional (\code{TRUE}, i.e. regression link)
#' or bidirectional (\code{FALSE}, i.e. covariance link).
#' }
#' 
#' @examples
#' library(lava)
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' categorical(m,labels=c("M","F","MF")) <- ~X1
#' findNewLink(m, rm.endo = FALSE)
#' findNewLink(m, rm.endo = TRUE)
#' findNewLink(m, exclude.var = "X1")
#' 
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' 
#' findNewLink(m, rm.endo = FALSE)
#' findNewLink(m, rm.endo = TRUE)
#' findNewLink(m, rm.endo = TRUE, output = "index")
#' findNewLink(m, type = "covariance")
#' findNewLink(m, type = "regression")
#' 
#' @concept modelsearch
#' @export
`findNewLink` <-
  function(object, ...) UseMethod("findNewLink")

## ** method findNewLink.lvm
#' @export
#' @rdname findNewLink
findNewLink.lvm <- function(object,                            
                            data = NULL,
                            type = "both",
                            exclude.var = NULL,
                            rm.latent_latent= FALSE,
                            rm.endo_endo= FALSE,
                            rm.latent_endo= FALSE,
                            output = "names",
                            ...){

    match.arg(output, choices = c("names","index"))
    match.arg(type, choices = c("both","covariance","regression"))
    if(is.null(data)){        
        data <- lava::sim(object, n = 1)
    }
    
    ## *** convertion to dummy variable name for categorical variables
    xF <- lava_categorical2dummy(object, data)
    AP <- with(lava::index(xF$x), A + t(A) + P)
    latent.xF <-  latent(xF$x)
    endogenous.xF <- endogenous(xF$x)
    exogenous.xF <- exogenous(xF$x)

    if(!is.null(exclude.var)){
        exclude.var <- var2dummy(xF, exclude.var)
    }
    
    if( any(exclude.var %in% colnames(AP) == FALSE) ){
        wrong.var <- exclude.var[exclude.var %in% colnames(AP) == FALSE]
        stop("unknown variable to exclude \n",
             "variable(s): \"",paste(wrong.var, collapse = "\" \""),"\"\n")
    }

    ## *** loop over links
    restricted <- c()
    directional <- c()
    for (i in seq_len(ncol(AP) - 1)){
        for (j in seq(i + 1, nrow(AP))){

            var.i <- rownames(AP)[i]
            var.j <- rownames(AP)[j]
            
            if(!is.null(exclude.var) && (var.i %in% exclude.var || var.j %in% exclude.var)){
                next
            }

            isLatent.i <- var.i %in% latent.xF
            isLatent.j <- var.j %in% latent.xF
            isEndogenous.i <- var.i %in% endogenous.xF
            isEndogenous.j <- var.j %in% endogenous.xF
            isExogenous.i <- var.i %in% exogenous.xF
            isExogenous.j <- var.j %in% exogenous.xF

            if(rm.latent_latent && isLatent.i && isLatent.j){
                next
            }
            if(rm.endo_endo && isEndogenous.i && isEndogenous.j){
                next
            }
            if(isExogenous.i && isExogenous.j){
                next
            }
            if(rm.latent_endo && ( (isLatent.i && isEndogenous.j) || (isEndogenous.i && isLatent.j) )){
                next
            }
            iDirectional <- (isExogenous.i+isExogenous.j)>0
            if(type == "regression" && iDirectional == FALSE){
                next
            }
            if(type == "covariance" && iDirectional == TRUE){
                next
            }
            
            if (AP[j, i] == 0){
                restricted <- rbind(restricted, c(i, j))
                directional <- c(directional, iDirectional)
            }
        }
    }
    ## *** export  
    out <- list(M.links = restricted,
                links = NULL,
                directional = directional)

    if(!is.null(restricted)){
        M.names <- cbind(rownames(AP)[restricted[,1]],
                         colnames(AP)[restricted[,2]])
        out$links <- paste0(M.names[,1], lava.options()$symbols[2-directional],M.names[,2])
        if(output == "names"){
            out$M.links <- M.names
        }
    }
    
    return(out)  
}

## * addLink
## ** doc addLink
#' @title Add a New Link Between Two Variables in a LVM
#' @rdname addLink
#' @description Generic interface to add links to \code{lvm} objects.
#' 
#' @param object a \code{lvm} object.
#' @param var1 [character or formula] the exogenous variable of the new link or a formula describing the link to be added to the lvm.
#' @param var2 [character] the endogenous variable of the new link. Disregarded if the argument \code{var1} is a formula.
#' @param all.vars [internal] a character vector containing all the variables of the \code{lvm} object.
#' @param covariance [logical] is the link is bidirectional? Ignored if one of the variables non-stochastic (e.g. exogenous variables).
#' @param warnings [logical] Should a warning be displayed when no link is added?
#' @param ... [internal] only used by the generic method and from \code{addLink.lvm.reduced} to \code{addLink.lvm}.
#'
#' @details
#' The argument \code{all.vars} is useful for \code{lvm.reduce} object where the command \code{vars(object)} does not return all variables. The command \code{vars(object, xlp = TRUE)} must be used instead.
#'
#' Arguments \code{var1} and \code{var2} are passed to \code{initVarlink}.
#'
#' @examples
#' library(lava)
#' set.seed(10)
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' m2 <- m
#' 
#' addLink(m, x1 ~ y1, covariance = FALSE)
#' addLink(m, y1 ~ x1, covariance = FALSE)
#' coef(addLink(m, y1 ~ y2, covariance = TRUE))
#' 
#' addLink(m2, "x1", "y1", covariance = FALSE)
#' addLink(m2, "y1", "x1", covariance = FALSE)
#' newM <- addLink(m, "y1", "y2", covariance = TRUE)
#' coef(newM)
#'
#' @concept setter
#' @export
`addLink` <-
    function(object, ...) UseMethod("addLink")

## ** method addLink.lvm
#' @export
#' @rdname addLink
addLink.lvm <- function(object,
                        var1,
                        var2,
                        covariance,
                        all.vars = lava::vars(object),
                        warnings = FALSE,
                        ...){

    res <- initVarLink(var1, var2, format = "list")
    var1 <- res$var1
    var2 <- res$var2
    endogenous.object <- endogenous(object)
    exogenous.object <- exogenous(object)
    latent.object <- latent(object)
    
    if(var1 %in% all.vars == FALSE){
        if(warnings){
            warning("addLink.lvm: var1 does not match any variable in object, no link is added \n",
                    "var1: ",var1,"\n")
        }
    }
    
    ####
    if(is.na(var2)){
        
        intercept(object) <- stats::as.formula(paste0("~", var1))
        
    }else{
        
        if(var1 == var2){
            if(warnings){
                warning("addLink.lvm: var1 equals var2, no link is added \n",
                        "var1/2: ",var1,"\n")
            }
        }
        
        
        ## if(var2 %in% all.vars == FALSE){
        ##     if(warnings){
        ##         warning("addLink.lvm: var2 does not match any variable in object, no link is added \n",
        ##                 "var2: ",var2,"\n")
        ##     }
        ## }
        if(var1 %in% endogenous.object && var2 %in% endogenous.object){
            if(missing(covariance)){
                covariance <- TRUE
            }else if(covariance == FALSE){
                stop("addLink.lvm: cannot add a link between two endogenous variable when argument \'covariance\' is FALSE \n")
            }                
        }
    
       
        if(covariance){
            covariance(object) <- stats::as.formula(paste(var1, var2, sep = "~"))  
        }else if(var1 %in% endogenous.object || var2 %in% exogenous.object){
            regression(object) <- stats::as.formula(paste(var1, var2,  sep = "~"))
        }else if(var2 %in% endogenous.object || var1 %in% exogenous.object){
            regression(object) <- stats::as.formula(paste(var2, var1, sep = "~"))
        }else {
            if(var1 %in% latent.object){
                regression(object) <- stats::as.formula(paste(var1, var2, sep = "~"))  
            }else if(var2 %in% latent.object){
                regression(object) <- stats::as.formula(paste(var2, var1, sep = "~"))  
            }else{
                stop("unknow configuration \n")
            }
            
        }
    }
    
    return(object)
}

## ** method addLink.lvm.reduced
#' @rdname addLink
addLink.lvm.reduced <- function(object, ...){
  return(addLink.lvm(object, all.vars = lava::vars(object, lp = FALSE, xlp = TRUE) , ...))
}

## * setLink
## ** Documentation - setLink
#' @title Set a Link to a Value
#' @name setLink
#' @description Generic interface to set a value to a link in a \code{lvm} object.
#' 
#' @param object a \code{lvm} object.
#' @param var1 [character or formula] the exogenous variable of the new link or a formula describing the link to be added to the lvm.
#' @param var2 [character] the endogenous variable of the new link. Disregarded if the argument \code{var1} is a formula.
#' @param value [numeric] the value at which the link should be set.
#' @param warnings [logical] should a warning be displayed if the link is not found in the \code{lvm} object.
#' @param ... [internal] only used by the generic method.
#'  
#' @examples
#' library(lava)
#' set.seed(10)
#' 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1 ~ y2
#' 
#' m1 <- setLink(m, y3 ~ u, value = 1)
#' estimate(m1, lava::sim(m,1e2))
#' # m1 <- setLink(m, u ~ y3, value = 1)
#' 
#' m2 <- setLink(m, y1 ~ y2, value = 0.5)
#' estimate(m2, lava::sim(m,1e2))
#'
#' @concept setter
#' @export
`setLink` <-
  function(object, ...) UseMethod("setLink")

## ** method setLink.lvm
#' @rdname setLink
#' @export
setLink.lvm <- function(object, var1, var2, value, warnings = FALSE, ...){

    res <- initVarLink(var1, var2)
    var1 <- res$var1
    var2 <- res$var2
    object.coef <- stats::coef(object)
    
  #### set the link
  if(is.na(var2)){
    intercept(object, stats::as.formula(paste0("~",var1))) <- value
  }else if(paste(var1, var2, sep = "~") %in% object.coef){
    regression(object, stats::as.formula(paste(var1,var2, sep = "~"))) <- value
  }else if(paste(var1,var2, sep = ",") %in% object.coef){
    covariance(object, stats::as.formula(paste(var1,var2, sep = "~"))) <- value
  }else if(paste(var2,var1, sep = ",") %in% object.coef){
    covariance(object, stats::as.formula(paste(var1,var2, sep = "~"))) <- value
  }else{
    if(warnings){
      warning("setLink.lvm: no link was found from var1 to var2, no link is set \n",
              "var1: ",var1,"\n",
              "var2: ",var2,"\n")
    }
  }
  
  return(object)
}




