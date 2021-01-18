## * initVarLink (Documentation)
#' @title Normalize var1 and var2
#' @name initVarLink
#' @description Convert var1 and var2 from formula or covariance to character.
#' 
#' @param var1 [character or formula] the exogenous variable of the new link or a formula describing the link.
#' @param var2 [character] the endogenous variable of the new link.
#' Disregarded if the argument \code{var1} is a formula.
#' @param rep.var1 [logical] should var1 be duplicated to match var2 length.
#' Only active if \code{format = "list"}.
#' @param format [character] should the name of the variable be returned (\code{format = "list"}),
#' a vector of character formula (\code{format = "txt.formula"}),
#' or a list of formula (\code{format = "formula"}).
#' @param Slink [character] the symbol for regression link.
#' @param Scov [character] the symbol for covariance link.
#' @param ... argument to be passed to \code{initVarLink}.
#'
#' @return See argument \code{format}.
#' 
#' @examples
#' initVarLink(y ~ x1)
#' initVarLink("y ~ x1")
#' initVarLink(y ~ x1 + x2)
#' initVarLink("y ~ x1 + x2")
#' initVarLink(y ~ x1 + x2, rep.var1 = TRUE)
#' initVarLink(y ~ x1 + x2, rep.var1 = TRUE, format = "formula")
#' initVarLink(y ~ x1 + x2, rep.var1 = TRUE, format = "txt.formula")
#' initVarLink("y", "x1", format = "formula")
#'
#' initVarLink("y ~ x1:0|1")
#'
#' initVarLinks(y ~ x1)
#' initVarLinks("y ~ x1")
#' initVarLinks(c("y ~ x1","y~ x2"))
#' initVarLinks(c(y ~ x1,y ~ x2))
#' initVarLinks(c("y ~ x1","y ~ x2"), format = "formula")
#' initVarLinks(c(y ~ x1,y ~ x2), format = "formula")
#' initVarLinks(c("y ~ x1","y~ x2"), format = "txt.formula")
#' initVarLinks(c(y ~ x1,y ~ x2), format = "txt.formula")

## * initVarLink
#' @rdname initVarLink
#' @export
initVarLink <- function(var1, var2, rep.var1 = FALSE, format = "list",
                         Slink = c(lava.options()$symbols[1],"~"),
                         Scov = lava.options()$symbols[2]){

    format <- match.arg(format, c("list","txt.formula","formula"))
    test.formula <- (class(var1) == "formula")
    test.covariance <- sapply(Scov,grepl,x=var1,fixed=TRUE)
    test.regression <- sapply(Slink,grepl,x=var1,fixed=TRUE)

    if(missing(var2)){
        if(test.formula){
            var2 <- selectRegressor(var1, format = "vars")
            var1 <- selectResponse(var1, format = "vars")
            sep <- if(format == "formula"){"~"}else{Slink}
        }else if(any(test.covariance)){ ## covariance
            Scov <- Scov[test.covariance][1]
            varSplit <- strsplit(var1, split = Scov)[[1]]
            var1 <- trimws(varSplit[1])
            var2 <- trimws(varSplit[2])
            sep <- if(format == "formula"){"~"}else{Scov}
        } else if(any(test.regression)){ ## regression
            Slink <- Slink[test.regression][1]
            varSplit <- strsplit(var1, split = Slink)[[1]]
            var1 <- trimws(varSplit[1])
            var2 <- trimws(varSplit[2])
            sep <- if(format == "formula"){"~"}else{Slink}
        } else {
            var1 <- var1
            var2 <- NA
        }
    }else{

        if(!is.character(var1) || !is.character(var2)){
            stop("\'var1\' and \'var2\' must be characters when both are specified \n")
        }
        sep <- if(format == "formula"){"~"}else{Slink}        
    }
  
  
#### convert to format
    if(format == "formula"){
        n.var2 <- length(var2)
        var1 <- rep(var1, times = n.var2)
        res <- sapply(1:n.var2, function(i){
            stats::as.formula(paste(var1[i], var2[i], sep = sep))
        })
    
  }else if(format == "txt.formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){paste(var1[i], var2[i], sep = sep)})
    
  }else if(format == "list"){
    if(rep.var1 && !missing(var1)){var1 <- rep(var1, length(var2))}
    res <- list(var1 = var1,
                var2 = if(!missing(var2)){var2}else{NULL} 
    )
  }
 
  ## export 
  return(res)
}

## * initVarLinks
#' @rdname initVarLink
#' @export
initVarLinks <- function(var1, format = "list",...){
        
    if("formula" %in% class(var1)){
        res <- initVarLink(var1, rep.var1 = TRUE, format = format,
                           ...)
    }else {
        res <- sapply(var1, function(x){
            initVarLink(x, rep.var1 = TRUE, format = format,
                        ...)
        })
        if(format == "list"){
            res <- list(var1 = unname(unlist(res["var1",])),
                        var2 = unname(unlist(res["var2",])))
        }else{
            res <- unname(unlist(res))
        }
    }

    return(res)
    
}

