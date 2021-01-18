## Turn a formula into a list of generators (glist)
##
#' @title Parse graphical model formula
#'
#' @description Parse graphical model formula to internal representation
#'
#' @param formula A right hand sided formula or a list.
#' @param varnames Specification of the variables.
#' @param marginal Possible specification of marginal (a set of
#'     variables); useful in connection with model specification
#'     shortcuts.
#' @param interactions The maximum order of interactions allowed;
#'     useful in connection with model specification shortcuts.
#'
#' @examples
#' vn <- c("me", "ve", "al", "an", "st")
#'
#' form1 <- ~me:ve:al + ve:al + an
#' form2 <- ~me:ve:al + ve:al + s
#' form3 <- ~me:ve:al + ve:al + anaba
#' parse_gm_formula(form1, varnames=vn)
#' parse_gm_formula(form2, varnames=vn)
#' ## parse_gm_formula(form3, varnames=vn)

#' parse_gm_formula(form1)
#' parse_gm_formula(form2)
#' parse_gm_formula(form3)
#'
#' ## parse_gm_formula(~.^1)
#' ## parse_gm_formula(~.^.)
#'
#' parse_gm_formula(~.^1, varnames=vn)
#' parse_gm_formula(~.^., varnames=vn)
#' parse_gm_formula(~.^., varnames=vn, interactions=3)
#' 
#' vn2 <- vn[1:3]
#' ## parse_gm_formula(form1, varnames=vn, marginal=vn2)
#' ## parse_gm_formula(form2, varnames=vn, marginal=vn2)
#' ## parse_gm_formula(form3, varnames=vn, marginal=vn2)
#' parse_gm_formula(~.^1, varnames=vn, marginal=vn2)
#' parse_gm_formula(~.^., varnames=vn, marginal=vn2)
#' 


#' @export
parse_gm_formula <- function (formula, varnames=NULL, marginal=NULL, interactions=NULL)
{

    varnames <- if (length(marginal) > 0) marginal else varnames 

    if (!is.atomic(varnames)) stop("'varnames' must be atomic\n")
    
    if (!inherits(formula, c("formula", "list"))) stop("Invalid formula specification")
        
    switch(class(formula),
           "formula"={
               glist <- .do.formula(formula, varnames)               
           },
           "list"={
               glist <- formula
           })

    
    glist <- remove_redundant(glist)    
    glist <- .check.glist(glist, varnames)  
    
    if (!is.null(interactions))
        glist <- .set.interactions(glist, interactions)
    
    value <- list(glist    = glist,
                  varNames = unique(unlist(glist)))
    value
}


.do.formula <- function(formula, varnames=NULL){

    ## 'formula' is a right hand sided formula.
    pow <- .extract.power(formula)
    ##cat(sprintf("A formula is given; power=%d\n", pow))

    if (is.na(pow)){
        return(rhsFormula2list(formula)) ##cat("A proper formula\n")
    }

    if (is.null(varnames))
        stop("'formula' is special, and 'varnames' is needed\n")
    
    if (identical(pow, -1L)){
        glist <- list(varnames)         ##cat("The saturated model\n")
    } else {
        if (identical(pow, 1L)){
            glist <- as.list(varnames)  ##cat("The independence model\n")
        } else {
            pow   <- min(c(pow, length(varnames)))
            glist <- combn_prim(varnames, pow, simplify=FALSE)
        }               
    }
    glist     
}

.check.glist <- function(glist, varnames){

    if (is.null(varnames)) return(glist)
        
    ## It is allowed to abbreviate variable names; they are matched
    ## against names in varnames
    if (any(is.na(pmatch(unlist(glist), varnames, duplicates.ok = TRUE)))) 
        stop("An invalid variable specification has been found\n")
    glist <- lapply(glist, function(x) {
        ii <- pmatch(x, varnames)
        varnames[ii]
    })
    
    modnames <- unique.default(unlist(glist))
    
    if (any(is.na(match(modnames, varnames))))
        stop("Variables in model not contained in the variable set. Perhaps a problem with 'marginal'?")
    
    glist
}


.set.interactions <- function(glist, interactions){
  zz <- lapply(glist, function(ss){
    if (length(ss) <= interactions){
        list(ss)
    } else {
        combn_prim(ss, interactions, simplify=FALSE)
    }
  })
  remove_redundant(unlist(zz, recursive=FALSE))
}



## .extract.power returns NA if the formula is "standard" or an
## integer if it has a "hat"
.extract.power <- function(form){

    if (length(form) == 3) stop("only rhs formula allowed")
    
    form.str <- deparse(form[[2]])

    has.hat <- length(grep("^\\.\\^", form.str)) > 0

    if (!has.hat) return (NA)

    rest <- gsub("\\.\\^", "", form.str)
    pp <- strsplit(rest, " ")[[1]][1]
    pow <- if (identical(pp, ".")) -1L
           else as.integer(pp)
    
    if (identical(pow, 0L))
        pow <- 1L
    
    pow
}

