#' @title Conditional probability tables based on logical dependencies
#' 
#' @description Generate conditional probability tables based on the
#'     logical expressions AND and OR.
#'
#' @name logical
#' 
#' @details Regarding the form of the argument \code{vpa}: To specify
#'     \eqn{P(a|b,c)} one may write \code{~a|b+c} or \code{~a+b+c} or
#'     \code{~a|b:c} or \code{~a:b:c} or \code{c("a","b","c")}.
#'     Internally, the last form is used. Notice that the \code{+} and
#'     \code{:} operator are used as separators only. The order of the
#'     variables is important so \code{+} and \code{:} DO NOT commute.
#' 
#' @aliases andtable ortable
#' @param vpa Node and two parents; as a formula or a character
#'     vector.
#' @param op A logical operator.
#' @param levels The levels (or rather labels) of v, see 'examples'
#'     below.
#'
#' @note \code{andtable} and \code{ortable} are aliases for
#'     \code{andtab} and \code{ortab} and are kept for backward
#'     compatibility.
#' @return An array.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{cptable}}
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities
#'
#' @rdname logical
#'
#' @examples
#'
#' ## Logical OR:
#'
#' ## A variable v is TRUE if either of its parents pa1 and pa2 are TRUE:
#' ortab( c("v", "pa1", "pa2") ) %>% ftable(row.vars="v")
#' ## TRUE and FALSE can be recoded to e.g. yes and no:
#' ortab( c("v", "pa1", "pa2"), levels=c("yes", "no") ) %>% ftable(row.vars="v")
#'
#' ## Logical AND:
#'
#' ## Same story here:
#' andtab(c("v", "pa1", "pa2") ) %>% ftable(row.vars="v")
#' andtab(c("v", "pa1", "pa2"), levels=c("yes", "no") ) %>% ftable(row.vars="v")
#'
#' ## Combined approach
#' 
#' booltab(c("v", "pa1", "pa2"), op=`&`) %>% ftable(row.vars="v") ## AND
#' booltab(c("v", "pa1", "pa2"), op=`|`) %>% ftable(row.vars="v") ## OR
#'
#' booltab(~v + pa1 + pa2, op=`&`) %>% ftable(row.vars="v") ## AND
#' booltab(~v + pa1 + pa2, op=`|`) %>% ftable(row.vars="v") ## OR
#'
## ## We use our own operator, for example "exclusive or" which is
## ## TRUE only if one but not both arguments are TRUE.

## xor <- function(e1,e2){!mapply(all, e1, e2) & mapply(any, e1, e2)}
## booltab(c("v", "pa1", "pa2"), op=xor) %>% ftable(row.vars="v") ## XOR
## 

#' @export
booltab <- function(vpa, levels=c(TRUE, FALSE), op=`&`){
    pa=c(TRUE,FALSE)
    vpa <- c(.formula2char(vpa))
    if (length(vpa) != 3)
        stop("Must have exactly two parents!")

    g <- expand.grid(pa1=pa, pa2=pa)
    ch <- op( g[[1]], g[[2]] )
    out <- c(ch, !ch) * 1
    
    uni <- list( levels, levels, levels )
    names(uni) <- vpa[c(2, 3, 1)]
    out <- array(out, dim=c(2, 2, 2), dimnames=uni)
    out <- aperm(out, c(3, 1, 2))
    out  
}

#' @export
#' @rdname logical 
andtab <- function(vpa, levels=c(TRUE,FALSE) ){
    booltab(vpa, levels, op=`&`)
}

#' @export
#' @rdname logical
ortab <- function(vpa, levels=c(TRUE,FALSE) ){
    booltab(vpa, levels, op=`|`)
}

## For backward compatibility

#' @export
#' @rdname logical
andtable <- andtab
#' @export
#' @rdname logical
ortable <- ortab

#' @title Mendelian segregation
#'
#' @description Generate conditional probability table for mendelian
#'     segregation.
#' @param allele A character vector.
#' @param names  Names of columns in dataframe.
#' @note No error checking at all on the input.
#' @examples
#' ## Inheritance of the alleles "y" and "g"
#' 
#' men <- mendel(c("y","g"), names=c("ch", "fa", "mo"))
#' men
#' 
#' @export
mendel <- function(allele, names=c("child", "father", "mother")){
    fa   <- unique(lapply(unlist(
        lapply(allele, function(s) {lapply(allele, c, s)}), recursive=FALSE),
                          sort))
    comb <- do.call(rbind, fa)[,2:1]
    gts  <- paste(comb[,1], comb[,2], sep="")
    tab  <- expand.grid(child=gts, mother = gts, father = gts,
                       stringsAsFactors=FALSE)
    tab$prob <- mapply(.tgprob, tab$child, tab$mother, tab$father)
    names(tab)[1:3] <- names
    tab
    }


## Thereses goodie!!
.tgprob <- function(child, mother, father){
    child <- strsplit(child, "")[[1]]
    mother <- strsplit(mother, "")[[1]]
    father <- strsplit(father, "")[[1]]
    ## Probability of inheriting allele a from genotype gt
    P <- function(a, gt)((a == gt[1]) + (a == gt[2]))/2
    if (child[1] != child[2]){
        P(child[1], mother) * P(child[2], father) +
            P(child[1], father) * P(child[2], mother)
    } else {
        P(child[1], mother) * P(child[2], father)
    }
}




