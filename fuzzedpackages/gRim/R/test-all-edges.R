
## Known issues
## It is not tested whether edges in edgeList are in the model.
## Should check if edgeList is NULL

#' @title Test edges in graphical models with p-value/AIC value
#' 
#' @description Test edges in graphical models with p-value/AIC
#'     value. The models must be \code{iModel}s.
#'
#' @name test-edges
#' 
#' @aliases testEdges testInEdges testOutEdges  
#'
#' @param object An \code{iModel} model object
#' @param edgeMAT A `p * 2` matrix with edges
#' @param criterion Either \code{"aic"} or \code{"test"} (for
#'     significance test)
#'
#' @details
#'
#' * testIn: Function which tests whether each edge in "edgeList" can
#'     be delete from model "object"
#'
#' * testOut: Is similar but in the other direction.
#' 
#' @param k Penalty term when \code{criterion="aic"}. Only k=2 gives
#'     genuine AIC.
#' @param ingraph If TRUE, edges in graph are tested; if FALSE, edges
#'     not in graph are tested.
#' @param alpha Critical value for deeming an edge to be significant/
#'     insignificant. When \code{criterion="aic"}, \code{alpha}
#'     defaults to 0; when \code{criterion="test"}, \code{alpha}
#'     defaults to 0.05.
#' @param headlong If TRUE then testing will stop once a model
#'     improvement has been found.
#' @param details Controls the level of printing on the screen.
#' @param \dots Further arguments to be passed on to \code{testdelete}
#'     (for \code{testInEdges}) and \code{testadd} (for
#'     \code{testOutEdges}).
#' @return A dataframe with test statistics (p-value or change in
#'     AIC), edges and logical telling if the edge can be deleted.
#' 
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{getEdges}}, \code{\link{testadd}},
#'     \code{\link{testdelete}}
#' @keywords models htest
#' @examples
#' 
#' data(math)
#' cm1 <- cmod(~me:ve + ve:al + al:an, data=math)
#' testEdges(cm1, ingraph=TRUE)
#' testEdges(cm1, ingraph=FALSE)
#' ## Same as
#' # testInEdges(cm1)
#' # testOutEdges(cm) 

#' @export
#' @rdname test-edges
testEdges <- function(object, edgeMAT=NULL, ingraph=TRUE, criterion="aic", k=2, alpha=NULL,
                      headlong=FALSE, details=1, ...){
  UseMethod("testEdges")
}

#' @export
testEdges.iModel <- function(object, edgeMAT=NULL, ingraph=TRUE, criterion="aic", k=2, alpha=NULL,
                      headlong=FALSE, details=1, ...){
  
    cl <- match.call()
    if (ingraph) cl[[1]] <- as.name("testInEdges")
    else cl[[1]] <- as.name("testOutEdges")
    eval(cl)
}

 
## FIXME : Need Generic testInEdges / testOutEdges

#' @export
#' @rdname test-edges
testInEdges <- function(object, edgeMAT=NULL, criterion="aic", k=2, alpha=NULL, headlong=FALSE, details=1, ...){

    criterion <- match.arg(criterion, c("aic", "test"))
    
    vn   <- getmi(object, "varNames")
    amat <- .glist2adjMAT(getmi(object, "glist"), vn=vn)

    if (is.null(alpha))
        alpha <- if (criterion=="aic") 0 else 0.05
    
    ff <- drop_func(criterion)
    
    
    if (is.null(edgeMAT)) edgeMAT <- getInEdgesMAT(amat)
    if (nrow(edgeMAT)==0) stop("There are no edges to test...\n")

    head.str <- if (headlong) "headlong" else "all" 
    .test_in_edges(object, edgeMAT, comp.op=ff$comp.op, crit.str=ff$crit.str, alpha=alpha, k=k, search=head.str, amat=amat, vn=vn, ...)
}

#' @export
#' @rdname test-edges
testOutEdges <- function(object, edgeMAT=NULL, criterion="aic", k=2, alpha=NULL, headlong=FALSE, details=1,...){

    criterion <- match.arg(criterion, c("aic", "test"))

    vn   <- getmi(object, "varNames")
    amat <- .glist2adjMAT(getmi(object, "glist"), vn=vn)
            
    if (is.null(alpha))
        alpha <- if (criterion=="aic") 0 else 0.05
    
    ff <- add_func(criterion)

        
    if (is.null(edgeMAT)) edgeMAT <- getOutEdgesMAT(amat)
    if (nrow(edgeMAT) == 0) stop("There are no missing edges to test...\n")

    head.str <- if (headlong) "headlong" else "all" 
    .test_out_edges(object, edgeMAT, comp.op=ff$comp.op, crit.str=ff$crit.str, alpha=alpha, k=k, search=head.str  ,amat=amat, vn=vn, ...)

}



.test_in_edges <- function(object, edgeMAT, comp.op=`<`, crit.str="aic", alpha=0, k=2, search="all", amat, vn, ...){

    headlong <- identical(search, "headlong")
    testMAT <- matrix(0, nrow=nrow(edgeMAT), ncol=4)
    colnames(testMAT) <- c("statistic", "df", "p.value", "aic")
    if (headlong)
        perm <- sample(nrow(edgeMAT))              
    else
        indic <- rep.int(0, nrow(edgeMAT))
    
    for (i in seq_len(nrow(edgeMAT))){
        if (headlong)
            uv          <- edgeMAT[perm[i], ]              ## The headlong part
        else 
            uv          <- edgeMAT[i, ]
        edgeTest    <- testdelete(object, uv, k=k, amat=amat, ...) ## amat       
        testMAT[i,] <- as.numeric(edgeTest[c("statistic", "df", "p.value", "aic")])
        curr.stat   <- edgeTest[[crit.str]]
        if (headlong){
            if (comp.op(curr.stat, alpha))  break  ## The headlong part
        } else {
            if (comp.op(curr.stat, alpha)) indic[i] <- 1             
        }
    }

    if (headlong){
        ans <- cbind(
            as.data.frame(testMAT[1:i, , drop=FALSE]),
            as.data.frame(edgeMAT[perm[1:i], , drop=FALSE],
                          stringsAsFactors=FALSE), ## The headlong part
            action=c(rep("-", i - 1), "+"))
    } else {
        ans <- cbind(
            as.data.frame(testMAT),
            as.data.frame(edgeMAT, stringsAsFactors=FALSE),
            action=c("-", "+")[indic + 1])
        }
    ans
}

.test_out_edges <- function(object, edgeMAT, comp.op=`>`, crit.str="aic", alpha=0, k=2, search="all", amat, ...)
{

    headlong <- identical(search, "headlong")
    testMAT <- matrix(0, nrow=nrow(edgeMAT), ncol=4)
    colnames(testMAT) <- c("statistic", "df", "p.value", "aic")
    
    if (headlong)
        perm <- sample(nrow(edgeMAT))
    else
        indic <- rep.int(0, nrow(edgeMAT))
    
    for (i in seq_len(nrow(edgeMAT))){
        if (headlong)
            uv          <- edgeMAT[perm[i], ]
        else 
            uv          <- edgeMAT[i,]
        edgeTest    <- testadd(object, uv, k=k, amat=amat, ...)
        testMAT[i,] <- as.numeric(edgeTest[c("statistic", "df", "p.value", "aic")])
        curr.stat   <- edgeTest[[crit.str]]
        if (headlong){
            if (comp.op( curr.stat, alpha)) break
        } else {
            if (comp.op(curr.stat, alpha)) indic[i] <- 1
        }
    }

    if (headlong){
        ans <- cbind(
            as.data.frame(testMAT[1:i, , drop=FALSE]),
            as.data.frame(edgeMAT[perm[1:i], , drop=FALSE],
                          stringsAsFactors=FALSE),
            action=c(rep("-", i - 1), "+"))
    } else {        
        ans <- cbind(
            as.data.frame(testMAT),
            as.data.frame(edgeMAT, stringsAsFactors=FALSE),
            action=c("-", "+")[indic + 1])
    }
    ans
}  




