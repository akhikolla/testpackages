## ######################################################################
#' @title Extract conditional probabilities and clique potentials from
#'     data.
#' @description Extract list of conditional probability tables and
#'     list of clique potentials from data.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name components_extract
## ######################################################################
#'
#' @details If \code{smooth} is non-zero then \code{smooth} is added
#'     to all cell counts before normalization takes place.
#' 
## #' @aliases extractCPT extractPOT extractMARG
#' 
#' @param data_ A named array or a dataframe.
#'
#' @param graph A \code{graphNEL} object or a list or formula which can be
#'     turned into a \code{graphNEL} object by calling \code{ug} or
#'     \code{dag}. For \code{extractCPT}, graph must be/define a DAG while for
#'     \code{extractPOT}, graph must be/define undirected triangulated graph.
#' 
#' @param smooth See 'details' below.
#' 
#' @return
#'   * \code{extractCPT}: A list of conditional probability tables.
#'   * \code{extractPOT}: A list of clique potentials.
#'   * \code{extractMARG}: A list of clique marginals. 
#'
## #' @details \code{extractCPT} is alias for \code{extractCPT}
## #'     \code{extractPOT} is alias for \code{extractPOT} and
## #'     \code{extractMARG} is alias for \code{extract_marg}; retained
## #'     for backward compatibility.
## #' 

#'
#' @seealso \code{\link{compileCPT}}, \code{\link{compilePOT}},
#'     \code{\link{grain}}
#'
#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords utilities
#' @examples
#' 
#' ## Extract cpts / clique potentials from data and graph
#' # specification and create network. There are different ways:
#'
#' data(lizard, package="gRbase")
#'
#' # DAG: height <- species -> diam
#' daG <- dag(~species + height:species + diam:species)
#'
#' # UG : [height:species][diam:species]
#' uG  <- ug(~height:species + diam:species)
#' 
#' pt <- extractPOT(lizard, ~height:species + diam:species) 
#' cp <- extractCPT(lizard, ~species + height:species + diam:species)
#'
#' pt
#' cp
#'
#' # Both specify the same probability distribution
#' tabListMult(pt) %>% as.data.frame.table
#' tabListMult(cp) %>% as.data.frame.table
#'
#' \dontrun{
#' # Bayesian networks can be created as
#' bn.uG   <- grain(pt)
#' bn.daG  <- grain(cp)
#'
#' # The steps above are wrapped into a convenience method which
#' # builds a network from at graph and data.
#' bn.uG   <- grain(uG, data=lizard)
#' bn.daG  <- grain(daG, data=lizard)
#' }

#' @rdname components_extract
#' @export 
extractCPT <- function(data_, graph, smooth=0){

    .is.valid.data(data_)
    if (inherits(graph, c("formula", "list"))) graph <- dag(graph)
    if (!is_dag(graph)) stop("'graph' not a DAG")
    
    vpa <- vpar(graph)
    out <- .extractCPT_primitive(data_, vpa=vpa, smooth=smooth)
    attr(out, "graph") <- graph
    class(out)         <- "cpt_rep"
    out
}


#' @export
#' @rdname components_extract
extractPOT <- function(data_, graph, smooth=0){
    
    .is.valid.data(data_)
    if (inherits(graph, c("formula", "list"))) graph <- ug(graph)    
    if (!is_tug(graph)) stop("'graph' not undirected and triangulated")

    rip_  <- rip(graph)
    out   <- .extractPOT_primitive(data_, rip_$cliques, rip_$sep, smooth=smooth)
    attr(out, "rip")   <- rip_
    attr(out, "graph") <- graph    
    class(out)         <- "pot_rep"
    out
}

#' @export 
#' @rdname components_extract
extractMARG <- function(data_, graph, smooth=0){

    .is.valid.data(data_)
    if (inherits(graph, c("formula", "list"))) graph <- ug(graph)    
    if (!is_tug(graph)) stop("'graph' not undirected and triangulated")
    
    rip_  <- rip(graph)
    out   <- .extractMARG_primitive(data_, rip_$cliques, rip_$sep, smooth=smooth)
    attr(out, "rip")   <- rip_
    attr(out, "graph") <- graph    
    class(out)         <- "marg_rep"
    out
}


#' @export 
#' @rdname components_extract
#' @param mg An object of class \code{marg_rep}
marg2pot <- function(mg){
    if (!inherits(mg, "marg_rep")) stop("'mg' not a marg_rep object\n")
    rip_ <- attr(mg, "rip")
    seps <- rip_$separators
    pt <- lapply(seq_along(rip_$cliques),
                 function(i){
                     if (length(seps[[i]]) == 0)
                         mg[[i]]
                     else
                         tabDiv0(mg[[i]], tabMarg(mg[[i]], seps[[i]]))               
                 })
    attr(pt, "rip") <- rip_
    class(pt) <- "pot_rep"
    pt
}

#' @export 
#' @rdname components_extract 
#' @param pt An object of class \code{pot_rep}
pot2marg <- function(pt){
    if (!inherits(pt, "pot_rep")) stop("'pt' not a pot_rep object\n")    
    mg <- pt
    rip_ <- attr(pt, "rip")
    seps <- rip_$separators
    par  <- rip_$parents
    
    for (i in 2:length(rip_$cliques)){
        if (par[i] > 0){
            mg[[i]] <- tabMult(mg[[i]], tabMarg(mg[[par[i]]], seps[[i]]))
        }
    }
    class(mg) <- "marg_rep"
    mg
}


## ##################################################################
##
## dot functions below here
##
## ##################################################################


.extractCPT_primitive <- function(data_, vpa, smooth=0){
        
    is.df <- is.data.frame(data_)
    out <- lapply(vpa, function(ss){.dataMarg(data_, ss, is.df)})
    
    ## FIXME : Get rid of this parray stuff (at least as a class)
    ## NOTE: Normalization takes place here
    out <- lapply(out, as.parray, normalize="first", smooth=smooth)
    
    
    chk <- unlist(lapply(out, function(zz) any(is.na(zz))))
    nnn <- names(chk)[which(chk)]
    if (length(nnn) > 0){
        cat(sprintf("NAs found in conditional probability table(s) for nodes: %s\n",
                    toString(nnn)))
        cat(sprintf("  ... consider using the smooth argument\n"))
    }
    out
}



.extractPOT_primitive <- function(data_, cliq, seps=NULL, smooth=0){        
    
    .normalize <- function(tt, sp){
        if (length(sp) > 0) tabDiv0(tt, tabMarg(tt, sp))
        else tt / sum(tt)        
    }
    
    out <- vector("list", length(cliq))
    is.df <- is.data.frame(data_)
    for ( i  in seq_along(cliq)){
        cq   <- cliq[[ i ]]
        sp   <- seps[[ i ]]
        t.cq <- .dataMarg(data_, cq, is.df) + smooth       
        ##str(list(cq=cq, sp=sp))
        out[[i]] <- .normalize(t.cq, sp)
    }
    out
}



    .extractMARG_primitive <- function(data_, cliq, seps=NULL, smooth=0){        
        out <- vector("list", length(cliq))
        is.df <- is.data.frame(data_)
        
        for (i in seq_along(cliq)){
            cq   <- cliq[[ i ]]
            t.cq <- .dataMarg(data_, cq, is.df) + smooth       
            out[[i]] <- t.cq / sum(t.cq)
        }
        out
    }


## #' @rdname components_extract
## data2cpt <- extractCPT

## #' @rdname components_extract
## data2pot <- extractPOT

## #' @rdname components_extract
## data2marg <- extractMARG





## helper function; can possibly be made faster
.dataMarg <- function(data_, cq, is.df=NULL){

    ## .dfMarg can possibly be made faster
    .dfMarg <- function(data_, cq){
        xtabs(~., data=data_[ , cq, drop=FALSE])
    }

    if (is.null(is.df))
        is.df <- is.data.frame(data_)

    if (is.df) .dfMarg(data_, cq)
    else tabMarg(data_, cq)
        
}



.is.valid.data <- function(data_){
    if (!(is.data.frame(data_) || is.named.array(data_)))
        stop("'data_' must be dataframe or array.")
}


## #' 
## #' plot(bn.uG)
## #' plot(bn.daG)
## #' 
## #' querygrain(bn.uG)
## #' querygrain(bn.daG)
## #' 
## #' # Sanity: Are the distributions identical?
## #' t1 <- querygrain(bn.uG, type="joint")
## #' t2 <- querygrain(bn.daG, type="joint")
## #' t1 %a/% t2
## #' 
## #' # At a lower level
## #' bn2.uG <- extractPOT(lizard, ~height:species + diam:species)  %>% grain
## #' bn2.daG <- extractCPT(lizard, ~species + height:species + diam:species)  %>% grain
## #' 
## #' plot(bn2.uG)
## #' plot(bn2.daG)
## #' 
## #' t1 <- querygrain(bn2.uG, type="joint")
## #' t2 <- querygrain(bn2.daG, type="joint")
## #' t1 %a/% t2
## #' 
