##################################################################
####
#### Update iModel
####
##################################################################

## cat("update.iModel\n"); cat("items: "); print( items )
## print(class(object))
## print( object$glist )
## print( object$varNames )

## items <<- items
## glist <<- object$glist
## vn <<- object$varNames

#' @export
update.iModel <- function(object, items, fit=TRUE, details=0, ...){

    ## cat("update.iModel: before \n"); print(.glist(object))
    glist            <- modify_glist(.glist(object), items)
    .glist(object)   <- glist
    ## cat("update.iModel: after \n"); print(.glist(object))
    
    switch(class(object)[1],
           "dModel"={
               upd <- .dModel_finalize(glist, object$varNames)
           },
           "cModel"={
               upd <- .cModel_finalize(glist, object$varNames)
           },
           "mModel"={
               upd <- .mModel_finalize(glist, object$varNames, object$datainfo)
           } )

    ##object[ names(upd) ] <- upd
    object$modelinfo <- upd    
    if (fit) fit(object) else object
}


.update.default <- function (object, formula., ..., evaluate = TRUE)
{
    if (is.null(call <- getCall(object)))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
        call$formula <- update.formula(formula(object), formula.)
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate)
        eval(call, parent.frame())
    else call
}


triangulate.dModel <- function(object, ...){
    cl <- getCall(object)
    cq <- getCliques(triangulate(.glist2adjMAT(terms(object))))
    ff <- list2rhsf(cq)
    cl$formula <- ff
    eval(cl, parent.frame())
}


##################################################################
####
#### Update generator list by adding/deleting edges and terms
####
#### FIXME: Perhaps add... should check if ee/term is in the list
#### in which case a special value should be returned
####
##################################################################

### modify_glist is the workhorse.
### Updates an entire glist with the elements (edges, terms) in items
###


#' @title Modify generating class for a graphical/hierarchical model
#' 
#' @description Modify generating class for a graphical/hierarchical model by 1)
#'     adding edges, 2) deleting edges, 3) adding terms and 4) deleting terms.
#'
#' @name modify_glist
#' 
#' @details
#' 
#' The \code{items} is a list with named entries as \code{list(add.edge=,
#' drop.edge=, add.term=, drop.term=)}
#' 
#' Not all entries need to be in the list. The corresponding actions are
#' carried out in the order in which they appear in the list.
#' 
#' See section 'examples' below for examples.
#' 
#' Notice that the operations do not in general commute: Adding an edge which
#' is already in a generating class and then removing the edge again does not
#' give the original generating class.
#' 
#' @param glist A list of vectors where each vector is a generator of the model.
#' @param items A list with edges / terms to be added and deleted. See section
#'     'details' below.
#' @param details Control the amount of output (for debugging purposes).
#' @return A generating class for the modified model. The elements of the list
#'     are character vectors.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{cmod}}, \code{\link{dmod}}, \code{\link{mmod}}
#' @keywords utilities
#' @examples
#' 
#' glist <- list(c(1, 2, 3), c(2, 3, 4))
#' 
#' ## Add edges
#' modify_glist(glist, items=list(add.edge=c(1, 4)))
#' modify_glist(glist, items=list(add.edge=~1:4))
#' 
#' ## Add terms
#' modify_glist(glist, items=list(add.term=c(1, 4)))
#' modify_glist(glist, items=list(add.term=~1:4))
#' 
#' ## Notice: Only the first term is added as the second is already 
#' ## in the model.
#' modify_glist(glist, items=list(add.term=list(c(1, 4), c(1, 3))))
#' modify_glist(glist, items=list(add.term=~1:4 + 1:3))
#' 
#' ## Notice: Operations are carried out in the order given in the
#' ## items list and hence we get different results: 
#' modify_glist(glist, items=list(drop.edge=c(1, 4), add.edge=c(1, 4)))
#' modify_glist(glist, items=list(add.edge=c(1, 4), drop.edge=c(1, 4)))
#' 
#' @export modify_glist
#' 

modify_glist <- function(glist, items, details=0){

    ## cat("modify_glist items (before): "); str(items)
    glist   <- lapply(glist, as.character)
    ## Here; whatever the input format is "taken apart into lists":
    
    action <- names( items )
    
  items <- lapply( items, .doInput )
  ## cat("modify_glist items (after .doInput): "); print(items)
  names( items ) <- action

  items   <- .parse.change.list(items, details)

  ## cat("modify_glist items (after parse change) : "); print(items)
  for (i in seq_along(items)){
    curr.action  <- action[ i ]
    curr.item    <- items[[ i ]]
    glist        <- .modify_glistPrim(glist, curr.action, curr.item, details)
  }
  glist
}

### Updates a glist (generating class) with the elements in
### curr.item. These can be of the type curr.action where valid
### choices are add.edge, drop.edge, add.term and drop.term
###
.modify_glistPrim <- function(glist, curr.action, curr.item, details=0){
    fname <- paste(".",curr.action,"_glist",sep="")
    #cat(sprintf("fname: %s\n", fname))
    ## .infoPrint(details,1,cat(sprintf("action: %s \n", curr.action)))

    for (k in seq_along(curr.item)){
        curr <- curr.item[[ k ]]
        ##cat(sprintf("action: %s item: %s\n", fname, paste(curr, collapse=" ")))
        glist <- do.call(fname, list(glist, curr))
        ##print(glist)
    }
    glist
}



.add.edge_glist <- function(glist, ee){
    extra <- list()
    count <- 1
    ss <- seq_along(glist)
    for (i in ss){
        if (ee[1] %in% glist[[i]]){
            for (j in ss){
                if (ee[2] %in% glist[[ j ]]){
                    zz <- intersectPrim(glist[[i]], glist[[ j ]])
                    extra[[ count ]] <- unique.default(c(ee, zz))
                    count <- count + 1
                }
            }
        }
    }
    remove_redundant( c(glist, extra) )
}






.doInput <- function( e ){
    cls <- class(e)
    if (cls == "data.frame" || cls == "matrix"){
        e <- as.matrix( e )
        if (ncol( e ) != 2)
            stop("Must have dimension p x 2\n")
        e <- rowmat2list( e )
    }
    rhsf2list( e )
}

### e1 <- c(1,4)
### e2 <- c(2,4)
### e3 <- ~1:4
### e4 <- ~1:4+2:4
### e5 <- rbind(e1,e2)
### e6 <- as.data.frame(e5)
### e7 <- list(e1, e2)
##
### .doInput( e1 )
### .doInput( e2 )
### .doInput( e3 )
### .doInput( e4 )
### .doInput( e5 )
### .doInput( e6 )
### .doInput( e7 )
##

.drop.edge_glist <- function(glist, ee){
  .drop.term_glist(glist, ee)
}

.add.term_glist <- function(glist, term){
    ##if (isin( glist, term ))
    if (is_inset(term, glist))
        glist
    else
        remove_redundant( c(list(term), glist) )
}

.drop.term_glist <- function(glist, term){
    #cat(".drop.term_glist\n")  #print(glist); print(term)
    extra <- list()
    count <- 1
    changed <- rep(0, length(glist))

    ## If the i'th generator 'gen.i' contains 'term' then gen.i will
    ## be marked with a 1, and otherwise with a 0.

    ## If gen.i and term are identical, then gen.i will be expanded to
    ## all terms one order lower; these will be included in the output
    ## whereas gen.i itself will not.

    for (i in seq_along(glist)){
        gen <- glist[[ i ]]
        ## cat("term:\n"); print(term); cat("gen:\n"); print(gen)
        if (subsetof( term, gen )){
            ##cat("term is subset of gen...\n")
            changed[ i ] <- 1
            lower <- combn_prim(gen, length(gen)-1, simplify=FALSE)
            ##cat("lower:\n"); print(lower)
            if ( length(term) == length(gen) ){
                extra[[ count ]] <- lower
            } else {
                keep   <- unlist(lapply(lower, function(s) !subsetof(term, s)), use.names=FALSE)
                ##print(keep)
                lower <- lower[ keep ]
                extra[[ count ]] <- lower
            }
            count <- count + 1
        }
    }

    glist.new <- c(glist[ changed==0 ], unlist(extra, recursive=FALSE, use.names=FALSE))
    remove_redundant( glist.new )
}


.aedge_glist <- .add.edge_glist
.dedge_glist <- .drop.edge_glist
.aterm_glist <- .add.term_glist
.dterm_glist <- .drop.term_glist


.preprocess <- function(e){
    if( class(e)=="data.frame" ){
        e <- as.matrix( e )
    }
    if (class(e)=="matrix"){
        if (ncol(e) !=2 )
            stop("matrix of edges must have two columns\n")
        e <- rowmat2list( e )
    }
    if (!is.list(e))
        e <- list(e)
    e
}

addEdge_glist <- function(glist, e){
    e <- .preprocess( e )
    for (i in seq_along(e))
        glist <- .add.edge_glist( glist, e[[i]] )
    glist
}

dropEdge_glist <- function(glist, e){
    e <- .preprocess( e )
    for (i in seq_along(e))
        glist <- .drop.edge_glist( glist, e[[i]] )
    glist
}

addTerm_glist <- function(glist, e){
    if (!is.list(e))
        e <- list(e)
    for (i in seq_along(e))
        glist <- .add.term_glist( glist, e[[i]] )
    glist
}

dropTerm_glist <- function(glist, e){
    if (!is.list(e))
        e <- list(e)

    for (i in seq_along(e))
        glist <- .drop.term_glist( glist, e[[i]] )
    glist
}





### An ad.list can have elements with names add.edge, drop.edge,
### add.term and drop.term These can be formulae, and
### .parse.change.list will transform these into appropritate lists.
###
.parse.change.list <- function(items,details=0){
    ##cat("In function: .parse.change.list:\n")
    .foo <- function(curr.action, curr.item){
        switch(curr.action,
               "add.edge"=,
               "drop.edge"=,
               "aedge"=,
               "dedge"={
                   zzz <- unlist(lapply(rhsf2list(curr.item), names2pairs),
                                 recursive=FALSE, use.names=FALSE)
               },
               "add.term"=,
               "drop.term"=,
               "aterm"=,
               "dterm"={
                   zzz <- rhsf2list(curr.item)
               })
        zzz
    }
    nam   <- names(items)
    valid <- c("add.edge","drop.edge","add.term","drop.term",
               "aedge","dedge","aterm","dterm")

    for (i in 1:length(items)){
        curr.action <- nam[i]
        aaa <- match(curr.action, valid)
        if (is.na(aaa))
            stop(sprintf("Item %i has name '%s' which is not valid\n",i, curr.action))
        curr.item <- items[[i]]
        .infoPrint(details,1, cat(sprintf("parsing action %s on item %s\n", curr.action, toString(curr.item))))
        items[[i]] <- .foo(curr.action, curr.item)
    }
    ##cat("On exit:\n"); print(items)
    items
}

















