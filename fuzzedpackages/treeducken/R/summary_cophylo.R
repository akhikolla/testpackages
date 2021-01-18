#' @param cophy Object of class `cophy`
#' @describeIn print.cophy Returns host tree of a cophylogenetic set
#' @export
host_tree <- function(cophy) UseMethod("host_tree")
#' @param cophy Object of class `cophy`
#' @describeIn print.cophy Returns host tree  of a cophylogenetic set
#' @export
host_tree.cophy <- function(cophy) cophy$host_tree
#' @param cophy Object of class multiCophy
#' @describeIn print.cophy Returns host tree of each member of a list
#'  of cophylogenetic sets
#' @export
host_tree.multiCophy <- function(cophy) {
    x <- unlist(cophy, recursive = FALSE)
    x[seq(from = 1, to = length(x), by = 4)]
}
#' @param cophy Object of class `cophy`
#' @describeIn print.cophy Returns symb tree of a cophylogenetic set
#' @export
symb_tree <- function(cophy) UseMethod("symb_tree")
#' @param cophy Object of class `cophy`
#' @describeIn print.cophy Returns symb tree of a cophylogenetic set
#' @export
symb_tree.cophy <- function(cophy) cophy$symb_tree

#' @param cophy Object of class `multiCophy`
#' @describeIn print.cophy Returns symb tree of each member of a
#'  list of cophylogenetic sets
#' @export
symb_tree.multiCophy <- function(cophy) {
    x <- unlist(cophy, recursive = FALSE)
    x[seq(from = 2, to = length(x), by = 4)]
}
#' @param cophy Object of class `cophy`
#' @describeIn print.cophy Returns association matrix of a cophylogenetic set
#' @export
association_mat <- function(cophy) UseMethod("association_mat")
#' @param cophy Object of class `cophy`
#' @describeIn print.cophy Returns association matrix of a cophylogenetic set
#' @export
association_mat.cophy <- function(cophy) cophy$association_mat

#' @param cophy Object of class `multiCophy`
#' @describeIn print.cophy Returns association matrix for each member of a
#'  list of cophylogenetic sets
#' @export
association_mat.multiCophy <- function(cophy) {
    x <- unlist(cophy, recursive = FALSE)
    x[seq(from = 2, to = length(x), by = 3)]
}

#' @param cophy Object of class `cophy`
#' @describeIn summary.cophy Returns event history of a cophylogenetic set
#' @export
event_history <- function(cophy) UseMethod("event_history")
#' @param cophy Cophylogenetic set
#' @describeIn summary.cophy Returns event history of a cophylogenetic set
#' @export
event_history.cophy <- function(cophy) cophy$event_history

#' @param cophy Object of class `multiCophy`
#' @describeIn print.cophy Returns event_history for each member of a list
#'  of cophylogenetic sets
#' @export
event_history.multiCophy <- function(cophy) {
    x <- unlist(cophy, recursive = FALSE)
    x[seq(from = 4, to = length(x), by = 4)]
}

#' Summarize a cophylogenetic set
#' @description Several utility functions for cophylogenetic set summarization.
#' Including functions for printing an entire summary, and a summary of each
#' part: host_tree, symb_tree, association_mat, and event_history.
#' 
#' @param object An object of class `cophy`
#' @param ... Further arguments used in generic classes
#'
#' @details The summary for a cophylogenetic set outputs a summary of the host
#'  tree and the symbiont tree.
#' The number of rows and columns of the association matrix, and a summary of
#'  the event_history.
#' 
#' @return Summary returns NULL.
#' @author Wade Dismukes and Emmanuel Paradis
#' @seealso sim_cophylo_bdp, summary for the generic, multiCophy, c.cophy
#' @examples
#' h_lambda <- 1.0
#' h_mu <- 0.3
#' c_lambda <- 0.0
#' s_lambda <- 1.0
#' s_mu <- 0.3
#' s_her <- 0.0
#' host_symb_sets <- sim_cophylo_bdp(hbr = h_lambda,
#'                                   hdr = h_mu,
#'                                   sbr = s_lambda,
#'                                   cosp_rate = c_lambda,
#'                                   sdr = s_mu,
#'                                   host_exp_rate = s_her,
#'                                   time_to_sim = 1.0,
#'                                   numbsim = 1)
#' summary(host_symb_sets[[1]])
#' @export
summary.cophy <- function(object, ...) {
        cat("Null tree set.")
    cat("\nSet of host and symbiont tree:", deparse(substitute(object)), "\n\n")
    # below is a modified version of ape's summary.phylo function
    if (is.null(object$host_tree))
        cat("\n Host tree is not a proper tree.\n")
    else {
        ht_tips <- length(object$host_tree$tip.label)
        ht_nodes <- object$host_tree$Nnode
        cat("\nHost tree: ")
        cat("  Number of tips:", ht_tips, "\n")
        cat("  Number of nodes:", ht_nodes, "\n")
        if (is.null(object$host_tree$edge.length))
            cat("  No branch lengths.\n")
        else {
            cat("  Branch lengths:\n")
            print(summary(object$host_tree$edge.length)[-4])
        }
        if (is.null(object$host_tree$root.edge))
            cat("  No root edge.\n")
        else
            cat("  Root edge:", object$host_tree$root.edge, "\n")
        if (ht_tips <= 10) {
            cat("  Tip labels:", object$host_tree$tip.label[1], "\n")
            cat(paste("             ",
                 object$host_tree$tip.label[-1]),
                 sep = "\n")
        }
        else {
            cat("  First ten tip labels:", object$host_tree$tip.label[1],
             "\n")
            cat(paste("                       ",
             object$host_tree$tip.label[2:10]),
             sep = "\n")
        }
        if (is.null(object$node.label))
            cat("  No node labels.\n")
        else {
            if (ht_nodes <= 10) {
                cat("  Node labels:",
                    object$host_tree$node.label[1],
                    "\n")
                cat(paste("              ",
                          object$host_tree$node.label[-1]),
                          sep = "\n")
            }
            else {
                cat("  First ten node labels:",
                    object$host_tree$node.label[1],
                    "\n")
                cat(paste("                        ",
                     object$host_tree$node.label[2:10]),
                     sep = "\n")

            }
        }
    }
    # symbiont tree summary
    if (is.null(object$symb_tree)) {
        cat("\n Host tree is not a proper tree.\n")
    }
    else {
        st_tips <- length(object$symb_tree$tip.label)
        st_nodes <- object$symb_tree$Nnode
        cat("\n\nSymb tree: ")
        cat("  Number of tips:", st_tips, "\n")
        cat("  Number of nodes:", st_nodes, "\n")
        if (is.null(object$symb_tree$edge.length))
            cat("  No branch lengths.\n")
        else {
            cat("  Branch lengths:\n")
            print(summary(object$symb_tree$edge.length)[-4])
        }
        if (is.null(object$symb_tree$root.edge))
            cat("  No root edge.\n")
        else
            cat("  Root edge:", object$symb_tree$root.edge, "\n")
        if (st_tips <= 10) {
            cat("  Tip labels:", object$symb_tree$tip.label[1], "\n")
            cat(paste("             ",
                object$symb_tree$tip.label[-1]),
                sep = "\n")
        }
        else {
            cat("  First ten tip labels:",
                object$symb_tree$tip.label[1],
                "\n")
            cat(paste("                       ",
                object$symb_tree$tip.label[2:10]),
                sep = "\n")
        }
        if (is.null(object$node.label))
            cat("  No node labels.\n")
        else {
            if (st_nodes <= 10) {
                cat("  Node labels:",
                    object$symb_tree$node.label[1],
                    "\n")
                cat(paste("              ",
                 object$symb_tree$node.label[-1]),
                  sep = "\n")
            }
            else {
                cat("  First ten node labels:",
                    object$symb_tree$node.label[1],
                    "\n")
                cat(paste("                        ",
                    object$symb_tree$node.label[2:10]),
                    sep = "\n")

            }
        }
    }

    cat("\n\nAssociation Matrix")
    cat("\n There are ",
        nrow(object$association_mat),
        " rows (i.e. extant symbionts.")
    cat("\n There are ",
        ncol(object$association_mat),
         " cols (i.e. extant hosts.")

    if (is.null(object$event_history)) {
        cat("\n No event history.\n")
    }
    else {
        cat("Event history summary: \n\n", summary(object$event_history))
    }
}
#' Print a cophylogenetic set
#' @description Prints a cophylogenetic set or a list of cophylogenetic sets.
#'
#' @param x An object of class `cophy` or class `multiCophy`
#' @param cophy An object of class `cophy`
#' @param ... Further arguments used in generic classes
#'
#' 
#' @return Print returns NULL. host_tree returns NULL, symb_tree returns NULL,
#' association_mat returns the dimensions of the matrix, event_history
#'  returns NULL.
#' @seealso sim_cophylo_bdp, print for the generic, multiCophy, c.cophy
#' @author Wade Dismukes, Ben Bolker, and Emmanuel Paradis
#' @examples
#' h_lambda <- 1.0
#' h_mu <- 0.3
#' c_lambda <- 0.0
#' s_lambda <- 1.0
#' s_mu <- 0.3
#' s_her <- 0.0
#' host_symb_sets <- sim_cophylo_bdp(hbr = h_lambda,
#'                                   hdr = h_mu,
#'                                   sbr = s_lambda,
#'                                   cosp_rate = c_lambda,
#'                                   sdr = s_mu,
#'                                   host_exp_rate = s_her,
#'                                   time_to_sim = 1.0,
#'                                   numbsim = 4)
#' print(host_symb_sets[[1]])
#' host_tree(host_symb_sets[[1]])
#' symb_tree(host_symb_sets[[1]])
#' association_mat(host_symb_sets[[1]])
#' event_history(host_symb_sets[[1]])
#' print(host_symb_sets)
#' @export
print.cophy <- function(x, ...) {
    cat("\n Host Tree:\n\n")
    ape::print.phylo(x$host_tree)
    cat("\n")
    cat("\n Symbiont Tree:\n\n")
    ape::print.phylo(x$symb_tree)
    cat("\n Association Matrix: \n\n")
    cat("\n There are ",
        nrow(x$association_mat),
        " rows (i.e. extant symbionts).")
    cat("\n There are ",
        ncol(x$association_mat),
        " cols (i.e. extant hosts).")
}
#' @describeIn print.cophy Prints a list of cophylogenetic sets
#' @param details A logical value, outputs brief summary of each
#' set in the list.
#' @export
print.multiCophy <- function(x, details = FALSE, ...) {
    num_cophys <- length(x)
    cat(num_cophys,
       "cophylogenetic",
        ifelse(num_cophys > 1, "sets.\n", "set.\n"))
    if (details) {
        for (i in 1:num_cophys) {
            cat("Cophylogenetic set", i, ":\n")
            cat("\tSymbiont tree has ",
                length(x[[i]]$symb_tree$tip.label),
                " tips.\n")
            cat("\tHost tree has ",
                length(x[[i]]$host_tree$tip.label),
                " tips.\n")
        }
    }
}


`$.multiCophy` <- function(x, name) x[[name]]

"[.multiCophy" <- function(x, i) {
    oc <- oldClass(x)
    class(x) <- NULL
    structure(x[i], host_tree = attr(x, "host_tree"),
              class = oc)
}
#' Retrieve the structure of a class multiCophy
#' @description Retrieves the structure of the class multiCophy
#' @param object An object of class multiCophy
#' @param ... Potential further arguments to generic str class
#' @return NULL value
#' @examples
#' h_lambda <- 1.0
#' h_mu <- 0.3
#' c_lambda <- 0.0
#' s_lambda <- 1.0
#' s_mu <- 0.3
#' s_her <- 0.0
#' host_symb_sets <-  sim_cophylo_bdp(hbr = h_lambda,
#'                                   hdr = h_mu,
#'                                   sbr = s_lambda,
#'                                   cosp_rate = c_lambda,
#'                                   sdr = s_mu,
#'                                   host_exp_rate = s_her,
#'                                   time_to_sim = 1.0,
#'                                   numbsim = 2)
#' str(host_symb_sets)
#'
#' @export
str.multiCophy <- function(object, ...)
{
    class(object) <- NULL
    cat('Class "multiCophy"\n')
    str(object, ...)
}


.c_cophy_single <- function(cophy) structure(list(cophy), class = "multiCophy")
#' Combine cophylogenetic sets into a multiCophy object
#' @description Combines cophylogenetic sets into a multiCophy object.
#' @param ... Values of class `cophy`
#' @return An object of class `multiCophy`
#' @seealso `c` generic function
#' @author Wade Dismukes and Emmanuel Paradis
#' @examples
#' h_lambda <- 1.0
#' h_mu <- 0.3
#' c_lambda <- 0.0
#' s_lambda <- 1.0
#' s_mu <- 0.3
#' s_her <- 0.0
#' host_symb_sets <-  sim_cophylo_bdp(hbr = h_lambda,
#'                                   hdr = h_mu,
#'                                   sbr = s_lambda,
#'                                   cosp_rate = c_lambda,
#'                                   sdr = s_mu,
#'                                   host_exp_rate = s_her,
#'                                   time_to_sim = 1.0,
#'                                   numbsim = 2)
#' host_symb_sets2 <- sim_cophylo_bdp(hbr = h_lambda,
#'                                   hdr = h_mu,
#'                                   sbr = s_lambda,
#'                                   cosp_rate = c_lambda,
#'                                   sdr = s_mu,
#'                                   host_exp_rate = s_her,
#'                                   time_to_sim = 1.0,
#'                                   numbsim = 2)
#' multi_host_symb <- c(host_symb_sets[[1]], host_symb_sets2[[2]])
#' multi_host_symb_alt <- c(host_symb_sets, host_symb_sets2)
#' @export
c.cophy <- function(...) {
    recursive <- TRUE
    obj <- list(...)
    classes <- lapply(obj, class)
    iscophylo <- sapply(classes, function(x) "cophy" %in% x)
    if (all(iscophylo)) {
        class(obj) <- "multiCophy"
        return(obj)
    }
    if (!recursive) return(obj)
    ismulti <- sapply(classes, function(x) "multiCophy" %in% x)
    if (all(iscophylo | ismulti)) {
        for (i in which(iscophylo)) obj[[i]] <- .c_cophy_single(obj[[i]])
        ## added by Klaus:
        for (i in which(ismulti)) obj[[i]] <- .uncompressSetTipLabel(obj[[i]])
        obj <- .makeMultiCophyFromObj(obj)
    } else {
        warning('some objects not of class "cophy" or "multiCophy": argument recursive=TRUE ignored')
    }
    obj
}

# modified from function by Klaus Schliep
.makeMultiCophyFromObj <- function(obj)
{
    n <- length(obj)
    N <- lengths(obj, FALSE)
    cs <- c(0, cumsum(N))
    x <- vector("list", cs[length(cs)])
    for (i in 1:n) {
        a <- cs[i] + 1L
        b <- cs[i + 1L]
        x[a:b] <- obj[[i]]
    }
    class(x) <- "multiCophy"
    x
}
#' @describeIn c.cophy Combines two multiCophy objects into one multiCophy object
#' @export
c.multiCophy <- function(...) {
    recursive <- TRUE
    obj <- list(...)
    if (!recursive) return(obj)
    classes <- lapply(obj, class)
    iscophy <- sapply(classes, function(x) "cophy" %in% x)
    ismulti <- sapply(classes, function(x) "multiCophy" %in% x)
    if (!all(iscophy | ismulti)) {
        warning('some objects not of class "cophy" or "multiPhylo": argument recursive=TRUE ignored')
        return(obj)
    }
    for (i in which(iscophy)) obj[[i]] <- .c_cophy_single(obj[[i]])
    ## added by Klaus
    for (i in which(ismulti)) obj[[i]] <- .uncompressSetTipLabel(obj[[i]])
    .makeMultiCophyFromObj(obj)
}

.uncompressSetTipLabel <- function(x) {
    HostLab <- attr(x$host_tree, "TipLabel")
    SymbLab <- attr(x$symb_tree, "TipLabel")
    if (is.null(HostLab)) return(x)
    if (is.null(SymbLab)) return(x)
    class(x) <- NULL
    for (i in seq_len(length(x))) {
        x[[i]]$host_tree$tip.label <- HostLab
        x[[i]]$symb_tree$tip.label <- SymbLab
    }
    class(x) <- "multiCophy"
    attr(x$host_Tree, "TipLabel") <- NULL
    attr(x$symb_tree, "TipLabel") <- NULL
    x
}

`[<-.multiCophy` <- function(x, ..., value) {
    ## recycling is allowed so no need to check: length(value) != length(..1)

    ## check that all elements in 'value' inherit class "phylo"
    test <- unlist(lapply(value, function(xx) !inherits(xx, "cophy")))
    if (any(test))
        stop("at least one element in 'value' is not of class 'cophy'.")

    oc <- oldClass(x)
    class(x) <- NULL

    if (is.null(attr(x, "host_tree"))) {
        x[..1] <- value
        class(x) <- oc
        return(x)
    }

    if (is.null(attr(x, "symb_tree"))) {
        x[..1] <- value
        class(x) <- oc
        return(x)
    }

    x[..1] <- 0L # in case x needs to be elongated
    class(x) <- oc
    j <- 1L
    for (i in ..1) {
        ## x is of class "multiPhylo", so this uses the operator below:
        x[[i]] <- value[[j]]
        j <- j + 1L
    }
    x
}

"[[.multiCophy" <- function(x, i) {
    class(x) <- NULL
    cophy <- x[[i]]
    if (!is.null(attr(x, "host_tree")))
        cophy$host_tree <- attr(x, "host_tree")
    if (!is.null(attr(x, "symb_tree")))
        cophy$symb_tree <- attr(x, "symb_tree")
    if (!is.null(attr(x, "association_mat")))
        cophy$association_mat <- attr(x, "association_mat")
    cophy
}

`[[<-.multiCophy` <- function(x, ..., value) {
    if (!inherits(value, "cophy"))
        stop('trying to assign an object not of class "cophy" into an object of class "multiCophy".')

    oc <- oldClass(x)
    class(x) <- NULL

    ht <- attr(x, "host_tree")

    if (!is.null(ht)) {
        n <- length(ht$tip.label)
        if (n != length(value$host_tree$tip.label))
            stop("tree with different number of tips than those in the list (which all have the same labels; maybe you want to uncompress them)")

        o <- match(value$host_tree$tip.label, ht$tip.label)
        if (any(is.na(o)))
            stop("tree tip labels do not match with those in the list; maybe you want to uncompress them.")
        value$host_tree$tip.label <- NULL
        ie <- match(o, value$host_tree$edge[, 2])
        value$host_tree$edge[ie, 2] <- 1:n
    }

    st <- attr(x, "host_tree")

    if (!is.null(st)) {
        n <- length(st$tip.label)
        if (n != length(value$symb_tree$tip.label))
            stop("Host tree with different number of tips than those in the list (which all have the same labels; maybe you want to uncompress them)")

        o <- match(value$symb_tree$tip.label, st$tip.label)
        if (any(is.na(o)))
            stop("tree tip labels do not match with those in the list; maybe you want to uncompress them.")
        value$symb_tree$tip.label <- NULL
        ie <- match(o, value$symb_tree$edge[, 2])
        value$symb_tree$edge[ie, 2] <- 1:n
    }
    x[[..1]] <- value
    class(x) <- oc
    x
}

`$<-.multiCophy` <- function(x, ..., value) {
    x[[..1]] <- value
    x
}