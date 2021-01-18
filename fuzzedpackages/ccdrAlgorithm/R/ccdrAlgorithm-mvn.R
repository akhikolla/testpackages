#
#  ccdrAlgorithm-mvn.R
#  ccdrAlgorithm
#
#  Created by Bryon Aragam (local) on 1/15/17.
#  Copyright (c) 2014-2017 Bryon Aragam. All rights reserved.
#

#' Generate data from a DAG
#'
#' Given a Gaussian DAG, generate data from the underlying distribution.
#' Equivalently, generate data from a multivariate normal distribution given
#' one of its SEM. Can generate both observational and intervention data.
#'
#' If \code{ivn = NULL}, then \code{n} observational samples are drawn. For each
#' component of \code{ivn} that is not \code{NULL}, interventional samples will
#' be drawn with the values of each node specified in the component.
#'
#' @param graph DAG in \code{\link{edgeList}} format.
#' @param params Vector of parameters. Last p elements correspond to variances (p = number of nodes in \code{graph}), initial elements correspond to edge weights.
#' @param n Number of samples to draw.
#' @param ivn List of interventions (see \code{\link[sparsebnUtils]{sparsebnData}}). Must be a \code{list} with exactly \code{n} components.
#' @param ivn.rand If \code{TRUE}, random N(0,1) values will be drawn for each intervention. Otherwise, these values need to supplied manually in \code{ivn}.
#'
#' @examples
#'
#' ### Generate observational data
#' gr <- sparsebnUtils::random.graph(5, 5) # use sparsebnUtils package to generate a random graph
#' gr.params <- runif(10) # there are 5 coefficients + 5 variances
#' data.obs <- ccdrAlgorithm::generate_mvn_data(graph = gr,
#'                                              n = 100,
#'                                              params = gr.params)
#'
#' ### Generate experimental data
#' ivn <- as.list(c(rep("V1", 50), rep("V2", 50))) # 50 interventions on V1, 50 interventions on V2
#' data.ivn <- ccdrAlgorithm::generate_mvn_data(graph = gr,
#'                                              n = 100,
#'                                              params = gr.params,
#'                                              ivn = ivn)
#'
#' ### Use pre-specified values for interventions
#' ###  In this toy example, we assume that all intervened nodes were fixed to
#' ###  to the value 1, although this can be any number of course.
#' ivn.vals <- lapply(ivn, function(x) sapply(x, function(x) 1)) # replace all entries with a 1
#' data.ivn <- ccdrAlgorithm::generate_mvn_data(graph = gr,
#'                                              n = 100,
#'                                              params = gr.params,
#'                                              ivn = ivn.vals,
#'                                              ivn.rand = FALSE)
#'
#' ### If ivn.rand = FALSE, you must specify values
#' ###  The code below will fail because ivn does not contain any values
#' ### (compare to ivn.vals above).
#' \dontrun{
#' data.ivn <- ccdrAlgorithm::generate_mvn_data(graph = gr,
#'                                              n = 100,
#'                                              params = gr.params,
#'                                              ivn = ivn,
#'                                              ivn.rand = FALSE)
#' }
#'
#' @export
generate_mvn_data <- function(graph, params, n = 1, ivn = NULL, ivn.rand = TRUE){
    ### This function requires the 'igraph' package to be installed
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("The igraph package is required for the method 'generate_mvn_data'. Please install it using install.packages(\"igraph\").", call. = FALSE)
    }

    stopifnot(sparsebnUtils::is.edgeList(graph))
    stopifnot(is.numeric(params))
    stopifnot(length(params) == sparsebnUtils::num.edges(graph) + sparsebnUtils::num.nodes(graph))

    if(is.null(names(graph))){
        stop("Input 'graph' requires node names!")
    }

    if(!is.null(ivn)){
        stopifnot(is.list(ivn))
        stopifnot(length(ivn) == n)

        ### Generate random intervention values
        if(ivn.rand){
            ivn <- lapply(ivn, function(x) sapply(x, function(x) rnorm(n = 1, mean = 0, sd = 1))) # assume standard normal
            # ivn <- lapply(ivn, function(x) sapply(x, function(x) 1)) # debugging
        } else{
            check_vals <- sparsebnUtils::check_list_class(ivn, c("NULL", "numeric")) # check to make sure list components are either numeric (ivn vals) or NULL (obs sample)
            check_names <- sapply(ivn, function(x) is.null(names(x))) # return TRUE if component has no names attribute (i.e. it is NULL)

            if(!check_vals || all(check_names)){
                err_msg <- paste0("ivn.rand set to FALSE with invalid input for ivn: ",
                                  "If ivn.rand = FALSE, you must pass explicit values ",
                                  "for each intervention used in your experiments. ",
                                  "Please check that the ivn argument is a list whose ",
                                  "arguments are named numeric vectors whose names ",
                                  "correspond to the node under intervention or NULL ",
                                  "if the corresponding row is observational.")
                stop(err_msg)
            }
        }
    }

    ### Need this to ensure the output has the same order as the input
    ###  after things get shuffled around
    original_node_order <- names(graph)

    ### Get topological sort
    ### Note that the check for the igraph pkg occurs in sparsebnUtils::to_igraph
    topsort <- names(igraph::topo_sort(sparsebnUtils::to_igraph(graph)))

    nnode <- length(original_node_order)
    vars <- utils::tail(params, nnode) # parameters associated with variances
    names(vars) <- original_node_order
    coefs <- params[1:(length(params) - nnode)] # parameters associated with edge weights
    sp <- sparsebnUtils::as.sparse(graph)
    sp$vals <- coefs # previous line leaves NAs for values in sparse object; need to fill these in
    edgelist <- sparse_to_edgeWeightList(sp, original_node_order)
    nodes <- names(edgelist) # this will be sorted according to the topological order

    ### The old way, efficient for obs data only
    # x <- replicate(n, generate_mvn_vector(edgelist, nodes, topsort, vars))
    # x <- t(x)[, original_node_order]

    x <- vector("list", length = n)
    for(i in 1:n){
        x[[i]] <- generate_mvn_vector(edgelist, nodes, topsort, vars, ivn = ivn[[i]])
    }
    x <- do.call("rbind", x)

    ### Permute columns back to original ordering
    x <- x[, original_node_order]
    x
}

generate_mvn_vector <- function(edgelist, nodes, topsort, vars = NULL, ivn = NULL){
    normal_seed <- sapply(vars, function(x) rnorm(n = 1, mean = 0, sd = sqrt(x)))
    gen_dag_vector_R(edgelist, nodes, topsort, seed = normal_seed, ivn = ivn)
}

#
# edgelist = graph information
# nodes = names of nodes in graph
# topsort = topological sort (indexed by node names)
# seed = random noise (Gaussian); bias term (binary)
# ivn = named vector of intervention values (do(child = x))
#
gen_dag_vector_R <- function(edgelist, nodes, topsort, seed, ivn = NULL){
    nnode <- length(edgelist)
    x <- numeric(nnode)
    names(x) <- nodes
    ivnnames <- names(ivn)

    for(j in seq_along(topsort)){
        child <- topsort[j]

        if(child %in% ivnnames){
            ### If node is intervened on, fix value according to input in 'ivn'
            x[child] <- ivn[child]
        } else{
            ### If no intervention, use DAG to determine value from parents
            parents <- edgelist[[child]]$parents
            weights <- edgelist[[child]]$weights
            nparents <- length(parents)
            if(nparents > 0){
                ### Iterate over parents and add associated effects
                for(i in seq_along(parents)){
                    this.par <- parents[i]
                    x[child] <- x[child] + weights[i] * x[this.par]
                    # x[child] <- x[child] + weights[i] * x[index[i]] # equivalent to above line
                }
            }

            ### Add noise: This is a crucial step. If nothing is added here, the
            ###            output will be all zeroes since the root node(s) will
            ###            have x[child] = 0 at this point.
            ###
            ### Gaussian model: This is random error ~ N(0, vars[j])
            ### Logistic model: This a (deterministic) bias term
            x[child] <- x[child] + seed[child]
        }

    }

    x
}

sparse_to_edgeWeightList <- function(x, nodes){
    stopifnot(sparsebnUtils::is.sparse((x)))
    # sp <- sparsebnUtils::as.sparse(x) # NOTE: no longer a bottleneck under sparsebnUtils v0.0.4

    # nodes <- colnames(x)
    stopifnot(x$dim[1] == x$dim[2])

    out <- lapply(vector("list", length = x$dim[1]), function(z) list(parents = character(0), index = integer(0), weights = numeric(0)))
    names(out) <- nodes
    for(j in seq_along(x$cols)){
        child <- x$cols[[j]]
        parent <- x$rows[[j]]
        weight <- x$vals[[j]]
        parents <- c(out[[child]]$parents, nodes[parent]) # !!! THIS IS SLOW
        index <- c(out[[child]]$index, parent) # !!! THIS IS SLOW
        weights <- c(out[[child]]$weights, weight) # !!! THIS IS SLOW
        out[[nodes[child]]] <- list(parents = parents, index = index, weights = weights)
    }

    out
}
