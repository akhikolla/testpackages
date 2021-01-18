#' generate_discrete_data
#'
#' data generating function
#'
#' @param graph A \code{\link[sparsebnUtils]{edgeList}} object.
#' @param params Coefficient list.
#' @param n Size of the data set, a scalar
#' @param ivn List of interventions.
#' @param ivn.rand If \code{TRUE}, random values will be drawn uniformly for each intervention. Otherwise, these values need to supplied manually in \code{ivn}.
#' @param n_levels A vector of number of levels for each node. Default is binary data.
#'
#' @examples
#'
#' ### generate observational data
#' gr <- sparsebnUtils::random.graph(5, 5) # use sparsebnUtils package to generate a random graph
#' names(gr) = c("V1", "V2", "V3", "V4", "V5")
#' nlevels <- c(3, 5, 2, 2, 3)
#' gr.params <- coef_gen(edge_list = gr, n_levels = nlevels)
#' data.obs <- discretecdAlgorithm::generate_discrete_data(graph = gr,
#'                                                         n = 100,
#'                                                         n_levels = nlevels,
#'                                                         params = gr.params)
#'
#' ### generate experimental data
#' ivn <- as.list(c(rep("V1", 50), rep("V2", 50))) # 50 interventions on V1, 50 interventions on V2
#' data.ivn <- discretecdAlgorithm::generate_discrete_data(graph = gr,
#'                                              n = 100,
#'                                              n_levels = nlevels,
#'                                              params = gr.params,
#'                                              ivn = ivn)
#'
#' ###  Use pre-specified values for interventions
#' ###  In this toy example, we assume that all intervened nodes were fixed to
#' ###  to the value 1, although this can be any number of course.
#' ivn.vals <- lapply(ivn, function(x) sapply(x, function(x) 1)) # replace all entries with a 1
#' data.ivn <- discretecdAlgorithm::generate_discrete_data(graph = gr,
#'                                              n = 100,
#'                                              n_levels = nlevels,
#'                                              params = gr.params,
#'                                              ivn = ivn.vals,
#'                                              ivn.rand = FALSE)
#'
#' @return data matrix
#' @export
generate_discrete_data <- function(graph,
                                   params,
                                   n,
                                   ivn = NULL,
                                   ivn.rand = TRUE,
                                   n_levels = NULL)
{
  datGen_call(edge_list = graph,
              dataSize = n,
              ivn = ivn,
              ivn_rand = ivn.rand,
              nlevels = n_levels,
              coef = params)
}

datGen_call <- function(edge_list,
                        dataSize,
                        ivn,
                        ivn_rand,
                        nlevels,
                        coef)
{
  # check input
  if(!sparsebnUtils::is.edgeList(edge_list)) stop("edge_list must be a edgeList object!")
  if(is.null(names(edge_list))) stop("edge_list must be named!")

  node_name <- names(edge_list)
  adj_matrix <- sparsebnUtils::get.adjacency.matrix(edge_list)
  dag_igraph <- igraph::graph.adjacency(as.matrix(adj_matrix))
  # dag_igraph <- sparsebnUtils::to_igraph(edge_list)

  edge_list <- as.list(edge_list)
  edge_list <- lapply(edge_list, as.integer)

  maxdeg <- max(sapply(edge_list, length))
  maxdeg <- as.integer(maxdeg)

  node <- length(edge_list)
  node<- as.integer(node)

  ts <- match(node_name, names(igraph::topo_sort(dag_igraph))) # for master brand
  # ts <- as.integer(as.vector(igraph::topo_sort(dag_igraph))) # for dev brand

  if (is.null(ts)) stop("Need topological sort for the graph!")
  if (length(ts) != node) stop("length of ts should be the same with node!")

  ordex <- sapply(edge_list, function(x, maxdeg){
    as.integer(c(x, rep(0, maxdeg-length(x))))
  }, maxdeg)
  ordex <- matrix(ordex, nrow = maxdeg)

  if(!is.numeric(dataSize) || length(dataSize) > 1) stop("data_size must be a scalar!")
  if(dataSize < 1) stop("data_size must be a positive integer!")
  dataSize <- as.integer(dataSize)

  if (is.null(nlevels)) {
    nlevels <- rep(2, node)
  }
  if(!is.vector(nlevels)) stop("n_levels must be a vector!")
  if(sum(nlevels<2)) stop("number of levels must be at least 2!")
  if(length(nlevels)!=node) stop("length of n_levels not compatible with edge_list!")
  nlevels <- as.integer(nlevels)

  ivn_vals <- as.list(rep(0, dataSize))

  if (is.null(ivn)) {
    ivn <- as.list(rep(0, dataSize))
  }
  else {
    if(!is.list(ivn)) stop("ivn must be a list!")
    if(ivn_rand) {
      if (sum(sapply(ivn, function(x){!is.character(x)})) && sum(sapply(ivn, function(x){!is.integer(x)}))) stop("ivn must be a list of all characters or all integers, or NA")
      if (sum(sapply(ivn, function(x){is.character(x)}))) {
        ivn <- lapply(ivn, function(x) match(x, node_name))
      }
    }
    else {
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
      ivn_vals <- lapply(ivn, as.integer)
      ivn <- lapply(ivn, function(x){names(x)})
      ivn <- lapply(ivn, function(x) match(x, node_name))
    }
  }
  if(length(ivn)!=dataSize) stop("length of ivn not compatible with data_size")
  ivn <- lapply(ivn, function(x){
    if (is.null(x)) {
      x <- 0
    }
    as.integer(x-1)})

  # if(is.null(coef)) stop("coef must have some value!")
  if(!is.list(coef)) stop("coef must be a list!")
  if(length(coef)!=node) stop("length of coef must be the same with number of nodes!")
  # check type of list element
  if (sum(sapply(coef, function(x){!is.matrix(x) && !is.null(x)}))) stop("element of coeff must be matrix!")
  # check dimension of the list element
  flag=FALSE
  for (i in 1:node) {
    if (length(edge_list[[i]])) {
      if (nrow(coef[[i]]) != nlevels[i])
        flag = TRUE
      if (ncol(coef[[i]]) != sum(nlevels[edge_list[[i]]]-1)+1) {
        flag = TRUE
      }
    } else {
      if (!is.null(coef[[i]])){
        flag = TRUE
      }
    }
  }
  if (flag == TRUE) stop("coef does not compatible with edge_list!")
  coef_list <- lapply(1:node, function(x, coef, edge_list, nlevels){
    out <- NULL
    if (length(edge_list[[x]])) {
      index <- rep(1:(length(edge_list[[x]])+1), c(1, nlevels[edge_list[[x]]]-1))
      m_to_list <- vector("list", length = length(edge_list[[x]])+1)
      for (i in 1:(length(edge_list[[x]])+1)) {
        m_to_list[[i]] <- coef[[x]][, index==i, drop = FALSE]
      }
      out = m_to_list
    }
    out
  }, coef, edge_list, nlevels)

  # coef_length
  coef_length <- sapply(coef_list, length)
  coef_length <- as.integer(coef_length)

  # call DatGen_cpp
  DatGen_cpp(maxdeg, node, ordex, ts, dataSize, ivn, ivn_vals, ivn_rand, nlevels, coef_list, coef_length)
}


# function that directly calls from cpp
# type check, no converting type of an input at this point

DatGen_cpp <- function(maxdeg,
                   node,
                   ordex,
                   ts,
                   dataSize,
                   ivn,
                   ivn_vals,
                   ivn_rand,
                   nlevels,
                   coef_list,
                   coef_length)
{
  # check for maxdeg
  if(!is.integer(maxdeg)) stop("maxdeg must be an integer!")
  if(maxdeg <= 0) stop("maxdeg must be a positive integer!")

  # check for node
  if(!is.integer(node)) stop("node must be an integer!")
  if(node <= 0) stop("node must be a positive integer!")

  # check for ordex
  if(!is.matrix(ordex)) stop("ordex must be a matrix!")
  if(sum(sapply(ordex, function(x){!is.integer(x)}))) stop ("ordex has to be a matrix with integer entries!")
  if(nrow(ordex) != maxdeg) stop("Incompatible size! Number of rows of ordex should maxdeg!")
  if(ncol(ordex) != node) stop("Incompatible size! Number of columns of ordex should be node!")

  # check for ts
  if(!is.vector(ts)) stop("ts must be a vector!")
  if(length(ts)!=node) stop("ts must have length node!")
  if(sum(!is.integer(ts))) stop("ts must be a vector of integer elements!")
  # if(sum(sort(ts, decreasing = FALSE)!=(1:node))) stop("ts is not a valid topological sort for node!")

  # check for dataSize
  if(!is.integer(dataSize)) stop("dataSize must be an integer!")
  if(dataSize <= 0) stop("dataSize must be a positive integer!")

  # check for ivn
  if(!is.list(ivn)) stop("ivn must be a list")
  if(sum(sapply(ivn, function(x){sum(!is.integer(x))}))) stop("ivn must be a list of integer vectors!")
  if(length(ivn)!=dataSize) stop("ivn must be a list of length dataSize!")

  # check for nlevels
  if(!is.vector(nlevels)) stop("nlevels must be a vector!")
  if(sum(!is.integer(nlevels))) stop("nlevels must be a vector of integers!")
  if(length(nlevels)!=node) stop("length of nlevels must be node!")

  # check for coef
  if(!is.list(coef_list)) stop("coef must be a list!")
  if(sum(sapply(coef_list, function(x){!is.list(x) && !is.null(x)}))) stop("coef_list must be a list of list")
  for (i in 1:node) {
    if (length(coef_list[[i]])) {
      for (j in 1:length(coef_list[[i]])) {
        if (!is.matrix(coef_list[[i]][[j]])) {
          # cat("error!")
          stop("element of coef_list must be NULL or a matrix!")
        }
      }
    }
  }

  if (length(coef_length)!= node) stop("length of coef_length should be node!")

  for (i in 1:node) {
    if (length(coef_list[[i]])) {
      temp_parent <- ordex[, i][which(ordex[, i]!=0)]
      if (length(coef_list[[i]])!=length(temp_parent)+1) stop("Length of coef_list sublist not compatible with number of parents!")
      if (nrow(coef_list[[i]][[1]])!=nlevels[i]) stop("number of rows for intercept should be the numebr of levels of the node!")
      if (ncol(coef_list[[i]][[1]])!=1) stop("number of columns for itercept should be 1!")
      for (j in 2:length(coef_list[[i]])) {
        if (nrow(coef_list[[i]][[j]])!=nlevels[i]) stop("number of rows for coefficient matrix for an edge should be the numebr of levels of the child!")
        if (ncol(coef_list[[i]][[j]])!=nlevels[temp_parent[j-1]]-1) stop("number of columns for coeffecient matrix for an edge should be the numebr of levels of parent minus 1!")
      }
    }
  }

  if (!is.integer(coef_length)) stop("coef_length must be an integer!")

  # call function from cpp
  DatGen(maxdeg, node, ordex, ts, dataSize, ivn, ivn_vals, ivn_rand, coef_length, nlevels, coef_list)
}

#' coef_gen
#'
#' coefficient generating function
#'
#' @param edge_list a \code{\link[sparsebnUtils]{edgeList}} object.
#' @param n_levels, a list of number of levels for each node.
#' @param FUN, a probability distribution to generate coefficients
#' @param flip, a bool parameter. If true, will randomly flip the sign of coefficients.
#' @return A list of coefficient matrix
#' @export
coef_gen <- function(edge_list, n_levels = NULL, FUN=NULL, flip=TRUE) {
  if (!sparsebnUtils::is.edgeList(edge_list)) stop("edge_list must be an edge_list object!")

  if (is.null(FUN)) {
    FUN <- function(n) {
      stats::runif(n, 1, 3)
    }
  }
  if(is.null(n_levels)) {
    n_levels <- rep(2, length(edge_list))
  }
  if(sum(n_levels <= 1)) stop("Number of levels for each node must be at least 2!")

  coef <- lapply(1:length(edge_list), function(x, edge_list, flip){
    if (length(edge_list[[x]])==0) {
      coef_matrix <- NULL;
    }
    else {
      ncol_coef <- 1
      for (i in 1:length(edge_list[[x]])) {
        ncol_coef = ncol_coef + n_levels[edge_list[[x]][i]]-1
      }
      n_coef <- n_levels[x]*ncol_coef
      coef_vector <- replicate(n_coef, FUN(n=1))
      if (flip) {
        flip_flag <- sample(c(-1, 1), size = n_coef, replace = TRUE)
        coef_vector <- coef_vector * flip_flag
      }
      coef_matrix <- matrix(coef_vector, nrow = n_levels[x])
    }
    return(coef_matrix)
  }, edge_list, flip)
  return(coef)
}

#' data_gen
#'
#' A function that generate discrete data set.
#'
#' @param graph a \code{\link[sparsebnUtils]{edgeList}} object.
#' @param n size of the data set, a scalar
#' @param ivn, a list of intervention for each data point.
#' @param n_levels, a list of number of levels for each node, default is binary data set.
#' @param params, coefficient list (optional).
#' @param FUN, a function to generate magnitude of influence (optional).
#' @param flip, a bool parameter. If true, when generating coefficients, will randomly flip the sign of coefficients.
#' @return data matrix
#' @export
data_gen <- function(graph,
                     n,
                     ivn = NULL,
                     n_levels = NULL,
                     params = NULL,
                     FUN = NULL,
                     flip = TRUE)
{
  # check graph
  if(!sparsebnUtils::is.edgeList(graph)) stop("graph must be an edgeList object!")

  # check n
  if (n<1) stop("n must be a postive number!")

  # check n_levels
  if(is.null(n_levels)) {
    n_levels <- rep(2, length(graph))
  }
  if(sum(n_levels <= 1)) stop("Number of levels for each node must be at least 2!")

  # check coef
  if(is.null(params)) {
    params <- coef_gen(graph, n_levels, FUN, flip)
  }

  # call data_gen
  generate_discrete_data(graph = graph, params = params, n = n, ivn = ivn, n_levels = n_levels)
}
