# to get a list of element for sparsebnFit object from the adjacency matrix.
# "graph" is the adjacency matrix
get.summary <- function(graph, dataSize, lambda, time) {
  lambda = lambda
  nedge <- sum(graph)
  pp <- ncol(graph)
  nn <- dataSize
  graph_list <- lapply(seq_len(ncol(graph)), function(i) graph[,i])
  edges <- lapply(graph_list, function(x){
    return(which(x==1))
  })
  edges <- sparsebnUtils::edgeList(edges)
  return(list(edges = edges, lambda = lambda, nedge = nedge, pp = pp, nn = nn, time = time))
}

# to convert an output solution path to a list of nlam elements
# each element is a list contains all required items for sparsebnFit
get.edgeList <- function(edgeMatrix, dataSize, lambda, time) {
  n_node <- ncol(edgeMatrix)
  n_graph <- nrow(edgeMatrix)/n_node
  edgeList <- vector("list", n_graph)
  for (i in 1:n_graph) {
    graph <- edgeMatrix[((i-1)*n_node+1):(i*n_node), ]
    edgeList[[i]] <- get.summary(graph, dataSize, lambda[i], time[i])
  }
  return(edgeList)
}

# to get obsIndex from a given intervention list (ivn).
# if for a node, it is under intervention for all observations, return 0.
get_obsIndex <- function(ivn, node) {
  obsIndex <- vector("list", length = node)
  ind <- 1:node
  obsIndex <- lapply(ind, function(x, ivn) {
    if_in <- sapply(ivn, function(y) {x %in% y})
    observation <- 1:length(ivn)
    obsIndex_one <-observation[!if_in]
    if (length(obsIndex_one)==0) {obsIndex_one = 0L}
    obsIndex_one
  }, ivn)
  return(obsIndex)
}

# to get adaptive weights if given L2 norm of parameters
get_adaptWeights <- function(beta_l2, max.weights = 10^5) {
  beta_l2 <- beta_l2 + t(beta_l2)
  cd.weights <- 1/beta_l2
  cd.weights[which(cd.weights==Inf)] = max.weights

  return(cd.weights)
}

dat_transform <- function(datbn) {
  data_matrix <- datbn$data
  for (i in 1:ncol(data_matrix)) {
    if (!is.numeric(data_matrix[[i]]) && !is.factor(data_matrix[[i]])) stop("data must be a data.frame object with integer or factor columns!")
  }
  data_matrix <- as.data.frame(sapply(data_matrix, function(x){as.integer(x)}))

  # test if the number of levels in data set is compatible with the number of true levels.
  dat_levels <- sapply(datbn$levels, length)
  max_levels <- sapply(data_matrix, max)
  if (sum(max_levels>dat_levels)) {
    data_matrix <- as.data.frame(sapply(data_matrix, function(x){as.integer(as.factor(x))}))
  }
  max_levels <- sapply(data_matrix, max)
  if (sum(max_levels>dat_levels)) stop("The number of levels and the data set is not compatible! Check data set and the input ivn list!")
  # adjust the levels to start from 0
  min_levels <- sapply(data_matrix, min)
  data_matrix <- as.data.frame(sapply(1:ncol(data_matrix), function(x, data_matrix, min_levels){as.integer(data_matrix[, x]-min_levels[x])}, data_matrix, min_levels))
  return(data_matrix)
}

get_bwlist <- function(bwlist, node_index){
  int_bwlist <- match(bwlist, node_index)
  int_bwlist <- matrix(int_bwlist, ncol = 2)
  return(int_bwlist)
}
