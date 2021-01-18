# discretize.jointly.R
#
# Created by Sajal Kumar
# Modified by Jiandong Wang
# Copyright (c) NMSU Song lab

#' Discretize Multivariate Continuous Data by a Cluster-Preserving Grid
#'
#' Discretize multivariate continuous data using a grid that captures the joint distribution via
#' preserving clusters in the original data
#'
#' @importFrom cluster silhouette
#' @importFrom stats kmeans
#' @importFrom fossil adj.rand.index
#' @importFrom stats dist
#' @importFrom dqrng dqsample
#' @importFrom Rdpack reprompt
#' @import Rcpp
#' @useDynLib GridOnClusters
#'
#' @param data a matrix containing two or more continuous variables.
#' Columns are variables, rows are observations.
#'
#' @param k either the number or range of clusters to be found on \code{data}.
#' The default is 2 to 10 clusters. If a range is specified, an optimal k in
#' the range is chosen to maximize the average silhouette width.
#' If \code{cluster_label} is specified, \code{k} is ignored.
#'
#' @param cluster_label a vector of user-specified cluster labels for each observation
#' in \code{data}. The user is free to choose any clustering.
#' If unspecified, k-means clustering is used by default.
#'
#' @param min_level the minimum number of levels along each dimension 
#'
#' @details 
#' 
#' The function implements algorithms described in \insertCite{Jwang2020BCB}{GridOnClusters}.
#' 
#'
#' @return
#'
#' A list that contains four items:
#' \item{\code{D}}{a matrix that contains the discretized version of the original \code{data}.
#' Discretized values are one(1)-based.}
#'
#' \item{\code{grid}}{a list of vectors containing decision boundaries for each variable/dimension.}
#'
#' \item{\code{clabels}}{a vector containing cluster labels for each observation in \code{data}.}
#'
#' \item{\code{csimilarity}}{a similarity score between clusters from joint discretization
#' \code{D} and cluster labels \code{clabels}. The score is the adjusted Rand index.}
#' 
#' @examples
#' # using a specified k
#' x = rnorm(100)
#' y = sin(x)
#' z = cos(x)
#' data = cbind(x, y, z)
#' discretized_data = discretize.jointly(data, k=5)$D
#'
#' # using a range of k
#' x = rnorm(1000)
#' y = log1p(abs(x))
#' z = tan(x)
#' data = cbind(x, y, z)
#' discretized_data = discretize.jointly(data, k=c(3:10))$D
#'
#' # using an alternate clustering method to k-means
#' library(cluster)
#' x = rnorm(1000)
#' y = log1p(abs(x))
#' z = sin(x)
#' data = cbind(x, y, z)
#'
#' # pre-cluster the data using partition around medoids (PAM)
#' cluster_label = pam(x=data, diss = FALSE, metric = "euclidean", k = 5)$clustering
#' discretized_data = discretize.jointly(data, cluster_label = cluster_label)$D
#'
#' @references
#' \insertAllCited{}
#' 
#' @author 
#' 
#' Jiandong Wang, Sajal Kumar and Mingzhou Song
#' 
#' @seealso
#'
#' See \link[Ckmeans.1d.dp]{Ckmeans.1d.dp} for discretizing univariate continuous data.
#'
#' @export
discretize.jointly = function(data, k=c(2:10), cluster_label=NULL, min_level = 2){

  # check if data provided is a matrix
  if( !("matrix" %in% class(data)) && !("data.frame" %in% class(data))){
    stop("'data' must be a matrix or data.frame.")
  }

  # check if all columns in 'data' are numeric or integer
  dim_class = apply(data, 2, class)
  if(!all(dim_class == "numeric" | dim_class == "integer")){
    stop("All columns in 'data' should be numeric or integer.")
  }

  # 'data' should have atleast 10 points
  dim_data = dim(data)
  if(dim_data[1] < 10){
    stop("'data' should have atleast 10 observations,")
  }

  # 'cluster_label' should either be null or integers (or numeric) matching nrow(data)
  if(!is.null(cluster_label) && length(cluster_label) != nrow(data)){
    stop("'cluster_label' should either be null or a vector with nrow(data) elements.")
  }

  if(!is.null(cluster_label) && !class(cluster_label) %in% c("numeric","integer")){
    stop("'cluster_label' should be either null or a numeric/integer vector.")
  }
  
  # 'min_level' should smaller then the max of k
  if(max(k)<min_level){
    k = c(min_level,min_level+5)
    warning("'min_level' should be in the range of k, k has been adapted")
  }

  # if no cluster labels are supplied, default to K-means
  if(is.null(cluster_label)){

    # is k a single number or a range
    if(length(k) == 1){

      # randomly generate centers
      centers = dqsample(which(!duplicated(data)), k)

      # get cluster information for 'k'
      cluster_info = kmeans(data, centers = as.matrix(data[centers,]))

      # only keep cluster centers and labels
      cluster_info = list(#centers = as.matrix(cluster_info$centers),
                          clusters = as.vector(cluster_info$cluster)-1,
                          data = as.matrix(data))

    } else {

      data_dist = dist(data) # distance matrix for data

      # compute cluster info and silhouette score the cluster range k
      cluster_info = lapply(k, function(i){
        centers = dqsample(which(!duplicated(data)), i)
        data_clust = kmeans(data, centers = as.matrix(data[centers,]))
        return(list(data_clust$cluster, data_clust$centers, mean(silhouette(data_clust$cluster, dist=data_dist)[,3])))
      })

      # find the max silhouette
      silhouette_scr = unlist(lapply(cluster_info, function(x){
        return(x[[3]])
      }), use.names = FALSE)
      max_silhouette = max(which.max(silhouette_scr)) # take the bigger k out of those with equal silhouette scores

      # only keep cluster centers and labels for the 'k' with max silhouette
      cluster_info = list(#centers = as.matrix(cluster_info[[max_silhouette]][[2]]),
                          clusters = as.vector(cluster_info[[max_silhouette]][[1]]-1),
                          data = as.matrix(data))
    }
  } else { # cluster labels are supplied

    cluster_label = as.numeric(as.factor(cluster_label)) # making labels consecutive

    # # find medoids
    # centers = data[medoids(D = dist(data), cl = cluster_label),,drop=FALSE]

    cluster_info = list(#centers = as.matrix(centers),
                        clusters = cluster_label-1,
                        data = as.matrix(data))
  }

  # get grid lines
  grid_lines = findgrid(cluster_info, length(unique(cluster_info$clusters)), nrow(data), ncol(data), min_level)

  # filter grid lines
  grid_lines = lapply(grid_lines, function(i){
    mx = max(i)
    return(i[-which(i == mx)])
  })

  # discretize data
  discr_data = discretize_data(data, grid_lines)

  # compute adjusted random index
  ndim_cluster_dist = discr_data[,1]
  for(i in 2:ncol(discr_data)){
    ndim_cluster_dist = paste0(ndim_cluster_dist,",",discr_data[,2])
  }
  cluster_similarity = adj.rand.index(as.numeric(as.factor(ndim_cluster_dist))-1, cluster_info$clusters)

  return(list(D=discr_data, grid=grid_lines, clabels=cluster_info$clusters+1, csimilarity=cluster_similarity))
}

# for internal use
# takes n-dimesional 'data' and gridlines to quantify each dimension
discretize_data = function(data, gridlines){
  
  # use gridlines for each dimension
  for(i in 1:ncol(data)){
    
    if(length(unique(gridlines))==0){
      discr = rep(1, nrow(data))
    }else{
      discr = rep(length(gridlines[[i]])+1, nrow(data))
      gridlines[[i]] = sort(gridlines[[i]])
      for(j in 1:length(gridlines[[i]])){ # determine discretization levels
        if(j == 1) {
          discr[data[,i] < gridlines[[i]][j]] = 1
        } else {
          discr[data[,i] < gridlines[[i]][j] & data[,i] >= gridlines[[i]][j-1]] = j
        }
      }
    }
    data[,i] = discr
  }
  
  return(data)
}