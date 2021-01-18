#' @useDynLib kmcudaR
#' @importFrom Rcpp evalCpp


#' @title K-Means Clustering using CUDA
#' @description Performs k-means clustering on a
#' numeric matrix using a NVIDIA GPU via CUDA
#' @param samples A numeric matrix
#' @param clusters the number of clusters
#' @param tolerance if the relative number of reassignments
#' drops below this value the algorithm stops
#' @param init A character vector or numeric matrix, sets the 
#' method for centroids initialization.  Options include 
#' "k-means++", "afk-mc2", "random" or numeric matrix of
#' shape [clusters, number of features]. Default = "kmeans++"
#' @param yinyang_t numeric value defining relative number of
#' cluster groups.  Usually 0.1 but 0 disables Yinyang refinement.
#' @param metric Character vector specifying distance metric to
#' use.  The default is Euclidean (L2), it can be changed to 
#' "cos" for Sphereical K-means with angular distance. NOTE - 
#' the samples must be normalized in the latter case.
#' @param average_distance logical indicating whether to calculate
#' the average distance between cluster elements and the
#' corresponding centroids.  Useful for finding the best 'K'.
#' Returned as third list element
#' @param seed random generator seed for reproducible results [deprecated]
#' @param device integer defining device to use. 1 = first device,
#' 2 = second device, 3 = first & second devices, 0 = use all devices.
#' Default = 0
#' @param verbosity Integer indicating amount of output to see.
#' 0 = silence, 1 = progress logging, 2 = all output
#' @return a list consisting of
#' \item{centroids}{Cluster centroids}
#' \item{assignments}{integer vector of sample-cluster associations}
#' \item{average_distance}{average distance between cluster elements}
#' @export
kmeans_cuda <- function(
	samples, 	clusters, tolerance = 0.01, 
	init = "k-means++", yinyang_t = 0.1, metric = "L2",
	average_distance = FALSE, seed = NULL, 
	device = 0L, verbosity = 0L)
{
	
	if(length(clusters) > 1){
		stop("clusters must be of length 1")
	}
	
	if(clusters <= 0 || !is.numeric(clusters)){
		stop("clusters must be a positive integer")
	}
	
	if(!is.numeric(tolerance) || tolerance < 0 || tolerance > 1){
		stop("tolerance must be a numeric value between [0, 1]")
	}
	if (yinyang_t < 0 || yinyang_t > 0.5) {
		stop("tolerance must be in [0, 0.5]");
	}
	
	if(is.list(samples)){
		samples <- do.call('rbind', samples)
	}
	
	result <- r_kmeans_cuda(samples, clusters, tolerance, init,
													yinyang_t, metric, average_distance, seed,
													device,	verbosity)
	
	return(result)
	
}


#' @title K-Nearest Neighbor Classification using CUDA
#' @description k-nearest neighbor classification using
#' a NVIDIA GPU via CUDA backend
#' @param k The number of neighbors to search for each sample
#' @param samples Numeric matrix
#' @param centroids Numeric matrix with precalculated clusters' 
#' centroids
#' @param assignments integer vector with sample-cluster 
#' associations. Indices start from 1.
#' @param metric character name of the distance metric to use.
#' The default is Euclidean (L2), it can be changed to 
#' "cos" for Sphereical K-means with angular distance. NOTE - 
#' the samples must be normalized in the latter case.
#' @param device integer defining device to use. 1 = first device,
#' 2 = second device, 3 = first & second devices, 0 = use all devices.
#' Default = 0
#' @param verbosity Integer indicating amount of output to see.
#' 0 = silence, 1 = progress logging, 2 = all output
#' @return Integer matrix with neighbor indices of shape [nsamp, k].
#' @export
knn_cuda <- function(
	k, samples, centroids, 
	assignments, metric = "L2", 
	device = 0, verbosity = 0)
{
	if(k <= 0 || !is.numeric(k)){
		stop("k must be a positive integer")
	}
	
	if(!is.numeric(centroids)){
		stop("centroids must be numeric")
	}
	
	if(length(assignments) != nrow(samples)){
		stop("invalid assignment's length")
	}
	
	if(is.list(samples)){
		samples <- do.call('rbind', samples)
	}
	
	result <- r_knn_cuda(k, samples, centroids,
											 assignments, metric,
											 device, verbosity)
	
	return(result)
	
}
