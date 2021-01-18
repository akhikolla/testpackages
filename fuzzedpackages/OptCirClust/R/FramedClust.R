# FramedCluster.R
#'
#'
#'
#' Framed Data Clustering
#'
#' Find a frame of given size, among all possible such frames
#' on the input data, to minimize the minimum within-cluster
#' sum of square distances.
#'
#' @import stats
#'
#' @param X a vector of data points to perform framed clustering
#' @param K the number of clusters in each frame
#' @param frame.size the number of points from X to be included in each frame.
#' It is not the width of the frame.
#' @param first.frame starting index of the first frame to be clustered.
#' The first point in the first frame is \code{X[first.frame]}.
#' @param last.frame starting index of the last frame to be clustered.
#' The first point in the first frame is \code{X[last.frame]}.
#' @param method the framed clustering method. See Details.
#'
#'
#' @details
#' The method option \code{"linear.polylog"} (default) performs
#'  fast optimal framed clustering. The runtime is
#'  \eqn{O(K N \log^2 N)}{O(K N log^2 N)}.
#'
#' The \code{"kmeans"} option repeatedly calling the heuristic
#'  k-means algorithm in all frames without any guarantee of
#'  cluster optimality.
#'
#' The method option \code{"Ckmeans.1d.dp"} performs optimal framed
#'  clustering by repeatedly finding the best clustering within
#'  each frame using the \code{"Ckmeans.1d.dp"} method. At a runtime
#'  of \eqn{O(K N^2)}{O(K N^2)}, the algorithm is slow but optimal.
#'  It is included to provide a baseline.
#'
#'
#' @return An object of class \code{"FramedClust"} which has a \code{plot}
#'  method. It is a list with the following components:
#'
#' \item{cluster}{a vector of clusters assigned to each element in x.
#' Each cluster is indexed by an integer from 1 to K.
#' NA represents points from X that are outside the optimal frame,
#' thus not part of any cluster.}
#'
#' \item{centers}{a numeric vector of the means for each cluster
#' in the frame.}
#'
#' \item{withinss}{a numeric vector of the within-cluster sum of
#' squared distances for each cluster.}
#'
#' \item{size}{a vector of the number of elements in each cluster.}
#'
#' \item{totss}{total sum of squared distances between each element
#'  and the sample mean. This statistic is not dependent on
#'  the clustering result.}
#'
#' \item{tot.withinss}{total sum of within-cluster squared
#' distances between each element and its cluster mean.
#' This statistic is minimized given the number of clusters.}
#'
#' \item{betweenss}{sum of squared distances between each
#' cluster mean and sample mean. This statistic is maximized
#' given the number of clusters.}
#'
#' \item{X_name}{a character string. The actual name of
#' the \code{X} argument.}
#'
#'@examples
#' N <- 100
#' X <- rnorm(N)
#' K <- 5
#' frame.size <- 60
#'
#' result <- FramedClust(X, K, frame.size)
#' plot(result, main="Example 1. Framed clustering on all frames")
#'
#' frame.size <- 40
#' first.frame <- 30
#' last.frame <- 50
#' method <- "linear.polylog"
#'
#' result <- FramedClust(X, K, frame.size, first.frame,
#'                       last.frame, method)
#' plot(result, main="Example 2. Framed clustering on a subset of frames")

#' @export
FramedClust <- function(
  X, K, frame.size,
  first.frame = 1,
  last.frame = length(X)-frame.size+1,
  method = c("linear.polylog", "kmeans", "Ckmeans.1d.dp")
)
{
  if(K < 1 || K > frame.size) {
    stop("Invalid number of clusters K")
  }

  if(first.frame < 1 ||
     first.frame > last.frame ||
     last.frame > length(X)-frame.size+1) {
    stop("Invalid first or last frame value")
  }

  X_name <- deparse(substitute(X))

  method <- match.arg(method)

  I <- order(X)

  X_sort <- X[I]

  prev_k_f = -1
  next_k_f = -1


  if(method == "linear.polylog")
  {
    result <- lin_polylog_framed_clust(
      as.double(X_sort), as.integer(K), as.integer(frame.size),
      as.integer(first.frame-1), as.integer(last.frame-1),
      as.integer(prev_k_f), as.integer(next_k_f))

    cluster <-  rep(1, (result$Border[1] - result$ID + 1))

    if(K > 1)
    {
      for (i in 2:K)
      {
        cluster <-
          c(cluster, rep(i, (result$Border[i] - result$Border[i - 1])))
      }

      cluster_new <- matrix( 0, nrow = 1, ncol = length(X_sort))


      cluster_new[(result$ID + 1):(result$ID + frame.size )] <- cluster


      cluster[I] <- cluster_new


      cluster[which(cluster==0)] <- NA


    }


  } else if(method == "kmeans"){

    result <- kmeans.framed.clust(X_sort, K, frame.size, as.integer(first.frame-1), as.integer(last.frame-1))

    cluster <- result$cluster

    cluster_new <- matrix( 0, nrow = 1, ncol = length(X_sort))


    cluster_new[(result$ID + 1):(result$ID + frame.size )] <- cluster


    cluster[I] <- cluster_new

    cluster[which(cluster==0)] <- NA

  } else if(method == "Ckmeans.1d.dp"){

    result <- quad.framed.clust(X_sort, K, frame.size, as.integer(first.frame-1), as.integer(last.frame-1))

    cluster <- result$cluster

    cluster_new <- matrix( 0, nrow = 1, ncol = length(X_sort))


    cluster_new[(result$ID + 1):(result$ID + frame.size )] <- cluster


    cluster[I] <- cluster_new

    cluster[which(cluster==0)] <- NA

  }


  # if(!is.null(first.frame_old))
  #  {
  #   cluster <- c(rep(NA,(first.frame_old-1)),cluster)
  # }

  # Border.mid <- (X[result$Border + 1] + X[result$Border + 2]) / 2

  df <-
    list(

      "cluster" = cluster,
      "centers" = result$centers + 1,
      "withinss" = result$withinss,
      "size" = result$size,
      "totss" = result$totss,
      "tot.withinss" = result$tot.withinss,
      "betweenss" = result$betweenss,
      "X_name" = X_name

    )

  class (df) <- "FramedClust"

  return(df)
}

