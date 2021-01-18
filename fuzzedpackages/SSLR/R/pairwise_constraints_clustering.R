

mustlink_cantlink <- function(y){
  mustLink = list()
  cantLink = list()

  i_mustkLink = 0
  i_cantLink = 0

  for(idx in 1:(length(y)-1)){

    if(!is.na(y[idx])){


      for(idx_dos in (idx + 1):length(y)){

        if(!is.na(y[idx_dos])){

          if(y[idx] == y[idx_dos]){

            i_mustkLink = i_mustkLink + 1

            mustLink[[i_mustkLink]] = c(idx, idx_dos)
          }
          else{
            i_cantLink = i_cantLink + 1

            cantLink[[i_cantLink]] = c(idx, idx_dos)
          }

        }



      }

    }

  }


  mustLink <- do.call(rbind, mustLink)
  cantLink <- do.call(rbind, cantLink)

  list("mustLink" = mustLink, "cantLink" = cantLink)
}



#' @title  General Interface COP K-Means Algorithm
#' @description Model from conclust \cr
#' This function takes an unlabeled dataset and two lists of must-link and cannot-link constraints
#' as input and produce a clustering as output.
#' @note This models only returns labels, not centers
#' @param n_clusters A number of clusters to be considered. Default is NULL (num classes)
#' @param mustLink A list of must-link constraints. NULL Default, constrints same label
#' @param cantLink A list of cannot-link constraints. NULL Default, constrints with different label
#' @param max_iter maximum iterations in KMeans. Default is 10
#' @references
#' Wagstaff, Cardie, Rogers, Schrodl\cr
#' \emph{Constrained K-means Clustering with Background Knowledge}\cr
#' 2001
#' @example demo/ckmeansSSLR.R
#' @importFrom conclust ckmeans
#' @export
ckmeansSSLR <- function(n_clusters = NULL, mustLink = NULL, cantLink = NULL,
                          max_iter = 10){


  train_function <- function(x, y) {

    #Load conclust
    load_conclust()


    x <- as.data.frame(x)
    y <- as.factor(y)

    if(is.null(n_clusters)){
      n_clusters <- length(levels(y))
    }

    if(is.null(mustLink) || is.null(cantLink)){
      list_mustlink_cantlink <- mustlink_cantlink(y)
      mustLink <- list_mustlink_cantlink$mustLink
      cantLink <-list_mustlink_cantlink$cantLink
    }


    labels <- ckmeans(x, k = n_clusters, mustLink = mustLink,
                    cantLink = cantLink,maxIter = max_iter)

    totss <- sum(scale(x, scale = FALSE)^2)

    result <- structure(list(cluster = labels, centers = NULL, totss = totss),
              class = "kmeans")

    ### Result ###
    result$classes = levels(y)
    #result$pred.params = c("class","raw", "prob")
    result$mode = "clustering"
    #class(result) <- "kmeans"

    return(result)

  }

  args <- list(
    n_clusters = n_clusters,
    mustLink = mustLink,
    cantLink = cantLink,
    max_iter = max_iter
  )

  new_model_sslr(train_function, "ckmeansSSLR", args)

}




#' @title  General Interface Pairwise Constrained Clustering By Local Search
#' @description Model from conclust \cr
#' This function takes an unlabeled dataset and two lists of must-link and cannot-link constraints
#' as input and produce a clustering as output.
#' @note This models only returns labels, not centers
#' @param n_clusters A number of clusters to be considered. Default is NULL (num classes)
#' @param mustLink A list of must-link constraints. NULL Default, constrints same label
#' @param cantLink A list of cannot-link constraints. NULL Default, constrints with different label
#' @param max_iter maximum iterations in KMeans. Default is 1
#' @param tabuIter Number of iteration in Tabu search
#' @param tabuLength The number of elements in the Tabu list
#' @references
#' Tran Khanh Hiep, Nguyen Minh Duc, Bui Quoc Trung\cr
#' \emph{Pairwise Constrained Clustering by Local Search}\cr
#' 2016
#' @example demo/cclsSSLR.R
#' @importFrom conclust ccls
#' @export
cclsSSLR <- function(n_clusters = NULL, mustLink = NULL, cantLink = NULL,
                      max_iter = 1, tabuIter = 100, tabuLength = 20){


  train_function <- function(x, y) {

    #Load conclust
    load_conclust()


    x <- as.data.frame(x)
    y <- as.factor(y)

    if(is.null(n_clusters)){
      n_clusters <- length(levels(y))
    }

    if(is.null(mustLink) || is.null(cantLink)){
      list_mustlink_cantlink <- mustlink_cantlink(y)
      mustLink <- list_mustlink_cantlink$mustLink
      cantLink <-list_mustlink_cantlink$cantLink
    }


    labels <- ccls(x, k = n_clusters, mustLink = mustLink,
                      cantLink = cantLink,maxIter = max_iter,
                      tabuIter = tabuIter, tabuLength = tabuLength)

    totss <- sum(scale(x, scale = FALSE)^2)

    result <- structure(list(cluster = labels, centers = NULL, totss = totss),
                        class = "kmeans")

    ### Result ###
    result$classes = levels(y)
    #result$pred.params = c("class","raw", "prob")
    result$mode = "clustering"
    #class(result) <- "kmeans"

    return(result)

  }

  args <- list(
    n_clusters = n_clusters,
    mustLink = mustLink,
    cantLink = cantLink,
    max_iter = max_iter,
    tabuIter = tabuIter,
    tabuLength = tabuLength
  )

  new_model_sslr(train_function, "cclsSSLR", args)

}



#' @title  General Interface MPC K-Means Algorithm
#' @description Model from conclust \cr
#' This function takes an unlabeled dataset and two lists of must-link and cannot-link constraints
#' as input and produce a clustering as output.
#' @note This models only returns labels, not centers
#' @param n_clusters A number of clusters to be considered. Default is NULL (num classes)
#' @param mustLink A list of must-link constraints. NULL Default, constrints same label
#' @param cantLink A list of cannot-link constraints. NULL Default, constrints with different label
#' @param max_iter maximum iterations in KMeans. Default is 10
#' @references
#' Bilenko, Basu, Mooney \cr
#' \emph{Integrating Constraints and Metric Learning in Semi-Supervised Clustering}\cr
#' 2004
#' @example demo/mpckmSSLR.R
#' @importFrom conclust mpckm
#' @export
mpckmSSLR <- function(n_clusters = NULL, mustLink = NULL, cantLink = NULL,
                        max_iter = 10){


  train_function <- function(x, y) {

    #Load conclust
    load_conclust()


    x <- as.data.frame(x)
    y <- as.factor(y)

    if(is.null(n_clusters)){
      n_clusters <- length(levels(y))
    }

    if(is.null(mustLink) || is.null(cantLink)){
      list_mustlink_cantlink <- mustlink_cantlink(y)
      mustLink <- list_mustlink_cantlink$mustLink
      cantLink <-list_mustlink_cantlink$cantLink
    }


    labels <- mpckm(x, k = n_clusters, mustLink = mustLink,
                      cantLink = cantLink,maxIter = max_iter)

    totss <- sum(scale(x, scale = FALSE)^2)

    result <- structure(list(cluster = labels, centers = NULL, totss = totss),
                        class = "kmeans")

    ### Result ###
    result$classes = levels(y)
    #result$pred.params = c("class","raw", "prob")
    result$mode = "clustering"
    #class(result) <- "kmeans"

    return(result)

  }

  args <- list(
    n_clusters = n_clusters,
    mustLink = mustLink,
    cantLink = cantLink,
    max_iter = max_iter
  )

  new_model_sslr(train_function, "mpckmSSLR", args)

}




#' @title  General LCVQE Algorithm
#' @description Model from conclust \cr
#' This function takes an unlabeled dataset and two lists of must-link and cannot-link constraints
#' as input and produce a clustering as output.
#' @note This models only returns labels, not centers
#' @param n_clusters A number of clusters to be considered. Default is NULL (num classes)
#' @param mustLink A list of must-link constraints. NULL Default, constrints same label
#' @param cantLink A list of cannot-link constraints. NULL Default, constrints with different label
#' @param max_iter maximum iterations in KMeans. Default is 2
#' @references
#' Dan Pelleg, Dorit Baras\cr
#' \emph{K-means with large and noisy constraint sets}\cr
#' 2007
#' @example demo/lcvqeSSLR.R
#' @importFrom conclust lcvqe
#' @export
lcvqeSSLR <- function(n_clusters = NULL, mustLink = NULL, cantLink = NULL,
                      max_iter = 2){


  train_function <- function(x, y) {

    #Load conclust
    load_conclust()


    x <- as.data.frame(x)
    y <- as.factor(y)

    if(is.null(n_clusters)){
      n_clusters <- length(levels(y))
    }

    if(is.null(mustLink) || is.null(cantLink)){
      list_mustlink_cantlink <- mustlink_cantlink(y)
      mustLink <- list_mustlink_cantlink$mustLink
      cantLink <-list_mustlink_cantlink$cantLink
    }


    labels <- lcvqe(x, k = n_clusters, mustLink = mustLink,
                    cantLink = cantLink,maxIter = max_iter)

    totss <- sum(scale(x, scale = FALSE)^2)

    result <- structure(list(cluster = labels, centers = NULL, totss = totss),
                        class = "kmeans")

    ### Result ###
    result$classes = levels(y)
    #result$pred.params = c("class","raw", "prob")
    result$mode = "clustering"
    #class(result) <- "kmeans"

    return(result)

  }

  args <- list(
    n_clusters = n_clusters,
    mustLink = mustLink,
    cantLink = cantLink,
    max_iter = max_iter
  )

  new_model_sslr(train_function, "lcvqeSSLR", args)

}
