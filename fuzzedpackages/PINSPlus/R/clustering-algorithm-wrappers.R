kmeansWrapper <- function(data, k, nstart = 20, iter.max = 1000, ...) {
    if (nstart > nrow(data)) nstart = nrow(data)
    kmeans(x = data, centers = k, nstart = nstart, iter.max = iter.max, ...)$cluster
}

pamWrapper <- function(data, k, ...){
    cluster::pam(data, k, cluster.only = T, ...)
}

hclustWrapper <- function(data, k, ...){
    cutree(hclust(dist(data), ...), k = k)
}