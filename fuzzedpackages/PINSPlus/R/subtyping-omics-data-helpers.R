FindMaxHeight <- function (hc, maxK) {
    height <- rev(diff(hc$height))[1:(maxK-1)]
    i=which.max(height)
    i+1
}

GetSimilarityFromGrouping <- function(g) {
    N = length(g)
    S = matrix(0, N, N)
    colnames(S) = rownames(S) = names(g)
    for (j in unique(g)) {
        X = rep(0, N);
        X[which(g == j)] = 1
        S = S + X %*% t(X)
    }
    S
}

CalcClusterAgreement <- function(l) {
    N <- length(l[[1]])
    A <- matrix(0, N, N)
    for (group in l) {
        A <- A + GetSimilarityFromGrouping(group)
    }
    A = A / (length(l))
    ret = (sum(A == 0) + sum(A == 1)) / (N ^ 2)
    
    ret
}

ClusteringFunWrapper <- function(kMax, groupings, FUN, k = NULL){
    clusters = list()
    
    if (is.null(k)){
        agreeS <- lapply(2 : kMax, FUN = function(i){
            cluster = FUN(i)
            clusters[[i]] <<- cluster
            
            CalcClusterAgreement(
                c(
                    groupings,
                    list(cluster)
                )
            )
        })
        
        ret <- list()
        ret$agreeS <- unlist(agreeS)
        ret$k <- which.max(agreeS) + 1
        ret$agree <- ret$agreeS[ret$k - 1]
        ret$cluster <- clusters[[ret$k]]
        return(ret)
    } else {
        ret <- list()
        ret$k = k
        ret$cluster <- FUN(k)
        ret$agree <- CalcClusterAgreement(
            c(
                groupings,
                list(ret$cluster)
            )
        )
        return(ret)
    }
}

ClusterUsingPAM <- function(orig, kMax, groupings, k = NULL) {
    ClusteringFunWrapper(kMax, groupings, FUN = function(k) pamWrapper(1 - orig, k, diss = T), k)
}

ClusterUsingHierarchical <- function(orig, kMax, groupings, k = NULL) {
    hcO <- hclust(as.dist(1 - orig), method = "average")
    ClusteringFunWrapper(kMax, groupings, FUN = function(k) cutree(hcO, k), k)
}
