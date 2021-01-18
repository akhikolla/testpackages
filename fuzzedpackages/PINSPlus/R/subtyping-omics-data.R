#' @title Subtyping multi-omics data
#' @description Perform subtyping using multiple types of data
#' 
#' @param dataList a list of data matrices. Each matrix represents a data type where the rows are items and the columns are features. The matrices must have the same set of items.
#' @param kMax The maximum number of clusters used for automatically detecting the number of clusters in \code{PerturbationClustering}. This paramter is passed to \code{PerturbationClustering} and does not affect the final number of cluster in \code{SubtypingOmicsData}. Default value is \code{5}.
#' @param kMin The minimum number of clusters used for automatically detecting the number of clusters in \code{PerturbationClustering}. This paramter is passed to \code{PerturbationClustering} and does not affect the final number of cluster in \code{SubtypingOmicsData}. Default value is \code{2}.
#' @param k The number of clusters. If k is set then kMin and kMax will be ignored.
#' @param agreementCutoff agreement threshold to be considered consistent. Default value is \code{0.5}.
#' @param ncore Number of cores that the algorithm should use. Default value is \code{1}.
#' @param verbose set it to \code{TRUE} of \code{FALSE} to get more or less details respectively.
#' @param ... these arguments will be passed to \code{PerturbationClustering} algorithm. See details for more information
#' 
#' @details
#' 
#' \code{SubtypingOmicsData} implements the Subtyping multi-omic data that are based on Perturbaion clustering algorithm of Nguyen, et al (2017) and Nguyen, et al (2019).
#' The input is  a list of data matrices where each matrix represents the molecular measurements of a data type. The input matrices must have the same number of rows. 
#' \code{SubtypingOmicsData} aims to find the optimum number of subtypes and location of each sample in the clusters from integrated input data \code{dataList} through two processing stages:
#' 
#' 1. Stage I: The algorithm first partitions each data type using the function \code{PerturbationClustering}.
#' It then merges the connectivities across data types into similarity matrices.
#' Both kmeans and similarity-based clustering algorithms - partitioning around medoids \code{pam} are used to partition the built similarity.
#' The algorithm returns the partitioning that agrees the most with individual data types.\cr
#' 2. Stage II: The algorithm attempts to split each discovered group if there is a strong agreement between data types,
#' or if the subtyping in Stage I is very unbalanced.
#'
#' @return
#' 
#' \code{SubtypingOmicsData} returns a list with at least the following components:
#' \item{cluster1}{A vector of labels indicating the cluster to which each sample is allocated in Stage I}
#' \item{cluster2}{A vector of labels indicating the cluster to which each sample is allocated in Stage II}
#' \item{dataTypeResult}{A list of results for individual data type. Each element of the list is the result of the \code{PerturbationClustering} for the corresponding data matrix provided in dataList.}
#' 
#' 
#' @references
#' 
#' 1. H Nguyen, S Shrestha, S Draghici, & T Nguyen. PINSPlus: a tool for tumor subtype discovery in integrated genomic data. Bioinformatics, 35(16), 2843-2846, (2019).
#' 
#' 2. T Nguyen, R Tagett, D Diaz, S Draghici. A novel method for data integration and disease subtyping. Genome Research, 27(12):2025-2039, 2017.
#' 
#' 3. T. Nguyen, "Horizontal and vertical integration of bio-molecular data", PhD thesis, Wayne State University, 2017.
#' 
#' @seealso \code{\link{PerturbationClustering}}
#' 
#' @examples
#' \donttest{
#' # Load the kidney cancer carcinoma data
#' data(KIRC)
#' 
#' # Perform subtyping on the multi-omics data
#' dataList <- list (as.matrix(KIRC$GE), as.matrix(KIRC$ME), as.matrix(KIRC$MI)) 
#' names(dataList) <- c("GE", "ME", "MI")
#' result <- SubtypingOmicsData(dataList = dataList)
#' 
#' # Change Pertubation clustering algorithm's arguments
#' result <- SubtypingOmicsData(
#'     dataList = dataList, 
#'     clusteringMethod = "kmeans", 
#'     clusteringOptions = list(nstart = 50)
#' )
#'
#' # Plot the Kaplan-Meier curves and calculate Cox p-value
#' library(survival)
#' cluster1=result$cluster1;cluster2=result$cluster2
#' a <- intersect(unique(cluster2), unique(cluster1))
#' names(a) <- intersect(unique(cluster2), unique(cluster1))
#' a[setdiff(unique(cluster2), unique(cluster1))] <- seq(setdiff(unique(cluster2), unique(cluster1))) 
#'                                                   + max(cluster1)
#' colors <- a[levels(factor(cluster2))]
#' coxFit <- coxph(
#'  Surv(time = Survival, event = Death) ~ as.factor(cluster2),
#'  data = KIRC$survival,
#'  ties = "exact"
#' )
#' mfit <- survfit(Surv(Survival, Death == 1) ~ as.factor(cluster2), data = KIRC$survival)
#' plot(
#'  mfit, col = colors,
#'  main = "Survival curves for KIRC, level 2",
#'  xlab = "Days", ylab = "Survival",lwd = 2
#' )
#' legend("bottomright", 
#'     legend = paste(
#'         "Cox p-value:", 
#'         round(summary(coxFit)$sctest[3], digits = 5), 
#'         sep = ""
#'     )
#' )
#' legend(
#'     "bottomleft",
#'     fill = colors,
#'     legend = paste(
#'         "Group ",
#'         levels(factor(cluster2)),": ", table(cluster2)[levels(factor(cluster2))], 
#'         sep =""
#'     )
#' )
#' 
#' }
#' @importFrom FNN knnx.index
#' @importFrom entropy entropy
#' @export
SubtypingOmicsData <- function (dataList, kMin = 2, kMax = 5, k = NULL, agreementCutoff = 0.5, ncore = 1, verbose = T, ...) {
    now = Sys.time()
    
    # defined log function
    mlog <- if(!verbose) function(...){} else function(...){
        message(...)
        flush.console()
    }
    
    dataListTrain <- NULL
    dataListTest <- NULL
    
    seed = round(rnorm(1)*10^6)
    
    if (nrow(dataList[[1]]) > 2000) {
        n_samples <- nrow(dataList[[1]])
        ind <- sample.int(n_samples, size = 2000)
        dataListTrain <- lapply(dataList, function(x) x[ind, ])
        dataListTest <- lapply(dataList, function(x) x[-ind, , drop=F])
        
        dataList <- dataListTrain
    }
    
    runPerturbationClustering <- function(dataList, kMin, kMax, stage = 1, forceSplit = FALSE, k = NULL){
        dataTypeResult <- lapply(dataList, function(data) {
            set.seed(seed)
            PerturbationClustering(data, kMin, kMax, ncore = ncore, verbose = verbose,...)
        })
        
        origList <- lapply(dataTypeResult, function(r) r$origS[[r$k]])
        orig = Reduce('+', origList)/length(origList)
        PW = Reduce('*', origList)
        agreement = (sum(orig == 0) + sum(orig == 1) - nrow(orig)) / (nrow(orig) ^ 2 - nrow(orig))
        
        pert =  Reduce('+', lapply(dataTypeResult, function(r) r$pertS[[r$k]]))/length(dataList)
        
        groups <- NULL
        
        mlog("STAGE : ", stage, "\t Agreement : ", agreement)
        
        if (agreement >= agreementCutoff | forceSplit){
            hcW <- hclust(dist(PW))
            
            maxK = min(kMax*2, dim(unique(PW, MARGIN = 2))[2] - (stage - 1))
            maxHeight = FindMaxHeight(hcW, maxK = min(2*maxK, 10))
            groups <- cutree(hcW, maxHeight)
            
            # if k is specific then only use that k if the number of groups > k
            # the max(groups) < k, this may cause small groups
            if (!is.null(k) && max(groups) > k){
                groups <- cutree(hcW, k)
            }
        }
        
        list(dataTypeResult = dataTypeResult, orig = orig, pert = pert, PW = PW, groups = groups, agreement = agreement)
    }
    
    pResult <- runPerturbationClustering(dataList, kMin, kMax, k = k)
    
    groups <- pResult$groups
    groups2 <- NULL
    
    if (!is.null(groups)) {
        groups2 <- groups
        
        if (is.null(k)){
            for (g in sort(unique(groups))) {
                miniGroup <- names(groups[groups == g])
                if (length(miniGroup) > 30) {
                    groupsM <- runPerturbationClustering(dataList = lapply(dataList, function(d) d[miniGroup, ]), kMin = kMin, kMax = min(kMax, 5), stage = 2)$groups
                    if (!is.null(groupsM))
                        groups2[miniGroup] <- paste(g, groupsM, sep = "-")
                }
            }
        } else {
            
            agreements <- rep(1, length(groups))
            names(agreements) <- names(groups)
            
            #  if k is specific then force further split
            tbl <- sort(table(groups2), decreasing = T)
            minGroupSize <- 30
            while (length(unique(groups2)) < k){
                if (all(tbl <= 30)) {
                    minGroupSize <- 10
                }
                
                # cannot split anymore
                if (all(tbl <= 10)) break()
                
                for (g in names(tbl)){
                    miniGroup <- names(groups2[groups2 == g])
                    if (length(miniGroup) > minGroupSize) {
                        splitRes <- runPerturbationClustering(dataList = lapply(dataList, function(d) d[miniGroup, ]), kMin = kMin, kMax = min(kMax, 5), stage = 2, forceSplit = T)
                        groupsM <- splitRes$groups
                        
                        if (!is.null(groupsM)){
                            groups2[miniGroup] <- paste(g, groupsM, sep = "-")
                            agreements[miniGroup] <- splitRes$agreement
                        }
                    }
                }
                tbl <- sort(table(groups2), decreasing = T)
            }
            
            # now after further splitting, the number of groups can be > k
            # need to merge cluster based on their aggrement
            agreements.unique = unique(agreements)
            for (aggr in sort(unique(agreements))){
                if (length(unique(groups2)) == k) break()
                merge.group <- agreements == aggr
                
                k.smallGroup <- length(unique(groups2[merge.group]))
                k.need <- k - (length(unique(groups2)) - k.smallGroup + 1) + 1
                
                groups2[merge.group] <- unlist(lapply(strsplit(groups2[merge.group], "-"), function(g){
                    paste0(g[1:(length(g)-1)], collapse = "-")
                }))
                
                if (k.need > 1){
                    splitRes <- runPerturbationClustering(dataList = lapply(dataList, function(d) d[miniGroup, ]), 
                                                          kMin = kMin, kMax = min(kMax, 5), stage = 2, forceSplit = T, 
                                                          k = k.need)
                    groupsM <- splitRes$groups
                    
                    groups2[merge.group] <- paste(groups2[merge.group], groupsM, sep = "-")
                }
                
                agreements[merge.group] <- 1
            }
        }
    } else{
        set.seed(seed)
        
        orig <- pResult$orig
        dataTypeResult <- pResult$dataTypeResult
        clusteringAlgorithm = GetClusteringAlgorithm(...)$fun
        
        groupings <- lapply(dataTypeResult, function(r) clusteringAlgorithm(data = r$origS[[r$k]], k = r$k))
        
        pGroups <- ClusterUsingPAM(orig = orig, kMax = kMax*2, groupings = groupings)
        hGroups <- ClusterUsingHierarchical(orig = orig, kMax = kMax*2, groupings = groupings)
        
        pAgree  = pGroups$agree; hAgree  = hGroups$agree;
        
        groups <- (if (pAgree > hAgree) pGroups else if (hAgree > pAgree) hGroups else {
            pAgree = ClusterUsingPAM(orig = pResult$pert, kMax = kMax, groupings = groupings)$agree
            hAgree = ClusterUsingHierarchical(orig = pResult$pert, kMax = kMax, groupings = groupings)$agree
            if (hAgree - pAgree >= 1e-3) hGroups else pGroups
        })$cluster
        
        names(groups) <- rownames(orig)
        groups2 <- groups
        
        if (is.null(k)){
            mlog("Check if can proceed to stage II")
            
            normalizedEntropy = entropy::entropy(table(groups)) / log(length(unique(groups)), exp(1))
            
            if (normalizedEntropy < 0.5) {
                for (g in sort(unique(groups))) {
                    miniGroup <- names(groups[groups == g])
                    #this is just to make sure we don't split a group that is already very small
                    if (length(miniGroup) > 30) {
                        #this is to check if the data types in this group can be split
                        groupsM <- runPerturbationClustering(dataList = lapply(dataList, function(d) d[miniGroup, ]),kMin = kMin, kMax = min(kMax, 5), stage = 2, forceSplit = T, k = NULL)$groups
                        if (!is.null(groupsM))
                            groups2[miniGroup] <- paste(g, groupsM, sep = "-")
                    }
                }
            }
        } else {
            if (length(unique(groups2)) > k){
                pGroups <- ClusterUsingPAM(orig = orig, kMax = kMax*2, groupings = groupings, k)
                hGroups <- ClusterUsingHierarchical(orig = orig, kMax = kMax*2, groupings = groupings, k)
                
                pAgree  = pGroups$agree; hAgree  = hGroups$agree;
                
                groups <- (if (pAgree > hAgree) pGroups else if (hAgree > pAgree) hGroups else {
                    pAgree = ClusterUsingPAM(orig = pResult$pert, kMax = kMax, groupings = groupings, k)$agree
                    hAgree = ClusterUsingHierarchical(orig = pResult$pert, kMax = kMax, groupings = groupings, k)$agree
                    if (hAgree - pAgree >= 1e-3) hGroups else pGroups
                })$cluster
                
                names(groups) <- rownames(orig)
                groups2 <- groups
            } else if (length(unique(groups2)) < k){
                
                # split like normal using entropy
                normalizedEntropy = entropy::entropy(table(groups)) / log(length(unique(groups)), exp(1))

                agreements <- rep(1, length(groups))
                names(agreements) <- names(groups)
                
                if (normalizedEntropy < 0.5) {
                    for (g in sort(unique(groups))) {
                        miniGroup <- names(groups[groups == g])
                        #this is just to make sure we don't split a group that is already very small
                        if (length(miniGroup) > 30) {
                            #this is to check if the data types in this group can be split
                            splitRes <- runPerturbationClustering(dataList = lapply(dataList, function(d) d[miniGroup, ]), kMin = kMin, kMax = min(kMax, 5), stage = 2, forceSplit = T)
                            groupsM <- splitRes$groups
                            
                            if (!is.null(groupsM)){
                                groups2[miniGroup] <- paste(g, groupsM, sep = "-")
                                agreements[miniGroup] <- splitRes$agreement
                            }
                        }
                    }
                }
                
                # if the number of group is still less than k, then force split
                if (length(unique(groups2)) < k){
                    
                    #  if k is specific then force further split
                    tbl <- sort(table(groups2), decreasing = T)
                    minGroupSize <- 30
                    while (length(unique(groups2)) < k){
                        if (all(tbl <= 30)) {
                            minGroupSize <- 10
                        }
                        
                        # cannot split anymore
                        if (all(tbl <= 10)) break()
                        
                        for (g in names(tbl)){
                            miniGroup <- names(groups2[groups2 == g])
                            if (length(miniGroup) > minGroupSize) {
                                splitRes <- runPerturbationClustering(dataList = lapply(dataList, function(d) d[miniGroup, ]), kMin = kMin, kMax = min(kMax, 5), stage = 2, forceSplit = T)
                                groupsM <- splitRes$groups
                                
                                if (!is.null(groupsM)){
                                    groups2[miniGroup] <- paste(g, groupsM, sep = "-")
                                    agreements[miniGroup] <- splitRes$agreement
                                }
                            }
                        }
                        tbl <- sort(table(groups2), decreasing = T)
                    }
                }
                
                # now after further splitting, the number of groups can be > k
                # need to merge cluster based on their aggrement
                agreements.unique = unique(agreements)
                for (aggr in sort(unique(agreements))){
                    if (length(unique(groups2)) == k) break()
                    merge.group <- agreements == aggr
                    
                    k.smallGroup <- length(unique(groups2[merge.group]))
                    k.need <- k - (length(unique(groups2)) - k.smallGroup + 1) + 1
                    
                    groups2[merge.group] <- unlist(lapply(strsplit(groups2[merge.group], "-"), function(g){
                        paste0(g[1:(length(g)-1)], collapse = "-")
                    }))
                    
                    if (k.need > 1){
                        splitRes <- runPerturbationClustering(dataList = lapply(dataList, function(d) d[merge.group, ]), 
                                                              kMin = kMin, kMax = min(kMax, 5), stage = 2, forceSplit = T, 
                                                              k = k.need)
                        groupsM <- splitRes$groups
                        
                        groups2[merge.group] <- paste(groups2[merge.group], groupsM, sep = "-")
                    }
                    
                    agreements[merge.group] <- 1
                }
            }
        }
    }
    
    train_y <- groups
    train_y2 <- groups2
    
    if(!is.null(dataListTest)) {
        set.seed(seed)
        RcppParallel::setThreadOptions(ncore)
        
        test_prob <- matrix(0, nrow = n_samples - 2000, ncol = length(unique(groups)))
        if(!is.null(train_y2))
        {
            test_prob2 <- matrix(0, nrow = n_samples - 2000, ncol = length(unique(groups2)))
        }
        for (i in 1:length(dataListTrain)) {
            train <- dataListTrain[[i]]
            test <- dataListTest[[i]]
            
            if(ncol(train)*nrow(train) > 2e7) {
                pca <- pResult$dataTypeResult[[i]]$pca
            } else {
                pca <- rpca.para(train, min(nrow(data), 20), scale = F)
            }

            train <- pca$x
            test <- predict.rpca.para(pca, test)
            
            nn_index <- FNN::knnx.index(train, test, k = 10)
            
            nn_group <- matrix(train_y[nn_index], ncol = 10)
            
            for (j in 1:nrow(test_prob)) {
                tmp <- table(nn_group[j, ])
                test_prob[j, as.numeric(names(tmp))] <- test_prob[j, as.numeric(names(tmp))] + as.numeric(tmp)/10
            }
            
            if(!is.null(train_y2))
            {
                nn_group2 <- matrix(train_y2[nn_index], ncol = 10)
                
                for (j in 1:nrow(test_prob)) {
                    tmp <- table(nn_group2[j, ])
                    test_prob2[j, as.numeric(names(tmp))] <- test_prob2[j, as.numeric(names(tmp))] + as.numeric(tmp)/10
                }
            }
        }
        
        test_y <- apply(test_prob, 1, which.max)
        
        groups <- rep(0, n_samples)
        groups[ind] <- train_y
        groups[-ind] <- test_y
        
        if(!is.null(train_y2)) {
            test_y2 <- apply(test_prob2, 1, which.max)
            groups2 <- rep(0, n_samples)
            groups2[ind] <- train_y2
            groups2[-ind] <- test_y2
        }
    }
    
    timediff = Sys.time() - now;
    mlog("Done in ", timediff, " ", units(timediff), ".\n")
    
    list(
        cluster1 = groups,
        cluster2 = groups2,
        dataTypeResult = pResult$dataTypeResult
    )
    
}