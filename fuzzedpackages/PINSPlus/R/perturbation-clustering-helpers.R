GetClusteringAlgorithm <- function(clusteringMethod = "kmeans", clusteringFunction = NULL, clusteringOptions = NULL, ...){
    name = clusteringMethod
    
    if (!is.function(clusteringFunction)) {
        switch(
            clusteringMethod,
            kmeans = {
                clusteringFunction <- kmeansWrapper
            },
            pam = {
                clusteringFunction <- pamWrapper
            },
            hclust = {
                clusteringFunction <- hclustWrapper
            },
            {
                stop("clusteringMethod not found. Please pass a clusteringFunction instead")
            }
        )
    }
    else {
        name = "Unknow"
    }
    
    list(
        fun = function(data, k) do.call(clusteringFunction, c(list(data = data, k = k), clusteringOptions)),
        name = name
    )
}

GetPerturbationAlgorithm <- function(data = data, perturbMethod = "noise", perturbFunction = NULL, perturbOptions = NULL, ...){
    name = perturbMethod
    
    if (!is.function(perturbFunction)) {
        switch(
            perturbMethod,
            noise = {
                if (is.null(perturbOptions))
                    perturbOptions <- list()
                
                noise = perturbOptions$noise
                if (is.null(noise))
                    noise <- GetNoise(data, noisePercent = perturbOptions$noisePercent)
                
                # get perturbed similarity
                perturbFunction <- function(data, ...) {
                    AddNoisePerturb(data = data, noise = noise)
                }
            },
            subsampling = {
                perturbFunction <- SubSampling
            },
            {
                stop("perturbMethod not found. Please pass a perturbFunction instead")
            }
        )
    }
    else {
        name = "Unknow"
    }
    
    list(
        fun = function(data) do.call(perturbFunction, c(list(data = data), perturbOptions)),
        name = name
    )
}

BuildConnectivityMatrix <- function(data, clus, clusteringAlgorithm) {
    rowNum <- nrow(data)
    S <- matrix(0L, rowNum, rowNum)
    cluster <- clusteringAlgorithm(data, clus)
    
    for (j in 1:clus) {
        idx <- cluster == j
        S[idx, idx] <- 1L
    }
    
    rownames(S) <- rownames(data)
    colnames(S) <- rownames(data)
    
    list(matrix = S, groups = cluster)
}

ClusterToConnectivity <- function(cluster){
    S <- matrix(0L, length(cluster), length(cluster))
    for (j in 1:max(cluster)) {
        idx <- cluster == j
        S[idx, idx] <- 1L
    }
    S
}

GetOriginalSimilarity <- function(data, clusRange, clusteringAlgorithm, showProgress = F, ncore) {
    groupings <- list()
    origS <- list()
    
    if (showProgress) {
        pb <- txtProgressBar(min = 0, max = length(clusRange), style = 3)
    }
    
    seeds = abs(round(rnorm(max(clusRange))*10^6))
    
    for (clus in clusRange){
        set.seed(seeds[clus])
        groupings[[clus]] <- clusteringAlgorithm(data, clus)
        origS[[clus]] <- ClusterToConnectivity(groupings[[clus]])
        rownames(origS[[clus]]) <- colnames(origS[[clus]]) <- rownames(data)
        if (showProgress) setTxtProgressBar(pb, clus)
    }
    
    list(origS = origS, groupings = groupings)
}

GetPerturbedSimilarity <- function(data, clusRange, iterMax, iterMin, origS, clusteringAlgorithm, perturbedFunction, stoppingCriteriaHandler, showProgress = F, ncore) {
    pertS <- list()
    currentIter <- rep(0,max(clusRange))

    jobs <- rep(clusRange, iterMax)
    maxJob = length(jobs)

    seeds = list()

    for (clus in clusRange) {
        pertS[[clus]] <- list()
        seeds[[clus]] <- round(rnorm(max(iterMax,1000))*10^6)
    }

    if (showProgress) pb <- txtProgressBar(min = 0, max = maxJob, style = 3)

    kProgress = rep(0, max(clusRange))
    step = iterMin
    if(.Platform$OS.type == "unix") doParallel::registerDoParallel(min(ncore, step))
    while (length(jobs) > 0) {
        jobLength = step*length(clusRange)
        step = 10

        if (jobLength > length(jobs)){
            jobLength =  length(jobs)
        }

        currentJobs = list()
        count = 1
        for (clus in jobs[1:jobLength]){
            kProgress[clus] = kProgress[clus] + 1
            currentJobs[[count]] <- list(iter = kProgress[clus], clus = clus)
            count = count + 1
        }

        if (jobLength < length(jobs)){
            jobs <- jobs[(jobLength+1):length(jobs)]
        }
        else jobs <- c()

        job <- NULL
        
        if(.Platform$OS.type == "unix") 
        {
            rets <- foreach(job = currentJobs) %dopar% {
                clus = job$clus
                
                set.seed(seeds[[clus]][job$iter])
                perturbedRet <- perturbedFunction(data = data)
                
                cMatrix <- BuildConnectivityMatrix(data = perturbedRet$data, clus, clusteringAlgorithm)
                connectivityMatrix = perturbedRet$ConnectivityMatrixHandler(connectivityMatrix = cMatrix$matrix, iter = job$iter, k = clus)
                
                list(
                    connectivityMatrix = connectivityMatrix,
                    clus = clus,
                    auc = CalcAUC(origS[[clus]], connectivityMatrix)$area
                )
            }
        } else {
            rets <- foreach(job = currentJobs) %do% {
                clus = job$clus
                
                set.seed(seeds[[clus]][job$iter])
                perturbedRet <- perturbedFunction(data = data)
                
                cMatrix <- BuildConnectivityMatrix(data = perturbedRet$data, clus, clusteringAlgorithm)
                connectivityMatrix = perturbedRet$ConnectivityMatrixHandler(connectivityMatrix = cMatrix$matrix, iter = job$iter, k = clus)
                
                list(
                    connectivityMatrix = connectivityMatrix,
                    clus = clus,
                    auc = CalcAUC(origS[[clus]], connectivityMatrix)$area
                )
            }
        }
        
        
        allStop = F
        for(ret in rets){
            if (kProgress[clus] != -1){
                clus = ret$clus
                currentIter[clus] <- currentIter[clus] + 1
                pertS[[clus]][[currentIter[clus]]] <- ret$connectivityMatrix
                
                stop <- stoppingCriteriaHandler(iter = currentIter[clus], k = clus, auc = ret$auc)
                if (stop == 2) {
                    allStop <- T
                }
                else if (stop == 1) {
                    if (showProgress) setTxtProgressBar(pb, getTxtProgressBar(pb) + length(jobs[jobs == clus]))
                    jobs <- jobs[jobs != clus]
                    kProgress[clus] <- -1
                }
            }
            if (showProgress) setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
        }
        if (allStop) break()

        countDone <- -1
        while(countDone != length(which(kProgress[clusRange] == -1))){
            countDone <- length(which(kProgress[clusRange] == -1))

            for (clus in clusRange){
                if (kProgress[clus] != -1){
                    stop <- stoppingCriteriaHandler(k = clus, type = 1)
                    if (stop == 1) {
                        if (showProgress) setTxtProgressBar(pb, getTxtProgressBar(pb) + length(jobs[jobs == clus]))
                        jobs <- jobs[jobs != clus]
                        kProgress[clus] <- -1
                    }
                }
            }
        }
    }
    if(.Platform$OS.type == "unix") doParallel::stopImplicitCluster()

    if (showProgress) {
        setTxtProgressBar(pb, maxJob)
        cat("\n")
    }

   

    pertS
}

CalcAUC <- function(orig, pert) {
    N <- nrow(orig)
    S <- abs(orig - pert)
    diag(S) <- 0
    # added -10^(-5) for visual purposes
    # A <- c(-10^(-5), sort(unique(as.numeric(S))))
    A <- c(-10^(-5), 0, sort(unique(as.numeric(S[S!=0]))))
    if (max(A) < 1)
        A <- c(A, 1)
    B <- NULL
    for (i in 1:length(A)) {
        B[i] <- sum(S <= A[i])/(N * N)
    }
    
    area <- 0
    for (i in (2:length(A))) {
        area <- area + B[i - 1] * (A[i] - A[i - 1])
    }
    
    list(area = area, entry = A, cdf = B)
}

CalcPerturbedDiscrepancy <- function(origS, pertS, clusRange) {
    diff <- NULL
    for (clus in clusRange) {
        diff[clus] <- sum(abs(origS[[clus]] - pertS[[clus]]))
    }
    
    AUC <- NULL
    entries <- list()
    cdfs <- list()
    
    for (clus in clusRange) {
        ret <- CalcAUC(origS[[clus]], pertS[[clus]]);
        
        entries[[clus]] <- ret$entry
        cdfs[[clus]] <- ret$cdf
        AUC[clus] <- ret$area
    }
    
    list(Diff = round(diff, digits = 10), Entry = entries, CDF = cdfs, AUC = round(AUC, digits = 10))
}