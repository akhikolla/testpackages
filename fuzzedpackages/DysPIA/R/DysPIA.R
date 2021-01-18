#' @title DysPIA: Dysregulated Pathway Identification Analysis
#' 
#' @description Runs Dysregulated Pathway Identification Analysis (DysPIA).The package 'DysPIAData' including the background data is needed to be loaded.
#' 
#' @param pathwayDB Name of the pathway database (8 databases:reactome,kegg,biocarta,panther,pathbank,nci,smpdb,pharmgkb). 
#'                  The default value is "kegg".
#' @param stats Named vector of CILP scores for each gene pair. Names should be the same as in pathways.
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm. 
#'              The default value is 10000.
#' @param minSize Minimal size of a gene pair set to test. All pathways below the threshold are excluded. 
#'                The default value is 15.
#' @param maxSize Maximal size of a gene pair set to test. All pathways above the threshold are excluded. 
#'                The default value is 1000.
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param DyspiaParam DysPIA parameter value, all gene pair-level status are raised to the power of `DyspiaParam`
#'                    before calculation of DysPIA enrichment scores.
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @return A table with DysPIA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'  \item pathway -- name of the pathway as in `names(pathway)`;
#'  \item pval -- an enrichment p-value;
#'  \item padj -- a BH-adjusted p-value;
#'  \item DysPS -- enrichment score, same as in Broad DysPIA implementation;
#'  \item NDysPS -- enrichment score normalized to mean enrichment of random samples of the same size;
#'  \item nMoreExtreme` -- a number of times a random gene pair set had a more extreme enrichment score value;
#'  \item size -- size of the pathway after removing gene pairs not present in `names(stats)`;
#'  \item leadingEdge -- vector with indexes of leading edge gene pairs that drive the enrichment.
#' }
#'
#' @export
#' @useDynLib DysPIA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom data.table data.table rbindlist setcolorder :=
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom fastmatch fmatch
#' @importFrom stats na.omit
#' @import DysPIAData
#' @examples
#' data(pathway_list,package="DysPIAData")
#' data(DysGPS_p53)
#' DyspiaRes_p53 <- DysPIA("kegg", DysGPS_p53, nperm = 100, minSize = 20, maxSize = 100)
#' 
DysPIA <- function(pathwayDB="kegg", stats, 
                   nperm=10000, minSize=15, maxSize=1000, 
                   nproc=0, DyspiaParam=1, BPPARAM=NULL) {

    pathwayDB <- tolower(pathwayDB)
    # Error if pathwayDB (pathway database name) is not correct
    {
        if (pathwayDB == "reactome")
            pathways <- pathway_list[[1]]
        else if (pathwayDB == "kegg")
            pathways <- pathway_list[[2]]
        else if (pathwayDB == "biocarta")
            pathways <- pathway_list[[3]]
        else if (pathwayDB == "panther")
            pathways <- pathway_list[[4]]
        else if (pathwayDB == "pathbank")
            pathways <- pathway_list[[5]]
        else if (pathwayDB == "nci")
            pathways <- pathway_list[[6]]
        else if (pathwayDB == "smpdb")
            pathways <- pathway_list[[7]]
        else if (pathwayDB == "pharmgkb")
            pathways <- pathway_list[[8]]
        else
            stop("The name of the pathway database you input is not correct!")
    }
    
    # Error if stats is not named
    if (is.null(names(stats))) {
        stop("stats should be named")
    }
    
    # Warning message for duplicate gene pair names
    if (any(duplicated(names(stats)))) {
        warning("There are duplicate gene pair names, DysPIA may produce unexpected results")
    }

    # Getting rid of check NOTEs
    leEs=leZero=geEs=geZero=leZeroSum=geZeroSum=NULL
    pathway=padj=pval=ES=NES=geZeroMean=leZeroMean=NULL
    nMoreExtreme=nGeEs=nLeEs=size=nLeZero=nGeZero=NULL
    leadingEdge=NULL
    
    
    granularity <- 1000
    permPerProc <- rep(granularity, floor(nperm / granularity))
    if (nperm - sum(permPerProc) > 0) {
        permPerProc <- c(permPerProc, nperm - sum(permPerProc))
    }
    seeds <- sample.int(10^9, length(permPerProc))

    BPPARAM <- setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)

    minSize <- max(minSize, 1)
    stats <- sort(stats, decreasing=TRUE)

    stats <- abs(stats) ^ DyspiaParam
    pathwaysFiltered <- lapply(pathways, function(p) { as.vector(na.omit(fmatch(p, names(stats)))) })
    pathwaysSizes <- sapply(pathwaysFiltered, length)

    toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
    m <- length(toKeep)

    if (m == 0) {
        return(data.table(pathway=character(),
                          pval=numeric(),
                          padj=numeric(),
                          DysPS=numeric(),
                          NDysPS=numeric(),
                          nMoreExtreme=numeric(),
                          size=integer(),
                          leadingEdge=list()))
    }

    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]

    DyspiaStatRes <- do.call(rbind,
                lapply(pathwaysFiltered, calcDyspiaStat,
                       stats=stats,
                       returnLeadingEdge=TRUE))

    leadingEdges <- mapply("[", list(names(stats)), DyspiaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(DyspiaStatRes[, "res"])


    pvals <- DyspiaSimpleImpl(pathwayScores, pathwaysSizes, pathwaysFiltered,
                             leadingEdges, permPerProc, seeds, m, stats, BPPARAM)
    if (nrow(pvals[is.na(pval)]) > 0){
        warning("There were ",
                paste(nrow(pvals[is.na(pval)])),
                " pathways for which P-values were not calculated properly due to ",
                "unbalanced gene pair-level statistic values")
    }

    pvals[, nLeZero := NULL]
    pvals[, nGeZero := NULL]
    pvals[, leZeroMean := NULL]
    pvals[, geZeroMean := NULL]
    pvals[, nLeEs := NULL]
    pvals[, nGeEs := NULL]

    setcolorder(pvals, c("pathway", "pval", "padj", "DysPS", "NDysPS",
                         "nMoreExtreme", "size", "leadingEdge"))
    # Makes pvals object printable immediatly
    pvals <- pvals[]
    pvals
}

#' @title calcDyspiaStat: Calculates DysPIA statistics
#' @description Calculates DysPIA statistics for a given query gene pair set.
#'
#' @param stats Named numeric vector with gene pair-level statistics sorted in decreasing order (order is not checked).
#' @param selectedStats Indexes of selected gene pairs in the `stats` array.
#' @param DyspiaParam DysPIA weight parameter (0 is unweighted, suggested value is 1).
#' @param returnAllExtremes If TRUE return not only the most extreme point, but all of them. Can be used for enrichment plot.
#' @param returnLeadingEdge If TRUE return also leading edge gene pairs.
#' @return Value of DysPIA statistic if both returnAllExtremes and returnLeadingEdge are FALSE.
#' Otherwise returns list with the folowing elements:
#' \itemize{
#' \item res -- value of DysPIA statistic
#' \item tops -- vector of top peak values of cumulative enrichment statistic for each gene pair;
#' \item bottoms -- vector of bottom peak values of cumulative enrichment statistic for each gene pair;
#' \item leadingEdge -- vector with indexes of leading edge gene pairs that drive the enrichment.
#' }
#' @export
#' 
calcDyspiaStat <- function(stats, selectedStats, DyspiaParam=1,
                           returnAllExtremes=FALSE,
                           returnLeadingEdge=FALSE) {
    
    S <- selectedStats
    r <- stats
    p <- DyspiaParam
    
    S <- sort(S)
    
    m <- length(S)
    N <- length(r)
    if (m == N) {
        stop("DysPS statistic is not defined when all gene pairs are selected.")
    }
    NR <- (sum(abs(r[S])^p))
    rAdj <- abs(r[S])^p
    if (NR == 0) {
        # this is equivalent to rAdj being rep(eps, m)
        rCumSum <- seq_along(rAdj) / length(rAdj)
    } else {
        rCumSum <- cumsum(rAdj) / NR
    }
    
    
    tops <- rCumSum - (S - seq_along(S)) / (N - m)
    if (NR == 0) {
        # this is equivalent to rAdj being rep(eps, m)
        bottoms <- tops - 1 / m
    } else {
        bottoms <- tops - rAdj / NR
    }
    
    maxP <- max(tops)
    minP <- min(bottoms)
    
    if(maxP > -minP) {
        genePairSetStatistic <- maxP
    } else if (maxP < -minP) {
        genePairSetStatistic <- minP
    } else {
        genePairSetStatistic <- 0
    }
    
    if (!returnAllExtremes && !returnLeadingEdge) {
        return(genePairSetStatistic)
    }
    
    res <- list(res=genePairSetStatistic)
    if (returnAllExtremes) {
        res <- c(res, list(tops=tops, bottoms=bottoms))
    }
    if (returnLeadingEdge) {
        leadingEdge <- if (maxP > -minP) {
            S[seq_along(S) <= which.max(bottoms)]
        } else if (maxP < -minP) {
            rev(S[seq_along(S) >= which.min(bottoms)])
        } else {
            NULL
        }
        
        res <- c(res, list(leadingEdge=leadingEdge))
    }
    res
}

#' @title DyspiaSimpleImpl
#' @description Runs dysregulated pathway identification analysis for preprocessed input data.
#'
#' @importFrom stats p.adjust
#' @param pathwayScores Vector with enrichment scores for the pathways in the database.
#' @param pathwaysSizes Vector of pathway sizes.
#' @param pathwaysFiltered Filtered pathways.
#' @param leadingEdges Leading edge gene pairs.
#' @param permPerProc  Parallelization parameter for permutations.
#' @param seeds Seed vector
#' @param toKeepLength  Number of `pathways` that meet the condition for `minSize` and `maxSize`.
#' @param stats Named vector of gene pair-level scores. Names should be the same as in pathways of `pathwayDB`.
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @return A table with DysPIA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'  \item pathway -- name of the pathway as in `names(pathway)`;
#'  \item pval -- an enrichment p-value;
#'  \item padj -- a BH-adjusted p-value;
#'  \item DysPS -- enrichment score, same as in Broad DysPIA implementation;
#'  \item NDysPS -- enrichment score normalized to mean enrichment of random samples of the same size;
#'  \item nMoreExtreme` -- a number of times a random gene pair set had a more extreme enrichment score value;
#'  \item size -- size of the pathway after removing gene pairs not present in `names(stats)`;
#'  \item leadingEdge -- vector with indexes of leading edge gene pairs that drive the enrichment.
#' }
DyspiaSimpleImpl <- function(pathwayScores, pathwaysSizes, pathwaysFiltered,
                            leadingEdges, permPerProc, seeds,
                            toKeepLength, stats, BPPARAM){
    K <- max(pathwaysSizes)
    universe <- seq_along(stats)

    counts <- bplapply(seq_along(permPerProc), function(i) {
        nperm1 <- permPerProc[i]
        leEs <- rep(0, toKeepLength)
        geEs <- rep(0, toKeepLength)
        leZero <- rep(0, toKeepLength)
        geZero <- rep(0, toKeepLength)
        leZeroSum <- rep(0, toKeepLength)
        geZeroSum <- rep(0, toKeepLength)
        if (toKeepLength == 1) {
            for (i in seq_len(nperm1)) {
                randSample <- sample.int(length(universe), K)
                randEsP <- calcDyspiaStat(
                    stats = stats,
                    selectedStats = randSample,
                    DyspiaParam = 1)
                leEs <- leEs + (randEsP <= pathwayScores)
                geEs <- geEs + (randEsP >= pathwayScores)
                leZero <- leZero + (randEsP <= 0)
                geZero <- geZero + (randEsP >= 0)
                leZeroSum <- leZeroSum + pmin(randEsP, 0)
                geZeroSum <- geZeroSum + pmax(randEsP, 0)
            }
        } else {
            aux <- calcDyspiaStatCumulativeBatch(
                stats = stats,
                DyspiaParam = 1,
                pathwayScores = pathwayScores,
                pathwaysSizes = pathwaysSizes,
                iterations = nperm1,
                seed = seeds[i])
            leEs = get("leEs", aux)
            geEs = get("geEs", aux)
            leZero = get("leZero", aux)
            geZero = get("geZero", aux)
            leZeroSum = get("leZeroSum", aux)
            geZeroSum = get("geZeroSum", aux)
        }
        data.table(pathway=seq_len(toKeepLength),
                   leEs=leEs, geEs=geEs,
                   leZero=leZero, geZero=geZero,
                   leZeroSum=leZeroSum, geZeroSum=geZeroSum
                   )
    }, BPPARAM=BPPARAM)

    counts <- rbindlist(counts)

    # Getting rid of check NOTEs
    leEs=leZero=geEs=geZero=leZeroSum=geZeroSum=NULL
    pathway=padj=pval=DysPS=NDysPS=geZeroMean=leZeroMean=NULL
    nMoreExtreme=nGeEs=nLeEs=size=nLeZero=nGeZero=NULL
    leadingEdge=NULL
    .="damn notes"

    pvals <- counts[, list(leZeroMean = sum(leZeroSum) / sum(leZero),
                           geZeroMean = sum(geZeroSum) / sum(geZero),
                           nLeZero = sum(leZero),
                           nGeZero = sum(geZero),
                           nLeEs = sum(leEs),
                           nGeEs = sum(geEs)),
                    by = .(pathway)]

    pvals[, DysPS := pathwayScores[pathway]]

    pvals[, NDysPS := as.numeric(NA)]
    pvals[(DysPS > 0 & geZeroMean != 0) | (DysPS <= 0 & leZeroMean != 0),
          NDysPS := DysPS / ifelse(DysPS > 0, geZeroMean, abs(leZeroMean))]

    pvals[, pval := as.numeric(NA)]
    pvals[!is.na(NDysPS), pval := pmin((1+nLeEs) / (1 + nLeZero),
                        (1+nGeEs) / (1 + nGeZero))]

    pvals[, padj := as.numeric(NA)]
    pvals[!is.na(pval), padj := p.adjust(pval, method = "BH")]

    pvals[, nMoreExtreme :=  ifelse(DysPS > 0, nGeEs, nLeEs)]
    pvals[, size := pathwaysSizes[pathway]]
    pvals[, pathway := names(pathwaysFiltered)[pathway]]
    pvals[, leadingEdge := .(leadingEdges)]

    pvals
}

#' @title setUpBPPARAM
#' @description Sets up parameter BPPARAM value.
#'
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @importFrom BiocParallel SnowParam MulticoreParam
#' @return parameter BPPARAM value
#' 
setUpBPPARAM <- function(nproc=0, BPPARAM=NULL){
    if (is.null(BPPARAM)) {
        if (nproc != 0) {
            if (.Platform$OS.type == "windows") {
                # windows doesn't support multicore, using snow instead
                result <- SnowParam(workers = nproc)
            } else {
                result <- MulticoreParam(workers = nproc)
            }
        } else {
            result <- bpparam()
        }
        return(result)
    }
    else {
        return(BPPARAM)
    }
}
