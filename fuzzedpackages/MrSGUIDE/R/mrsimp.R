#' MrSGUIDE variable importance
#'
#' @title Variable importance for predictive variables
#' @author Peigen Zhou
#'
#' @name MrSImp
#'
#' @description Variable importance in subgroup identification for predictive variables.
#'
#' @param dataframe train data frame
#' @param role role follows GUIDE role
#' @param B bootstrap number default = 100
#' @param bestK number of covariates in the regression model
#' @param maxDepth maximum tree depth
#' @param minTrt minimum treatment and placebo sample in each node
#' @param minData minimum sample in each node
#' @param batchNum related with exhaustive search for numerical split variable
#' @param faster related with tree split searching
#' @param display Whether display tree in the end
#' @param treeName yaml file for save the tree
#' @param nodeName file same for each node
#' @param impName important variable file name
#'
#'
#' @return A list contains importance score variable names and roles
#' \item{imp}{Importance score data frame}
#' \item{role}{Role for each variable}
#' \item{Settings}{Settings used to build the tree}
#'
#'
#' @importFrom utils read.table write.table
#'
#' @examples
#' library(MrSGUIDE)
#' set.seed(1234)
#'
#' N = 200
#' np = 3
#'
#' numX <- matrix(rnorm(N * np), N, np) ## numerical features
#' gender <- sample(c('Male', 'Female'), N, replace = TRUE)
#' country <- sample(c('US', 'UK', 'China', 'Japan'), N, replace = TRUE)
#'
#' z <- sample(c(0, 1), N, replace = TRUE) # Binary treatment assignment
#'
#' y1 <- numX[, 1] + 1 * z * (gender == 'Female') + rnorm(N)
#' y2 <- numX[, 2] + 2 * z * (gender == 'Female') + rnorm(N)
#'
#' train <- data.frame(numX, gender, country, z, y1, y2)
#' role <- c(rep('n', 3), 'c', 'c', 'r', 'd', 'd')
#'
#' mrsobj <- MrSImp(dataframe = train, role = role, B = 10)
#' mrsobj$imp
#'
#'
#' @export
MrSImp <- function(dataframe, role, B = 100, bestK = 1,
                   maxDepth = 5,
                   minTrt = 5, minData = max(c(minTrt * maxDepth, NROW(Y) / 20)),
                   batchNum = 1L,
                   faster = FALSE, display = FALSE,
                   treeName = paste0("tree_", format(Sys.time(), "%m%d"), ".yaml"),
                   nodeName = paste0("node_", format(Sys.time(), "%m%d"), ".txt"),
                   impName = paste0("imp_", format(Sys.time(), "%m%d"), ".txt")) {
    t1 = Sys.time()

    if(display) cat("Start data processing: \n")
    process_list = .data_process(dataframe, role)
    Y = process_list[['Y']]
    cLevels = process_list[['cLevels']]
    cXL = process_list[['cXL']]
    nX = process_list[['nX']]
    tLevels = process_list[['tLevels']]
    TrtL = process_list[['TrtL']]
    numVarName = process_list[['numVarName']]
    catVarName = process_list[['catVarName']]
    newVar = process_list[['newVar']]
    nX = process_list[['nX']]
    fitIndex = process_list[['fitIndex']]
    splitIndex = process_list[['splitIndex']]
    holdIndex = process_list[['holdIndex']]
    non_miss = process_list[['non_miss']]
    varName = process_list[['varName']]


    if(display) cat("Finish processing, start call main function. ", Sys.time() - t1, "s\n")
    t2 = Sys.time()
    stopifnot(NROW(Y) > 0)
    CVFolds = 0
    CVSE = 0
    bootNum = 0
    alpha = 0.05
    treeRes = .mrs.pure(nX = nX, cX = cXL$intX,
                       Y = Y, Trt = TrtL$intX,
                       splitIndex = splitIndex, fitIndex = fitIndex,
                       holdIndex = holdIndex, bk = bestK, maxDepth = maxDepth,
                       minTrt = minTrt, minData = minData, batchNum = batchNum,
                       CVFolds = CVFolds, CVSE = CVSE, bootNum = bootNum,
                       alpha = alpha, faster = faster, display = FALSE, varName = newVar,
                       treeName = treeName, nodeName = nodeName,
                       impName = impName)
    #impRes <- cbind(read.table(impName), newVar)
    impO <- read.table(impName)
    #colnames(impRes) = c('Importance', 'Feature')
    file.remove(nodeName, treeName, impName)

    imp_baseline <- matrix(0, length(newVar), B)
    N = NROW(Y)

    for (i in 1:B) {
        idx = sample(1:N, N)
        #if (NCOL(Y) > 1) Yb = Y[idx,]
        #else Yb = Y[idx]
        Yb = as.matrix(Y[idx,])
        fitB <- .mrs.pure(nX = nX, cX = cXL$intX,
                         Y = Yb, Trt = TrtL$intX,
                         splitIndex = splitIndex, fitIndex = fitIndex,
                         holdIndex = holdIndex, bk = bestK, maxDepth = maxDepth,
                         minTrt = minTrt, minData = minData, batchNum = batchNum,
                         CVFolds = CVFolds, CVSE = CVSE, bootNum = bootNum,
                         alpha = alpha, faster = faster, display = FALSE, varName = newVar,
                         treeName = treeName, nodeName = nodeName,
                         impName = impName)
        imp_baseline[,i] = read.table(impName)$V1
        file.remove(nodeName, treeName, impName)
    }

    imp = pmax(impO$V1 - rowMeans(imp_baseline, na.rm=TRUE), 0)
    impRes <- data.frame(Importance = imp, Feature = newVar)
    impRes <- impRes[order(impRes[,1]), ]

    Settings <- list(CVFolds = 0, CVSE = 0, bestK = bestK, maxDepth = maxDepth,
                     minData = minData, minTrt = minTrt)

    res <- list(imp = impRes,
                role = role, varName = varName,
                Settings = Settings)
    if(display) {
        cat("Finish importance. ", difftime(Sys.time(), t2, units = "secs"), "s\n")
    }

    class(res) <- c('guideImp')
    return(res)
}
