.predictNode <- function(node, dataframe, cLevels) {
    if (node$Type == 'Terminal') {
        return(node$ID)
    } else {
        val = dataframe[[node$SplitVar]]
        if (node$Role == 'num') {
            if (is.na(val)) {
                if (node$MisDirection %in% c('A', 'L')) {
                    return(.predictNode(node$Left, dataframe, cLevels))
                } else {
                    return(.predictNode(node$Right, dataframe, cLevels))
                }
            } else {
                if (val <= node$Threshold) {
                    return(.predictNode(node$Left, dataframe, cLevels))
                } else {
                    return(.predictNode(node$Right, dataframe, cLevels))
                }
            }
        } else {
            if (val %in% node$ThreshSet) {
                return(.predictNode(node$Left, dataframe, cLevels))
            } else {
                return(.predictNode(node$Right, dataframe, cLevels))
            }
        }
    }
}

#' @importFrom stats predict.lm
#' @noRd
.predictY <- function(node, newx) {
    models = node[['model']]
    ny = length(models)

    y = rep(NA, ny)
    for (i in 1:ny) {
        y[i] = predict.lm(models[[i]], newx)
    }
    y
}

.processTrt <- function(nodeMap, yname, trtname, tLevels) {
    termNodes = names(nodeMap)
    ny = length(nodeMap[[termNodes[1]]][['Trts']])
    nt = length(nodeMap[[termNodes[1]]][['Trts']][[1]])
    trt = data.frame()

    for (term in termNodes) {
        tmp = sapply(nodeMap[[term]][['Trts']], FUN = function(x){x})
        if (ny == 1) tmp = t(tmp)
        colnames(tmp) = paste0(trtname, tLevels)
        trt = rbind(trt,
                    data.frame(nodeId = gsub('term', '', term),
                               node = term, y = yname, tmp))
    }
    trt
}

#' Predict the node id of MrSGUIDE regression tree
#'
#' @param mrsobj MrSGUIDE object
#' @param dataframe data used for prediction
#' @param type node id
#'
#' @return A data frame of each object node id and outcome
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
#' mrsobj <- MrSFit(dataframe = train, role = role)
#' newX = train[1:10,]
#' predictTree(mrsobj, newX, type='outcome')
#' predictTree(mrsobj, newX, type='node')
#'
#' @export
predictTree <- function(mrsobj, dataframe, type = 'node') {
    stopifnot(class(mrsobj) == 'guide')
    stopifnot(type %in% c("node", "outcome"))
    n = NROW(dataframe)
    yp = mrsobj$yp
    tp = mrsobj$tp
    if (is.null(dataframe)) return(mrsobj[['node']])
    node <- sapply(1:n, FUN = function(i) {.predictNode(mrsobj$treeRes, dataframe[i, ], mrsobj$cLevels)})
    if (type=='node') return(node)

    nodemap = mrsobj[['nodeMap']]
    tn = paste0('term', node)

    y = sapply(1:n, FUN = function(i) {.predictY(nodemap[[tn[i]]], dataframe[i,])})
    y = t(y)
    ans = data.frame(node = node, y)
    colnames(ans) = c('node', mrsobj$ynames)
    return(ans)
}
