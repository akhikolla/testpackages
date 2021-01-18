.getTrt <- function(termMap, Yname, tLevels) {
    .helper <- function(nodeTrt, nodeSE) {
        res <- c()
        for(i in seq_along(Yname)) {
            trt = nodeTrt[[i]]
            trtse = nodeSE[[i]]
            Effect = c()
            SE = c()
            for (j in seq_along(trt)) {
                Effect = c(Effect, trt[j])
                SE = c(SE, trtse[j])
            }
            resT <- data.frame(Estimate = Effect, SE = SE, Assignment = tLevels, Outcome = Yname[i])
            res <- rbind(res, resT[-1, ])
        }
        return(res)
    }
    nodeTrt = lapply(termMap, FUN = function(x) {x$Trts})
    nodeSE = lapply(termMap, FUN = function(x) {x$SEs})
    nodeName = names(termMap)
    trtData <- c()
    for (i in seq_along(nodeName)) {
        node = nodeName[i]
        tmp = .helper(nodeTrt[[i]], nodeSE[[i]])
        tmp$Node = gsub('term', '', node)
        trtData = rbind(trtData, tmp)
    }
    trtData
}

.node_df <- function(node, treatNode, level = 1,
                     digits = 3, font.size = 16,
                     nodesPopSize = FALSE, fixed = FALSE,
                     minNodeSize = 15,
                     maxNodeSize = 30) {

    if (nodesPopSize) {
        minNodeSize = minNodeSize
        maxNodeSize = maxNodeSize
    } else {
        minNodeSize = (minNodeSize + maxNodeSize) / 2
        maxNodeSize = minNodeSize
    }
    settings <- data.frame(font.size = font.size,
                           font.color="black",
                           fixed = fixed, font.multi = TRUE,
                           scaling.min = minNodeSize, scaling.max = maxNodeSize)
    .node_helper <- function(node, treatNode, level, digits) {
        if (node$Type == 'Terminal') {
            trt = treatNode[treatNode$Node == node$ID, 'Estimate']
            fillcolor = ifelse(mean(trt) > 0, "#FC4E07", "#00AFBB")
            ndf <- data.frame(id = node$ID, level = level,
                              type = node$Type,
                              label = paste0("Node ", node$ID, ": ", node$Size),
                              color = fillcolor,
                              shape = 'square')
            return(ndf)
        } else {
            text = paste("Node ", node$ID, ":\n", node$SplitVar, sep = '')
            ndf <- data.frame(id = node$ID, level = level,
                              type = node$Type,
                              label = text,
                              color = "lightblue",
                              shape = 'dot')
            ndf <- rbind(ndf,
                         .node_helper(node$Left, treatNode, level = level + 1,
                                      digits = digits),
                         .node_helper(node$Right, treatNode, level = level + 1,
                                      digits = digits))
        }
    }
    ndf <- .node_helper(node, treatNode, level)
    return(cbind(ndf, settings))
}

.edge_df <- function(node, node_df, clevels, digits = 3,
                     color = "#8181F7", font.size = 16, font.align = "horizontal") {
    settings <- data.frame(color = color, font.size = font.size, font.align = font.align,
                           smooth.enabled = TRUE, smooth.type = "cubicBezier",
                           smooth.roundess = 0.5)

    .edge_helper <- function(node, node_df, clevels, digits) {
        if (node$Type == 'Terminal') {
            return()
        } else {
            id = node$ID
            lid = id * 2
            rid = id * 2 + 1
            if (node$Role == 'num') {
                if (node$MisDirection != 'A') {
                    ldis = paste0('<=', ifelse(node$MisDirection == 'L', '* ', ' '), round(node$Threshold, digits))
                    rdis = paste0('>', ifelse(node$MisDirection == 'L', ' ', '* '), round(node$Threshold, digits))
                } else {
                    ldis = paste0(' = NA')
                    rdis = paste0(' != NA')
                }
            } else {
                if (node$MisDirection != 'A') {
                    ldis = paste0("{ ", paste0(node$ThreshSet, collapse = ', '),
                                  ifelse(node$MisDirection == 'L', ', NA }', ' }'))
                    varLevel = clevels[[node$SplitVar]]
                    rdis = paste0("{ ", paste0(varLevel[which(!varLevel %in% node$ThreshSet)], collapse = ', '),
                                  ifelse(node$MisDirection == 'L', ' }', ', NA }'))
                } else {
                    ldis = paste0(' = NA')
                    rdis = paste0(' != NA')
                }
            }
            edf <- data.frame(from = c(id, id), to = c(lid, rid),
                              label = c(ldis, rdis))
            edf <- rbind(edf,
                         .edge_helper(node$Left, node_df, clevels, digits),
                         .edge_helper(node$Right, node_df, clevels, digits))
        }
    }

    edf <- .edge_helper(node, node_df, clevels, digits)
    if (is.null(edf)) return(NULL)

    return(cbind(edf, settings))
}


#' Plot MrSGUIDE regression tree
#'
#' @param mrsobj MrSGUIDE object
#' @param digits digits for split threshold
#' @param height figure height
#' @param width figure width
#' @param nodefontSize node font size
#' @param edgefontSize edge font size
#' @param minNodeSize minimal node size
#' @param maxNodeSize maximum node size
#' @param nodeFixed whether you can drag node
#' @param edgeColor edge color
#' @param highlightNearest choose node will highlight nearby
#' @param collapse list, collapse or not using double click on a node
#' @param alphaInd 1 is original alpha, 2 is individual level alpha, 3 is overall alpha
#'
#' @return A list contains plot figure
#' \item{treeplot}{The tree plot uses \code{\link[visNetwork]{visNetwork}} function.}
#' \item{nodeTreat}{A data frame contain each elements used for tree plot.}
#' \item{trtPlot}{A treatment effects plot of each node.}
#'
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
#' plotObj <- plotTree(mrsobj)
#' #plotObj$treePlot
#' plotObj$nodeTreat ## node information
#' plotObj$trtPlot ## treatment effect plot
#'
#' @rdname plotTree
#' @export
#'
plotTree <- function(mrsobj, digits = 3, height = "600px", width = "100%",
                       nodefontSize = 16, edgefontSize = 14,
                       minNodeSize = 15, maxNodeSize = 30,
                       nodeFixed = FALSE, edgeColor = "#8181F7",
                       highlightNearest =  list(enabled = TRUE,
                                                degree = list(from = 50000, to = 0), hover = FALSE,
                                                algorithm = "hierarchical"),
                       collapse = list(enabled = FALSE, fit = TRUE, resetHighlight = TRUE,
                                       clusterOptions = list(fixed = TRUE, physics = FALSE)),
                       alphaInd = 3) {

    treatNode <- .getTrt(mrsobj$nodeMap, mrsobj$ynames, paste0(mrsobj$trtname, '.', mrsobj$tLevels[[1]]))
    ndf1 <- .node_df(mrsobj$treeRes, treatNode, digits = digits, font.size = nodefontSize,
                     nodesPopSize = FALSE, fixed = nodeFixed,
                     minNodeSize = minNodeSize,
                     maxNodeSize = maxNodeSize)

    edf1 <- .edge_df(mrsobj$treeRes, ndf1, mrsobj$cLevels, digits = digits, color = edgeColor,
                     font.size = edgefontSize, font.align = "horizontal")

    palpha <- 1.96
    if (!is.null(mrsobj$bootAlpha)) {
        palpha <- abs(stats::qnorm(mrsobj$bootAlpha[alphaInd]))
    }

    treatNode$Quantity <- paste0('Node: ', treatNode$Node,'; ', treatNode$Outcome)
    treatNode$ymin <- treatNode$Estimate - palpha * treatNode$SE
    treatNode$ymax <- treatNode$Estimate + palpha * treatNode$SE
    treatNode$zalpha <- palpha

    if (requireNamespace("visNetwork", quietly = TRUE)) {
        tree <- visNetwork::visNetwork(nodes = ndf1, edges = edf1, height = height, width = width) %>%
            visNetwork::visHierarchicalLayout(direction = 'UD') %>%
            visNetwork::visPhysics(barnesHut = list(avoidOverlap = 1)) %>%
            # visOptions(highlightNearest =  highlightNearest, collapse = collapse) %>%
            visNetwork::visInteraction(dragNodes = !nodeFixed, selectConnectedEdges = FALSE,
                                       tooltipStyle = 'position: fixed;visibility:hidden;padding: 5px;
                                   white-space: nowrap;
                                   font-family: cursive;font-size:12px;font-color:purple;background-color: #E6E6E6;
                                   border-radius: 15px;') %>%
            visNetwork::visEdges(scaling = list(label = list(enabled = FALSE))) %>%
            visNetwork::visEvents(type = "once", stabilized = "function() {
                              this.setOptions({layout:{hierarchical:false}, physics:{solver:'barnesHut', enabled:true, stabilization : false}, nodes : {physics : false, fixed : true}});
}") %>%
            visNetwork::visExport()
    } else {
        warning("Package \"visNetwork\" needed for this function to work better. Please install it.")
        tree <- NULL
    }

    if (requireNamespace("ggplot2", quietly = TRUE)) {
        trtPlot <- ggplot2::ggplot(data = treatNode,
                                   ggplot2::aes_string(x = 'Quantity', y = 'Estimate',
                                              color = 'Node', group = 'Node',
                                              shape = 'Assignment')) +
            ggplot2::xlim(rev(treatNode$Quantity)) +
            ggplot2::geom_point(size=6) +
            ggplot2::labs(y = 'Treatment Effect') +
            ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "red") +
            ggplot2::coord_flip() +
            ggplot2::theme_bw(base_size = 20)

        if (!is.null(mrsobj$bootAlpha)) {
            trtPlot <- trtPlot +
                ggplot2::geom_errorbar(ggplot2::aes_string(x = "Quantity",
                                                           ymin = "ymin",
                                                           ymax = "ymax"), data = treatNode, width=0.1, size=1)
        }
    } else {
        warning("Package \"ggplot2\" needed to plot treatment effect with node. Please install it.")
        trtPlot <- NULL
    }

    list(treeplot = tree, nodeTreat = treatNode, trtPlot = trtPlot)
}

# plot.guideImp <- function(mrsobj) {
#
# }
