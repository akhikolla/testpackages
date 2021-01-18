.termNode <- function(tree, term = list()) {
    if(tree$Type == 'Terminal') {
        term[[paste0('term', tree$ID)]] = tree
        return(term)
    } else {
        term = .termNode(tree$Left, term)
        term = .termNode(tree$Right, term)
    }
    return(term)
}

## Missing data support
#' @importFrom stats as.formula lm
.node <- function(termList, data, ynames, trtname) {
    model <- list()
    yp = length(ynames)
    fitIndex <- termList$FitIndex
    for (i in 1:yp) {
        yname = ynames[i]
        fitvar = fitIndex[[i]]
        if (length(fitvar) == 0) {
            reg_for = paste0(yname, '~ -1 + ', trtname)
        } else {
            reg_for = paste0(yname, '~ -1 + ', paste(fitvar, collapse = '+'), ' + ', trtname)
        }
        model[[yname]] <- lm(as.formula(reg_for), data = data)
    }
    return(model)
}

.node.guide <- function(treeRes, nodeID, dataframe, ynames, trtname) {
    uniID <- unique(nodeID)
    termNodeMap = .termNode(treeRes)
    for (id in uniID) {
        naID = paste0('term', id)
        ind = nodeID == id
        treeTmp = termNodeMap[[naID]]
        stopifnot(treeTmp$Type == 'Terminal')
        termNodeMap[[naID]][['model']] = .node(treeTmp, dataframe[ind, ], ynames, trtname)
    }
    return(termNodeMap)
}

## tmp = node.guide(res_subguide$treeRes, res_subguide$node, train_data, yname, "Z")
