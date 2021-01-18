globalVariables(c("nd", "hazard", "Survival")) ## global variables for plot.rocTree

#' Plotting an \code{rocTree} object
#'
#' Plots an \code{rocTree} object. The function returns a \code{dgr_graph} object that is rendered in the RStudio Viewer or survival/hazard estimate at terminal nodes.
#'   
#' @param x an object of class "\code{rocTree}", usually returned by the rocTree function.
#' @param output a string specifying the output type; graph (the default) renders the graph using the \code{grViz} function, and \code{visNetwork} renders the graph using the visnetwork function.
#' @param digits the number of digits to print.
#' @param rankdir is a character string specifying the direction of the tree flow. The available options are top-to-bottom ("TB"), bottom-to-top ("BT"), left-to-right ("LR"),
#' and right-to-left ("RL"); the default value is "TB".
#' @param shape is a character string specifying the shape style.
#' Some of the available options are "ellipse", "oval", "rectangle", "square", "egg", "plaintext", "diamond", and "triangle". The default value is "ellipse".
#' @param nodeOnly is a logical value indicating whether to display only the node number; the default value is "TRUE".
#' @param savePlot is a logical value indicating whether the plot will be saved (exported); the default value is "FALSE".
#' @param file_name is a character string specifying the name of the plot when "savePlot = TRUE". The file name should include its extension. The default value is "pic.pdf"
#' @param file_type is a character string specifying the type of file to be exported. Options for graph files are: "png", "pdf", "svg", and "ps". The default value is "pdf".
#' @param type is an optional character string specifying the type of plots to produce. The available options are "tree" for plotting survival tree (default),
#' "survival" for plotting the estimated survival probabilities for the terminal nodes, and "hazard" for plotting the estimated hazard for the terminal nodes.
#' The \code{dgr_graph} options are available only when \code{type = tree}.
#' @param tree is an optional integer specifying the \eqn{n^{th}} tree in the forest to print.
#' @param ... arguments to be passed to or from other methods.
#'
#' @seealso See \code{\link{rocTree}} for creating \code{rocTree} objects.
#' @importFrom DiagrammeR export_graph render_graph
#' @importFrom data.tree SetGraphStyle SetNodeStyle
#' @importFrom ggplot2 labs geom_smooth
#' 
#' @export
#' @example inst/examples/ex_plot_rocTree.R
plot.rocTree <- function(x, output = c("graph", "visNetwork"),
                         digits = 4, tree = 1L,
                         rankdir = c("TB", "BT", "LR", "RL"),
                         shape = "ellipse",
                         nodeOnly = FALSE,
                         savePlot = FALSE, 
                         file_name = "pic.pdf",
                         file_type = "pdf",
                         type = c("tree", "survival", "hazard"),...) {
    if (!is.rocTree(x)) stop("Response must be a \"rocTree\" object")
    output <- match.arg(output)
    type <- match.arg(type)
    rankdir <- match.arg(rankdir)
    if (x$ensemble) {
        if (!is.wholenumber(tree)) stop("Tree number must be an integer.")
        if (tree > length(x$trees)) stop("Tree number exceeded the number of trees in forest.")    
        Frame <- x$Frame[[tree]]
    } else {
        Frame <- x$Frame
        if (type == "survival") {
            tmp <- data.frame(x$data$.Y0, x$data$.X0)[x$data$.D0 > 0,]
            colnames(tmp) <- c(x$rName, x$vNames)
            t0 <- tmp[,1]
            atTerm <- lapply(split(tmp, x$nodeLabel), function(xx) {
                xx <- xx[findInt(t0, xx[,1]),]
                xx[,1] <- t0
                rownames(xx) <- NULL
                return(xx)
            })
            atTerm <- lapply(atTerm, function(xx)
                data.frame(Time = t0, Survival = predict(x, newdata = xx)$survFun(t0)))
            atTerm <- do.call(rbind, atTerm)
            atTerm$nd <- as.factor(rep(sort(unique(x$nodeLabel)) + 1, each = length(t0)))
            ## atTerm$nd <- as.factor(rep(sort(unique(x$nodeLabel)) + 1, table(x$nodeLabel)))
            rownames(atTerm) <- NULL
            gg <- ggplot(atTerm, aes(x = Time, y = Survival, col = nd)) + geom_step(lwd = I(1.1)) +
                xlab("Time") + ylab("Survival probabilities") + labs(col = "Node")
            return(gg)
        }
        if (type == "hazard") {
            tmp <- data.frame(x$data$.Y0, x$data$.X0)[x$data$.D0 > 0,]
            colnames(tmp) <- c(x$rName, x$vNames)
            t0 <- tmp[,1]
            atTerm <- lapply(split(tmp, x$nodeLabel), function(xx) {
                xx <- xx[findInt(t0, xx[,1]),]
                xx[,1] <- t0
                rownames(xx) <- NULL
                return(xx)
            })
            atTerm <- lapply(atTerm, function(xx)
                data.frame(Time = t0, hazard = predict(x, newdata = xx, type = "haz")$hazFun(t0)))
            atTerm <- do.call(rbind, atTerm)
            atTerm$nd <- as.factor(rep(sort(unique(x$nodeLabel)) + 1, each = length(t0)))
            rownames(atTerm) <- NULL
            gg <- ggplot(atTerm, aes(x = Time, y = hazard, col = nd)) +
                geom_smooth(method = "loess", se = FALSE) + 
                xlab("Time") + ylab("Hazard estimates") + labs(col = "Nodes")
            return(gg)
        }
    }
    if (type == "tree") {    
        ## create data.tree
        root <- Node$new("Root", type = "root", decision = "", nd = 1)
        for (i in 2:nrow(Frame)) {
            if (i <= 3) parent <- "root"
            if (i > 3) parent <- paste0("Node", Frame$nd[i] %/% 2)
            if (Frame$is.terminal[i] == 2) {
                type <- "terminal"
                display <- with(
                    Frame, paste0(nd[i], ") ", tree.split.names(nd[i], nd, p, cutVal, x$vNames, digits), "*"))
            } else {
                type <- "interior"
                display <- with(
                    Frame, paste0(nd[i], ") ", tree.split.names(nd[i], nd, p, cutVal, x$vNames, digits)))
            }
            eval(parse(text = paste0("Node", Frame$nd[i], "<-", parent,
                                     "$AddChild(display, type = type, nd = Frame$nd[i])")))
        }
        SetGraphStyle(root, rankdir = rankdir)
        if (nodeOnly) {
            GetNodeLabel <- function(node) {
                switch(node$type,
                       root = "Root", 
                       terminal = paste0("Node ", node$nd),
                       interior = paste0("Node ", node$nd))
            }
            SetNodeStyle(root, fontname = 'helvetica', label = GetNodeLabel, shape = shape)
        } else {
            SetNodeStyle(root, fontname = 'helvetica', shape = shape)
        }
        ## x$graph <- ToDiagrammeRGraph(x, direction = "climb", pruneFun = NULL)
        plot.Node(root, output = output, control = list(savePlot = savePlot,
                                                        file_name = file_name, file_type = file_type))
    }
}

## ---------------------------------------------------------------------------------------
## From data.tree package that are not exported by them
## ---------------------------------------------------------------------------------------

#' @importFrom data.tree ToDiagrammeRGraph
#' @keywords internal
#' @noRd
plot.Node <- function(x, ..., output = "graph", control) {
    graph <- ToDiagrammeRGraph(x, direction = "climb", pruneFun = NULL)
    if (control$savePlot) {
        export_graph(graph, file_name = control$file_name, file_type = control$file_type)
        render_graph(graph, output = output)
    } else {
        render_graph(graph, output = output)
    }
}
