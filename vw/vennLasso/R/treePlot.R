



#' plotting function to investigate hierarchical structure of selection
#'
#' @param object fitted vennLasso object
#' @param s lambda value for the predictions. Only one can be specified at a time
#' @param type type of plot to make. Currently only "d3.tree" and "igraph.tree" available
#' @param ... other graphical parameters for the plot
#' @importFrom graphics plot
#' @importFrom visNetwork visNetwork visIgraphLayout visLegend visOptions visInteraction visNodes
#' @import igraph
#' @export
#' @examples
#' set.seed(123)
#'
#' dat.sim <- genHierSparseData(ncats = 3, nvars = 25, nobs = 200)
#'
#' fit <- vennLasso(x = dat.sim$x, y = dat.sim$y, groups = dat.sim$group.ind)
#'
#' plotSelections(fit, s = fit$lambda[32])
#'
#'
plotSelections <- function(object, s = NULL,
                           type = c("d3.tree"),
                           ...)
{

    if (class(object)[1] == "cv.vennLasso")
    {
        if (is.null(s))
        {
            s = object$lambda.min
        } else if (is.character(s))
        {
            if (s == "lambda.min")
            {
                s <- object$lambda.min
            } else if (s == "lambda.1se")
            {
                s <- object$lambda.1se
            } else stop("Invalid specification for s (lambda)")
        }
        object <- object$vennLasso.fit
    }

    type         <- match.arg(type)

    nbeta        <- object$beta
    combin.mat   <- object$condition.combinations
    combin.names <- object$combin.names
    n.conditions <- object$n.conditions
    M            <- object$n.combinations
    N            <- object$nobs
    p            <- object$nvars

    dimnames(nbeta) <- list(NULL, NULL, NULL)
    if(!is.null(s))
    {
        if (length(s) > 1)
        {
            s <- s[1]
            warning("multiple s values given, using first value only")
        }
        
        if (s < 0) stop("lambda must be positive")
        
        lambda  <- object$lambda
        lamlist <- lambdaInterp(lambda, s)
        nbeta   <- nbeta[,,lamlist$left,drop=TRUE]  * lamlist$frac + 
                   nbeta[,,lamlist$right,drop=TRUE] * (1 - lamlist$frac)
        rownames(nbeta) <- combin.names
        colnames(nbeta) <- object$var.names
        nlam <- length(s)
    } else 
    {
        stop("must specify s, the lambda value for which to plot estimates")
    }

    direct.above.idx.list <- above.idx.list <- vector(mode = "list", length = M)

    for (c in 1:M) 
    {
        indi <- which(combin.mat[c,] == 1)

        # get all indices of g terms which are directly above the current var
        # in the hierarchy, ie. for three categories A, B, C, for the A group,
        # this will return the indices for AB and AC. For AB, it will return
        # the index for ABC. for none, it will return the indices for
        # A, B, and C, etc
        inner.loop <- (1:(M))[-c]
        for (j in inner.loop) 
        {
            diffs.tmp <- combin.mat[j,] - combin.mat[c,]
            if (all( diffs.tmp >= 0 )) 
            {
                if (sum( diffs.tmp == 1 ) == 1) 
                {
                    direct.above.idx.list[[c]] <- c(direct.above.idx.list[[c]], j)
                }
                above.idx.list[[c]] <- c(above.idx.list[[c]], j)
            }
        }
        above.idx.list[[c]] <- c(above.idx.list[[c]], c)
    }
    rsc <- rowSums(combin.mat)
    if (any(rsc == 0)) {
        above.idx.list[[which(rsc == 0)]] <- which(rsc == 0)
    }



    plotMat <- plotMatZero <- array(0, dim = rep(M, 2))
    colnames(plotMat) <- rownames(plotMat) <- combin.names


    edges <- data.frame(array(NA, dim = c(length(unlist(direct.above.idx.list)), 3)))
    colnames(edges) <- c("from", "to", "value")
    nodes <- NULL
    coefs <- varnames <- list()

    e.ct <- 0
    for (c in 1:M)
    {
        num.2.plot <- sum(nbeta[c, ] != 0)
        above.idx.cur <- direct.above.idx.list[[c]]
        for (k in 1:length(direct.above.idx.list))
        {
            if (!is.null(above.idx.cur[k]) & num.2.plot)
            {
                e.ct <- e.ct + 1

                edges[e.ct, ] <- c(above.idx.cur[k], c, num.2.plot)
                if (num.2.plot)
                {
                    nodes <- c(nodes, c, above.idx.cur[k])
                }
            }
        }
    }


    plotMat <- plotMatZero <- array(0, dim = rep(M, 2))
    colnames(plotMat) <- rownames(plotMat) <- combin.names
    for (c in 1:M)
    {
        num.2.plot <- sum(nbeta[c, ] != 0)
        above.idx.cur <- direct.above.idx.list[[c]]
        for (k in 1:length(direct.above.idx.list))
        {
            plotMat[above.idx.cur[k], c] <- num.2.plot
        }
    }


    sums.inds <- sapply(combin.names, function(x) sum(as.numeric(strsplit(x, ",")[[1]])) )
    ord.names.decr <- order(sums.inds, decreasing = TRUE)

    ## all the observed numbers of conditions
    unique.cond.nums <- unique(sums.inds)
    unique.cond.nums <- unique.cond.nums[order(unique.cond.nums, decreasing = TRUE)]

    zero.sm.idx <- which(unique.cond.nums == 0)
    ## the number of times each number of conditions will appear in the graph
    num.nodes.per.row <- sapply(unique.cond.nums, function(x) sum(sums.inds == x))


    zero.idx <- which(sums.inds[ord.names.decr] == 0)
    if (length(zero.idx))
    {
        plotMat.ordered <- plotMat[ord.names.decr[-zero.idx], ord.names.decr[-zero.idx]]
        plotMatZero     <- plotMatZero[ord.names.decr[-zero.idx], ord.names.decr[-zero.idx]]
        names           <- combin.names[ord.names.decr][-zero.idx]
        positions       <- num.nodes.per.row[-zero.sm.idx]
    } else
    {
        plotMat.ordered <- plotMat[ord.names.decr, ord.names.decr]
        names           <- combin.names[ord.names.decr]
        positions       <- num.nodes.per.row
    }

    cs <- colSums(plotMat.ordered)
    rs <- rowSums(plotMat.ordered)
    csplusrs <- cs + rs

    nonzero.keep <- csplusrs != 0

    plotMat.keep <- plotMat.ordered[nonzero.keep, nonzero.keep]


    edges <- na.omit(edges)
    nodes <- unique(nodes)
    nodes <- nodes[!is.na(nodes)]
    nodes.df <- data.frame(id = nodes, label = combin.names[nodes], stringsAsFactors = FALSE)

    num.nz.per.strata  <- apply(nbeta, 1, function(x) sum(x != 0))
    
    # save nonzero coefficients for each subpopulations
    beta.nz.per.strata <- lapply(apply(nbeta, 1, function(x) list(x[x != 0])  ), "[[", 1)

    nodes.df <- nodes.df[nodes.df$label %in% combin.names[num.nz.per.strata > 0], ]
    nodes.df$value <- num.nz.per.strata[match(nodes.df$label, names(num.nz.per.strata))]
    
    # reorder subpopulations to align with how nodes are ordered
    beta.nz.per.strata <- beta.nz.per.strata[match(nodes.df$label, names(num.nz.per.strata))]

    sums.inds <- sapply(combin.names, function(x) sum(as.numeric(strsplit(x, ",")[[1]])) )
    zero.idx <- which(sums.inds == 0)

    if (length(zero.idx))
    {
        edges <- edges[edges$to != zero.idx,]
        edges <- edges[edges$from != zero.idx,]
    }

    #visNetwork(nodes.df, edges)


    node.sizes <- num.nz.per.strata[match(colnames(plotMat.keep), names(num.nz.per.strata))]
    

    plotMat.keep.sqrt <- plotMat.keep
    plotMat.keep.sqrt[plotMat.keep.sqrt != 0] <- sqrt(plotMat.keep.sqrt[plotMat.keep.sqrt != 0])



    edges2    <- edges
    nodes.df2 <- nodes.df
    nodes.df2$group <- 1

    unique.ids <- unique(c(edges$from, edges$to, nodes.df$id))

    new.ids <- 0:(length(unique.ids) - 1)

    edges2$from <- new.ids[match(edges2$from, unique.ids)]
    edges2$to   <- new.ids[match(edges2$to, unique.ids)]

    edges2 <- edges2[order(edges2$from),]
    colnames(edges2)[1:2] <- c("source", "target")


    nodes.df2$id <- new.ids[match(nodes.df2$id, unique.ids)]

    nodes.df2 <- nodes.df2[order(nodes.df2$id),]

    grp.memberships <- data.frame(t(sapply(strsplit(nodes.df2$label, ","), function(x) as.numeric(x) )),
                                  stringsAsFactors = FALSE)
    colnames(grp.memberships) <- toupper(letters[1:ncol(grp.memberships)])

    if (type == "d3.network")
    {
        #forceNetwork(Links = edges2, Nodes = nodes.df2, Source = "source",
        #             Target = "target", Value = "value", NodeID = "label",
        #             Group = "group", opacity = 0.8, opacityNoHover = 0.5,
        #             fontSize = 12)
    } else if (type == "d3.sankey")
    {
        #sankeyNetwork(Links = edges2, Nodes = nodes.df2, Source = "source",
        #              Target = "target", Value = "value", NodeID = "label",
        #              fontSize = 16, nodePadding = 20)
    } else if (type == "igraph.tree")
    {
        ### igraph

        # nodes.df2b <- cbind.data.frame(nodes.df2, grp.memberships)
        # 
        # net <- graph.data.frame(edges2, directed=FALSE, nodes.df2b)
        # E(net)$width <- E(net)$value/4
        # V(net)$size  <- V(net)$value/1.5
        # 
        # net.h <- net - E(net)[E(net)$type=="mention"]
        # #plot(net.h, vertex.color="orange", main="Tie: Hyperlink")
        # l <- layout.fruchterman.reingold(net, niter = 100000)
        # 
        # root.nodes <- nodes.df2b$id[which(rowSums(grp.memberships) == 1)]
        # 
        # l <- layout.reingold.tilford(net)
        # grp.member.list <- lapply(grp.memberships, function(x) which(x == 1))
        # 
        # 
        # 
        # plot(net.h, vertex.color="darkgoldenrod1", edge.color = "grey40",
        #      layout=l, main="Selection Patterns",
        #      mark.groups = grp.member.list, mark.border = NA, mark.expand = 35,
        #      vertex.frame.color = NA,
        #      vertex.label.degree = 0)
    } else if (type == "d3.tree")
    {
        
        node.texts <- sapply(beta.nz.per.strata, 
                             function(bt) 
                                 paste(paste0("<b>", names(bt), "</b>: ", round(bt, 5), "<br>"), collapse = " ") )
        
        nodes.df$title = paste0("<p>Num vars selected: ", nodes.df$value,"</p>",
                                "<div class='rPartvisNetworkTooltipShowhim' style='color:blue;'> 
                                <U>Coefficients</U><div class='rPartvisNetworkTooltipShowme' 
                                style='color:black;overflow-y: scroll;height: 150px;'>",
                                node.texts, "</div></div>")

        visNetwork(nodes.df, edges) %>%
            visIgraphLayout(layout  = "layout_as_tree", 
                            physics = FALSE, 
                            type    = "full",
                            smooth  = TRUE, 
                            flip.y  = FALSE) %>%
            visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1),
                       nodesIdSelection = TRUE) %>%
            visInteraction(tooltipDelay = 0) %>%
            visNodes(font = list(size = 25), mass = 1, physics = FALSE)
    } else if (type == "ggnet.network")
    {

        # net %v% "num.vars" <- node.sizes
        # net %v% "names" <- colnames(plotMat.keep)
        # net %e% "num.vars" <- plotMat.keep
        # net %e% "edge.sizes" <- plotMat.keep.sqrt
        # 
        # 
        # 
        # ggnet2(net, size = "num.vars",
        #        max_size = 16, label = TRUE,
        #        edge.label = "num.vars",
        #        edge.size = "edge.sizes") ##, mode = "hall
    } else if (type == "diagram.tree")
    {

        # num.vars.per.node <- colSums(plotMat.ordered)
        # 
        # box.size <- 0.06
        # 
        # plotmat(plotMat.ordered, main = "Selection Patterns",
        #         pos = positions, curve = 0,
        #         name = names,
        #         arr.lwd = plotMat.ordered / (0.25 * max(plotMat.ordered)), arr.type = "circle",
        #         arr.width = plotMatZero, arr.length = plotMatZero,
        #         shadow.size = 0, box.size = box.size, box.cex = 0.9, cex = 0.9)
    }

}




#' plotting function to investigate estimated coefficients
#'
#' @param object fitted vennLasso object
#' @param s lambda value for the predictions. Only one can be specified at a time
#' @param ... other graphical parameters for the plot
#' @export
#' @examples
#' set.seed(123)
#'
#' dat.sim <- genHierSparseData(ncats = 3, nvars = 25, nobs = 200)
#'
#' fit <- vennLasso(x = dat.sim$x, y = dat.sim$y, groups = dat.sim$group.ind)
#'
#' plotCoefs(fit, s = fit$lambda[22])
#'
#'
plotCoefs <- function(object, s = NULL, ...)
{
    
    if (class(object)[1] == "cv.vennLasso")
    {
        if (is.null(s))
        {
            s = object$lambda.min
        } else if (is.character(s))
        {
            if (s == "lambda.min")
            {
                s <- object$lambda.min
            } else if (s == "lambda.1se")
            {
                s <- object$lambda.1se
            } else stop("Invalid specification for s (lambda)")
        }
        object <- object$vennLasso.fit
    }
    
    nbeta        <- object$beta
    combin.mat   <- object$condition.combinations
    combin.names <- object$combin.names
    n.conditions <- object$n.conditions
    M            <- object$n.combinations
    N            <- object$nobs
    p            <- object$nvars
    
    dimnames(nbeta) <- list(NULL, NULL, NULL)
    if(!is.null(s))
    {
        if (length(s) > 1)
        {
            s <- s[1]
            warning("multiple s values given, using first value only")
        }
        
        if (s < 0) stop("lambda must be positive")
        
        lambda  <- object$lambda
        lamlist <- lambdaInterp(lambda, s)
        nbeta   <- nbeta[,,lamlist$left,drop=TRUE]  * lamlist$frac + 
            nbeta[,,lamlist$right,drop=TRUE] * (1 - lamlist$frac)
        rownames(nbeta) <- combin.names
        colnames(nbeta) <- object$var.names
        nlam <- length(s)
    } else 
    {
        stop("must specify s, the lambda value for which to plot estimates")
    }
    
    direct.above.idx.list <- above.idx.list <- vector(mode = "list", length = M)
    
    for (c in 1:M) 
    {
        indi <- which(combin.mat[c,] == 1)
        
        # get all indices of g terms which are directly above the current var
        # in the hierarchy, ie. for three categories A, B, C, for the A group,
        # this will return the indices for AB and AC. For AB, it will return
        # the index for ABC. for none, it will return the indices for
        # A, B, and C, etc
        inner.loop <- (1:(M))[-c]
        for (j in inner.loop) 
        {
            diffs.tmp <- combin.mat[j,] - combin.mat[c,]
            if (all( diffs.tmp >= 0 )) 
            {
                if (sum( diffs.tmp == 1 ) == 1) 
                {
                    direct.above.idx.list[[c]] <- c(direct.above.idx.list[[c]], j)
                }
                above.idx.list[[c]] <- c(above.idx.list[[c]], j)
            }
        }
        above.idx.list[[c]] <- c(above.idx.list[[c]], c)
    }
    rsc <- rowSums(combin.mat)
    if (any(rsc == 0)) {
        above.idx.list[[which(rsc == 0)]] <- which(rsc == 0)
    }
    
    
    
    plotMat <- plotMatZero <- array(0, dim = rep(M, 2))
    colnames(plotMat) <- rownames(plotMat) <- combin.names
    
    
    edges <- data.frame(array(NA, dim = c(length(unlist(direct.above.idx.list)), 3)))
    colnames(edges) <- c("from", "to", "value")
    nodes <- NULL
    coefs <- varnames <- list()
    
    e.ct <- 0
    for (c in 1:M)
    {
        num.2.plot <- sum(nbeta[c, ] != 0)
        above.idx.cur <- direct.above.idx.list[[c]]
        for (k in 1:length(direct.above.idx.list))
        {
            if (!is.null(above.idx.cur[k]) & num.2.plot)
            {
                e.ct <- e.ct + 1
                
                edges[e.ct, ] <- c(above.idx.cur[k], c, num.2.plot)
                if (num.2.plot)
                {
                    nodes <- c(nodes, c, above.idx.cur[k])
                }
            }
        }
    }
    
    
    plotMat <- plotMatZero <- array(0, dim = rep(M, 2))
    colnames(plotMat) <- rownames(plotMat) <- combin.names
    for (c in 1:M)
    {
        num.2.plot <- sum(nbeta[c, ] != 0)
        above.idx.cur <- direct.above.idx.list[[c]]
        for (k in 1:length(direct.above.idx.list))
        {
            plotMat[above.idx.cur[k], c] <- num.2.plot
        }
    }
    
    
    sums.inds <- sapply(combin.names, function(x) sum(as.numeric(strsplit(x, ",")[[1]])) )
    ord.names.decr <- order(sums.inds, decreasing = TRUE)
    
    ## all the observed numbers of conditions
    unique.cond.nums <- unique(sums.inds)
    unique.cond.nums <- unique.cond.nums[order(unique.cond.nums, decreasing = TRUE)]
    
    zero.sm.idx <- which(unique.cond.nums == 0)
    ## the number of times each number of conditions will appear in the graph
    num.nodes.per.row <- sapply(unique.cond.nums, function(x) sum(sums.inds == x))
    
    
    zero.idx <- which(sums.inds[ord.names.decr] == 0)
    if (length(zero.idx))
    {
        plotMat.ordered <- plotMat[ord.names.decr[-zero.idx], ord.names.decr[-zero.idx]]
        plotMatZero     <- plotMatZero[ord.names.decr[-zero.idx], ord.names.decr[-zero.idx]]
        names           <- combin.names[ord.names.decr][-zero.idx]
        positions       <- num.nodes.per.row[-zero.sm.idx]
    } else
    {
        plotMat.ordered <- plotMat[ord.names.decr, ord.names.decr]
        names           <- combin.names[ord.names.decr]
        positions       <- num.nodes.per.row
    }
    
    cs <- colSums(plotMat.ordered)
    rs <- rowSums(plotMat.ordered)
    csplusrs <- cs + rs
    
    nonzero.keep <- csplusrs != 0
    
    plotMat.keep <- plotMat.ordered[nonzero.keep, nonzero.keep]
    
    
    edges <- na.omit(edges)
    nodes <- unique(nodes)
    nodes <- nodes[!is.na(nodes)]
    nodes.df <- data.frame(id = nodes, label = combin.names[nodes], stringsAsFactors = FALSE)
    
    num.nz.per.strata  <- apply(nbeta, 1, function(x) sum(x != 0))
    
    # save nonzero coefficients for each subpopulations
    beta.nz.per.strata <- lapply(apply(nbeta, 1, function(x) list(x[x != 0])  ), "[[", 1)
    
    nodes.df <- nodes.df[nodes.df$label %in% combin.names[num.nz.per.strata > 0], ]
    nodes.df$value <- num.nz.per.strata[match(nodes.df$label, names(num.nz.per.strata))]
    
    # reorder subpopulations to align with how nodes are ordered
    beta.nz.per.strata <- beta.nz.per.strata[match(nodes.df$label, names(num.nz.per.strata))]
    
    all.selected.vars  <- unique(unlist(sapply(beta.nz.per.strata, names)))
    
    sums.inds <- sapply(combin.names, function(x) sum(as.numeric(strsplit(x, ",")[[1]])) )
    zero.idx  <- which(sums.inds == 0)
    
    if (length(zero.idx))
    {
        edges <- edges[edges$to != zero.idx,]
        edges <- edges[edges$from != zero.idx,]
    }
    
    for (i in 1:length(all.selected.vars))
    {
        nodes.df.cur <- nodes.df
        edges.cur    <- edges
        
        coefs.cur <- sapply(1:length(beta.nz.per.strata), 
                            function(idx) beta.nz.per.strata[[idx]][match(all.selected.vars[i], names(beta.nz.per.strata[[idx]]))])
        coefs.cur[is.na(coefs.cur)] <- 0
        
        
        nodes.df.cur$value <- unname(coefs.cur)
        nodes.df.cur$group <- all.selected.vars[i]
        
        is.zero <- nodes.df.cur$value == 0
        
        edges.cur <- edges.cur[!(edges.cur$from %in% nodes.df.cur$id[is.zero]),]
        edges.cur <- edges.cur[!(edges.cur$to %in% nodes.df.cur$id[is.zero]),]
        
        
        if (i == 1)
        {
            nodes.df.all  <- nodes.df.cur[!is.zero,,drop=FALSE]
            edges.all     <- edges.cur
        } else 
        {
            
            nodes.df.cur$id <- nodes.df.cur$id + (i - 1) * nrow(nodes.df)
            nodes.df.all <- rbind(nodes.df.all, nodes.df.cur[!is.zero,,drop=FALSE])
            
            edges.cur[,c("from", "to")] <- edges.cur[,c("from", "to")] + (i - 1) * nrow(nodes.df)
            
            edges.all <- rbind(edges.all, edges.cur)
        }
        
    }
    

    
    #visNetwork(nodes.df, edges)
    
    
    
    #################################################################################
    
    
    #################################################################################
        
    node.texts <- sapply(beta.nz.per.strata, 
                         function(bt) 
                             paste(paste0("<b>", names(bt), "</b>: ", round(bt, 5), "<br>"), collapse = " ") )
    
    
    
    edges.all$value <- 1
    
    nodes.df.all$coef <- nodes.df.all$value
    nodes.df.all$value <- abs(nodes.df.all$value)
    nodes.df.all$col <- ifelse(nodes.df.all$value >= 0, "blue", "red")
    
    nodes.df.all$title = paste0("<p>Coef: ", round(nodes.df.all$coef, 6), "</p>")
    
    visNetwork(nodes.df.all, edges.all) %>%
        visIgraphLayout(layout = "layout_as_tree", physics = FALSE, type = "full",
                        smooth = TRUE, flip.y = FALSE) %>%
        visLegend(main = "Variable") %>%
        visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1),
                   selectedBy = list(variable  = "group")
                   ) %>%
        visInteraction(tooltipDelay = 0, navigationButtons = TRUE) %>%
        visNodes(font = list(size = 25), mass = 1, physics = FALSE)
    
}


