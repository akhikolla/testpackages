############################
## Methods for "beam-class"
############################

#' @rdname beam-class
#' @aliases print
#' @param x An object of class \code{beam-class}
#' @param ... further arguments passed to or from other methods.
setMethod(
	f = "print",
	signature = "beam",
	definition = function(x,...){
		cat("Class \"beam\"\n")
	}
)

#' @rdname beam-class
#' @aliases show
#' @param object An object of class \code{beam-class}
setMethod(
	f = "show",
	signature = "beam",
	definition = function(object) print(object)
)

#' @rdname beam-class
#' @aliases summary
#' @param object An object of class \code{beam-class}
setMethod(
  f = "summary",
  signature = "beam",
  definition = function(object, ...){
    cat("\n")
    cat("Object of class \"beam\"\n\n")
    cat("+ INPUT DATA \n")
    cat("number of observations:", object@dimX[1], "\n")
    cat("number of variables:", object@dimX[2], "\n\n")
    
    cat("+ ESTIMATION \n")
    cat("Shrinkage parameter alpha =", object@alphaOpt, "\n\n")
    
    labs <- colnames(object@table)
    if(!(object@type == 'conditional')){ # could be either 'marginal' or 'both'
      cat("+ SUMMARY MARGINAL DEPENDENCIES")
      smry <- NULL
      if("m_cor" %in% labs){
        smry <- rbind(smry, format(round(quantile(object@table[,"m_cor"], probs=c(0,0.5,1)), 3), nsmall=3))
      }
      if("m_logBF" %in% labs){
        smry <- rbind(smry, format(round(quantile(object@table[,"m_logBF"], probs=c(0,0.5,1)), 3), nsmall=3))
      }
      if("m_tail_prob" %in% labs){
        smry <- rbind(smry, format(quantile(object@table[,"m_tail_prob"], probs=c(0,0.5,1)), digits=3))
      }
      rownames(smry) <- c(" marginal correlations ", " log(Bayes factors)", "tail probabilities")[c("m_cor", "m_logBF", "m_tail_prob")%in%labs]
      colnames(smry) <- c("Min.", "Median", "Max.")
      print(kable(smry, align=rep('r', 3)))
      cat("\n")
    }
    if(!(object@type == 'marginal')){  # could be either 'conditional' or 'both'
      cat("+ SUMMARY CONDITIONAL DEPENDENCIES")
      smry <- NULL
      if("p_cor"%in%labs){
        smry <- rbind(smry, format(round(quantile(object@table[,"p_cor"], probs=c(0,0.5,1)), 3), nsmall=3))
      }
      if("p_logBF"%in%labs){
        smry <- rbind(smry, format(round(quantile(object@table[,"p_logBF"], probs=c(0,0.5,1)), 3), nsmall=3))
      }
      if("p_tail_prob"%in%labs){
        smry <- rbind(smry, format(quantile(object@table[,"p_tail_prob"], probs=c(0,0.5,1)), digits=3))
      }
      rownames(smry) <- c(" partial correlations ", " log(Bayes factors)", "tail probabilities")[c("p_cor", "p_logBF", "p_tail_prob")%in%labs]
      colnames(smry) <- c("Min.", "Median", "Max.")
      print(kable(smry, align=rep('r', 3)))
      cat("\n")
    }
    cat("+ ESTIMATED COMPUTING TIME (IN SECONDS): ", object@time, "\n\n")
  }
)

#' @rdname beam-class
#' @aliases marg
#' @param object An object of class \code{beam-class}
setMethod(
  f = "marg",
  signature = "beam",
  definition = function(object){
    
    # Check input
    assert_that(object@type != 'conditional', msg="no information available in input beam object")
    
    # Extract table
    matidxs <- .upperTriIdxs(object@dimX[2])
    rowcol <- data.frame(row = matidxs[,1], col = matidxs[,2])
    df <- object@table
    marg_cols <- c('m_cor','m_logBF','m_tail_prob')
    df <- df[, intersect(marg_cols, colnames(df))]
    return(cbind(rowcol,df))
  }
)

#' @rdname beam-class
#' @aliases cond
#' @param object An object of class \code{beam-class}
setMethod(
  f = "cond",
  signature = "beam",
  definition = function(object){
    
    # Check input
    assert_that(object@type != 'marginal', msg="no information available in input beam object")
    
    # Extract table
    df <- object@table
    matidxs <- .upperTriIdxs(object@dimX[2])
    rowcol <- data.frame(row = matidxs[,1], col = matidxs[,2])
    cond_cols <- c('p_cor','p_logBF','p_tail_prob')
    df <- df[,intersect(cond_cols, colnames(df))]
    
    return(cbind(rowcol,df))
  }
)

#' @rdname beam-class
#' @aliases mcor
#' @param object An object of class \code{beam-class}
setMethod(
  f = "mcor",
  signature = "beam",
  definition = function(object){
    
    # Check input
    assert_that('m_cor' %in% colnames(object@table), msg="'m_cor' not available in input beam object")
    
    # Extract marginal correlation matrix
    p <- object@dimX[2]
    idxs <- .upperTriIdxs(p)
    matcor <- as.matrix(Matrix::sparseMatrix(i=c(idxs[,1],idxs[,2]), j=c(idxs[,2],idxs[,1]), x=rep(object@table[,"m_cor"],2), dims=c(p,p)))
    diag(matcor) <- 1
    
    return(matcor)
  }
)

#' @rdname beam-class
#' @aliases pcor
#' @param object An object of class \code{beam-class}
setMethod(
  f = "pcor",
  signature = "beam",
  definition = function(object){
    
    # Check input
    assert_that('p_cor' %in% colnames(object@table), msg="'p_cor' not available in input beam object")
    
    # Extract partial correlation matrix
    p <- object@dimX[2]
    idxs <- .upperTriIdxs(p)
    matpcor <- as.matrix(Matrix::sparseMatrix(i=c(idxs[,1],idxs[,2]), j=c(idxs[,2],idxs[,1]), x=rep(object@table[, "p_cor"],2), dims=c(p, p)))
    diag(matpcor) <- 1
    
    return(matpcor)
  }
)

#' @rdname beam-class
#' @aliases postExpSigma
#' @param object An object of class \code{beam-class}
#' @param vars.method method of shrinkage estimation for variances 
setMethod(
  f = "postExpSigma",
  signature = "beam",
  definition = function(object, vars.method="eb"){
    
    # Check input
    assert_that('m_cor' %in% colnames(object@table), msg="'m_cor' not available in input beam object")
    assert_that(is.character(vars.method))
    assert_that(length(vars.method)==1, msg="vars.method must be of length 1")
    assert_that(vars.method %in% c("eb", "mean", "median", "none", "scaled"), msg="unkown vars.method")
    
    # Posterior expectation
    covmat <- mcor(object)
    
    # Rescale or not
    if(vars.method != "scaled"){
      s2 <- .shrinkvars(object@s, object@dimX[1], method=vars.method)
      covmat <- covmat * tcrossprod(sqrt(s2))
    }
    
    return(covmat)
  }
)


#' @rdname beam-class
#' @aliases postExpOmega
#' @param object An object of class \code{beam-class}
#' @param vars.method method of shrinkage estimation for variances 
setMethod(
  f = "postExpOmega",
  signature = "beam",
  definition = function(object, vars.method="eb"){

    # Check input
    assert_that('p_cor' %in% colnames(object@table), msg="'p_cor' not available in input beam object")
    assert_that(is.character(vars.method))
    assert_that(length(vars.method)==1, msg="vars.method must be of length 1")
    assert_that(vars.method %in% c("eb", "mean", "median", "none", "scaled"), msg="unkown vars.method")
    
    # Posterior expectation
    icovmat <- - pcor(object)
    diag(icovmat) <- 1
    icovmat <- icovmat * tcrossprod(object@TinvStdev)
    icovmat <- (object@dimX[1] + object@deltaOpt) * icovmat
    
    # Rescale or not
    if(vars.method != "scaled"){

      s2 <- .shrinkvars(object@s, object@dimX[1], method=vars.method)
      icovmat <- icovmat * tcrossprod(1/sqrt(s2))
      
    }
    
    return(icovmat)
  }
)

#' @rdname beam-class
#' @aliases plotML
setMethod(
  f = "plotML",
  signature = "beam",
  definition = function(object, ...){
    plot(object@gridAlpha[,2], object@gridAlpha[,3], type="l", xlab=expression(alpha), ylab="log-marginal likelihood", lwd=2, ...)
    abline(v=object@alphaOpt, col="black", lty=2)
    abline(h=object@valOpt, col="black", lty=2)
  }
)

#' @rdname beam-class
#' @aliases plotCor
#' @param type character. Type of correlation to be displayed (marginal, conditional or both)
#' @param order character. "original" or "clust"
#' @param by character. When order = 'clust', this argument specifies whether the clustering has to be performed on marginal or partial correlations
setMethod(
  f = "plotCor",
  signature = "beam",
  definition = function(object, type=object@type, order = "original", by = 'marginal'){
    if(!type%in%c("marginal", "conditional", "both")){
      stop("type must be 'marginal', 'conditional' or 'both' ")
    }else{
      if(!order%in%c('original', 'clust')){
        stop("order must be equal to 'original' or 'clust'")
      }else{
        
        # Extract correlations
        if(type=="marginal"){
          if(object@type %in% c('both', 'marginal')){
            if(order == "clust"){
              bgr <- bgraph(object)
              clust <- igraph::cluster_louvain(bgr)
              memb <- igraph::membership(clust)
              ordered.idxs <- order(as.vector(memb))
              themat <- mcor(object)
              themat <- themat[ordered.idxs,ordered.idxs]
            }else{
              themat <- mcor(object)
            }
          }else{
            stop('In beam call type = "conditional". No information about marginal')
          }
        }
  
        if(type=="conditional"){
          if(object@type %in% c('both', 'conditional')){
            if(order == "clust"){
              ugr <- ugraph(object)
              clust <- igraph::cluster_louvain(ugr)
              memb <- igraph::membership(clust)
              ordered.idxs <- order(as.vector(memb))
              themat <- pcor(object)
              themat <- themat[ordered.idxs,ordered.idxs]
            }else{
              themat <- pcor(object)
            }
          }else{
            stop('In beam call type = "marginal". No information about conditional')
          }
        }
  
        if(type=="both"){
          if(object@type == 'both'){
            if(order == "clust"){
              if(by == 'marginal'){
                bgr <- bgraph(object)
                clust <- igraph::cluster_louvain(bgr)
                memb <- igraph::membership(clust)
                ordered.idxs <- order(as.vector(memb))
                mymatcor <- mcor(object)
                mymatcor <- mymatcor[ordered.idxs,ordered.idxs] # reorder columns according to the clusters, as specified in membership
                mymatpcor <- pcor(object)
                mymatpcor <- mymatpcor[ordered.idxs,ordered.idxs] # reorder columns according the order imposed by clustering on marginal associations
                themat <- mymatcor
                themat[upper.tri(themat, diag=FALSE)] <- mymatpcor[upper.tri(mymatpcor, diag=FALSE)]
              }else{    # by = 'conditional'
                ugr <- ugraph(object)
                clust <- igraph::cluster_louvain(ugr)
                memb <- igraph::membership(clust)
                ordered.idxs <- order(as.vector(memb))
                mymatpcor <- pcor(object)
                mymatpcor <- mymatpcor[ordered.idxs,ordered.idxs] # reorder columns according the order imposed by clustering on marginal associations
                mymatcor <- mcor(object)
                mymatcor <- mymatcor[ordered.idxs,ordered.idxs] # reorder columns according the order imposed by clustering on conditional associations
                themat <- mymatcor
                themat[upper.tri(themat, diag=FALSE)] <- mymatpcor[upper.tri(mymatpcor, diag=FALSE)]
              }
            }else{
              mymatcor <- mcor(object)
              mymatpcor <- pcor(object)
              themat <- mymatcor
              themat[upper.tri(themat, diag=FALSE)] <- mymatpcor[upper.tri(mymatpcor, diag=FALSE)]
            }
          }else{
            stop('In beam call type is not "both"')
          }
        }
        if(order == "clust"){
          cat('Number of clusters identified with algorithm "cluster_louvain":', max(memb), '\n')
          cat('Number of nodes in each cluster:', igraph::sizes(clust), '\n')
          cat('Modularity =', igraph::modularity(clust))
        }
  
        # Function needed
        mirror <- function(mymat){
          xx <- as.data.frame(mymat);
          xx <- rev(xx);
          xx <- as.matrix(xx);
          xx;
        }
  
        # Colors
        #f <- colorRampPalette(c("#ff0000", "#FFFFFF", "#0000ff"))
        f <- colorRampPalette(c("red","white","blue"))
        rg <- f(200)
  
        diag(themat) <- NA
        
        # Plot
        par(mar=c(2, 2, 2, 2) + 0.1)
        #image(1:nrow(themat), 1:ncol(themat), mirror(themat), zlim=c(-1,1), col=rg, xlab="", ylab="", main="", xaxt="n", yaxt="n")
        image(1:nrow(themat), 1:ncol(themat), mirror(themat), zlim=c(min(themat, na.rm=TRUE), max(themat, na.rm=TRUE)), col=rg, xlab="", ylab="", main="", xaxt="n", yaxt="n")
        #box()
      }
    }
  }
)

#' @rdname beam-class
#' @aliases bgraph
setMethod(
  f = "bgraph",
  signature = "beam",
  definition = function(object){

    # Check input
    assert_that('m_cor' %in% colnames(object@table), msg="'m_cor' not available in input beam object")
    
    # Get bidirected graph
    p <- object@dimX[2]
    idxs <- .upperTriIdxs(p)
    edges <- as.data.frame(idxs)
    edges <- cbind(edges, abs(object@table[,"m_cor"]))
    colnames(edges) <- c('node1','node2','weight')
    myigraph <- igraph::graph_from_data_frame(d=edges, directed=FALSE)
    if(length(object@varlabs)>0){
      myigraph <- igraph::set.vertex.attribute(myigraph, "name", value=object@varlabs)
    }
    
    return(myigraph)
  }
)

#' @rdname beam-class
#' @aliases ugraph
setMethod(
  f = "ugraph",
  signature = "beam",
  definition = function(object){

    # Check input
    assert_that('p_cor' %in% colnames(object@table), msg="'p_cor' not available in input beam object")
    
    # Get undirected graph
    p <- object@dimX[2]
    idxs <- .upperTriIdxs(p)
    edges <- as.data.frame(idxs)
    edges <- cbind(edges, abs(object@table[,"p_cor"]))
    colnames(edges) <- c('node1','node2','weight')
    myigraph <- igraph::graph_from_data_frame(d=edges, directed=FALSE)
    if(length(object@varlabs)>0){
      myigraph <- igraph::set.vertex.attribute(myigraph, "name", value=object@varlabs)
    }
    
    return(myigraph)

  }
)
