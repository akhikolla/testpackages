# Copyright 2015-2019 Venelin Mitov
#
# This file is part of POUMM.
#
# POUMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# POUMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with POUMM  If not, see <http://www.gnu.org/licenses/>.

# Tree-processing functions

NULL


#' Check if the POUMM version correpsonds to a dev release
#' @importFrom utils packageDescription
#' @description We define a dev release as having a sub-release, eg 0.9.15.5 is
#' one whereas 0.9.16 is not.
#' @return a logical
#' @export
POUMMIsADevRelease <- function() {
  # !is.na( packageDescription("POUMM") ) &&
  #   length(strsplit(packageDescription("POUMM")$Version, "\\.")[[1]]) >= 3
  TRUE
}
printPOUMMLikelihoodMainLoop <- function(tree) {
  pruneInfo <- pruneTree(tree)
  M <- pruneInfo$M
  endingAt <- pruneInfo$endingAt
  nodesVector <- pruneInfo$nodesVector
  nodesIndex <- pruneInfo$nodesIndex
  nLevels <- pruneInfo$nLevels
  unVector <- pruneInfo$unVector
  unIndex <- pruneInfo$unIndex
  
  unJ <- 1
  N <- length(tree$tip.label)
  edge <- rbind(pruneInfo$edge, c(0,1))
  
  for(i in 1:(nLevels + 1)) {
    if(i <= nLevels) {
      es <- endingAt[nodesVector[(nodesIndex[i] + 1):nodesIndex[i + 1]]]  
    } else {
      es <- M
    }
    
    if(es[1] != M) {
      cat("Processing edges ending at nodes:\n")
      print(edge[es,2])
    } else {
      cat("Processing root node.\n")
    }
    
    if(es[1] != M) {
      lenUnAll <- 0L
      
      while(lenUnAll != length(es)) {
      # while(length(es)>0) {
        un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
        unJ <- unJ+1
        cat("Updating parent-nodes:\n")
        
        print(cbind(edge[es[un], 1], " <- ", edge[es[un], 2]))
        
        lenUnAll <- lenUnAll + length(un)
        
        #es <- es[-un]
      }  
    }
    
  }
}

#' @importFrom graphics par plot
plotPOUMMLikelihoodMainLoop <- function(tree, x.lim=c(-2,16.4)) {
  pruneInfo <- pruneTree(tree)
  M <- pruneInfo$M
  endingAt <- pruneInfo$endingAt
  nodesVector <- pruneInfo$nodesVector
  nodesIndex <- pruneInfo$nodesIndex
  nLevels <- pruneInfo$nLevels
  unVector <- pruneInfo$unVector
  unIndex <- pruneInfo$unIndex
  
  unJ <- 1
  N <- length(tree$tip.label)
  
  edge <- rbind(pruneInfo$edge, c(0,1))
  
  edgeColor <- rep("black", nrow(edge))
  tipBg <- "white"
  tipColor <- "red"
  tipFrame <- "circle"
  
  nodeBg <- rep("black", length(tree$node.label))
  nodeCol <- rep("white", length(tree$node.label))
  
  names(nodeBg) <- names(nodeCol) <- tree$node.label
  
  for(i in 1:(nLevels + 1)) {
    if(i <= nLevels) {
      es <- endingAt[nodesVector[(nodesIndex[i] + 1):nodesIndex[i + 1]]]  
    } else {
      es <- M
    }
    
    edgeColor[es] <- "red"
    
    plot(tree, type="cladogram", show.tip.label=FALSE, edge.color = edgeColor, edge.width=2, x.lim=x.lim)
    
    if(tipBg == "darkgrey") {
      parColDefault <- par(col="darkgrey")
    } else {
      parColDefault <- par(col="red")
    }
    tiplabels(text = tree$tip.label, 
              frame = tipFrame, 
              bg = tipBg, 
              col = tipColor, cex = 1.4)
    par(col=parColDefault)
    
    # update these for the next iteration
    tipBg <- "darkgrey"
    tipColor <- "darkgrey"
    
    if(i > 1) {
      # not at tip-level
      nodeBg[as.character(edge[es, 2] - 1)] <- "white"
      nodeCol[as.character(edge[es, 2] - 1)] <- "red"
    }
    
    if(any(nodeCol == "red")) {
      parColDefault <- par(col="red")
      nodelabels(text = tree$node.label[nodeCol=="red"],
                 node = which(nodeCol == "red") + N,
                 frame = "circle",
                 bg = nodeBg[nodeCol == "red"],
                 col=nodeCol[nodeCol == "red"],
                 cex = 1)
      par(col=parColDefault)
    }
    if(any(nodeCol == "darkgrey")) {
      parColDefault <- par(col="darkgrey")
      nodelabels(text = tree$node.label[nodeCol=="darkgrey"],
                 node = which(nodeCol == "darkgrey") + N,
                 frame = "circle",
                 bg = nodeBg[nodeCol == "darkgrey"],
                 col = nodeCol[nodeCol == "darkgrey"],
                 cex = 1)
      par(col=parColDefault)
    }
    if(any(nodeCol == "white")) {
      parColDefault <- par(col="black")
      nodelabels(text = tree$node.label[nodeCol=="white"],
                 node = which(nodeCol == "white") + N,
                 frame = "circle",
                 bg = nodeBg[nodeCol == "white"],
                 col = nodeCol[nodeCol == "white"],
                 cex = 1)
      par(col=parColDefault)
    }
    
    if(i > 1) {
      # not at tip-level
      nodeBg[as.character(edge[es, 2] - 1)] <- "darkgrey"
      nodeCol[as.character(edge[es, 2] - 1)] <- "darkgrey"
    }
    
    edgeColor[es] <- "darkgrey"
    
    if(es[1] != M) {
      cat("Processing edges ending at nodes:\n")
      print(edge[es,2])
    } else {
      cat("Processing root node.\n")
    }
    
    if(es[1] != M) {
      queueNo <- 1
      letters <- c('a','b','c')
      
      #lenUnAll <- 0L
      
      #while(lenUnAll != length(es)) {
      while(length(es)>0) {
        un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
        unJ <- unJ+1
        cat("Updating parent-nodes:\n")
        
        print(cbind(edge[es[un], 1], " <- ", edge[es[un], 2]))
        parColDefault <- par(col="yellow")
        edgelabels(paste0(if(i!=3) "   " else "", letters[queueNo], if(i < nLevels) " " else ""), es[un], cex = 1, 
                   frame="none", col="red", 
                   adj=c(0.6,1.4))
        par(col=parColDefault)
        
        #lenUnAll <- lenUnAll + length(un)
        
        es <- es[-un]
        queueNo <- queueNo + 1
      }  
    }
    
  }
}

#' Writes verbose messages of the order of tree traversal during likelihood calculation
#' @param tree A phylo object.
#' @return Nothing
#' @export
simulatePOUMMLikelihoodMainLoop <- function(tree) {
  pruneInfo <- pruneTree(tree)
  M <- pruneInfo$M
  endingAt <- pruneInfo$endingAt
  nodesVector <- pruneInfo$nodesVector
  nodesIndex <- pruneInfo$nodesIndex
  nLevels <- pruneInfo$nLevels
  unVector <- pruneInfo$unVector
  unIndex <- pruneInfo$unIndex
  edge <- pruneInfo$edge
  unJ <- 1
  
  
  for(i in 1:nLevels) {
    es <- endingAt[nodesVector[(nodesIndex[i] + 1):nodesIndex[i + 1]]]
    
    cat("Processing edges ending at nodes:\n")
    print(edge[es,2])
    
    while(length(es)>0) {
      un <- unVector[(unIndex[unJ]+1):unIndex[unJ+1]]
      unJ <- unJ+1
      cat("Updating parent-nodes:\n")
      print(cbind(edge[es[un], 1], " <- ", edge[es[un], 2]))
    
      es <- es[-un]
    }
  }
}

#' Node indices of the direct descendants of n in the phylogeny.
#' 
#' @param tree an object of class phylo
#' @param n an index of a node (root, internal or tip) in tree
#' @return An integer vector.
chld <- function(tree, n) {
  as.integer(tree$edge[tree$edge[, 1]==n, 2])
}

#' Edge indices of the edges in tree starting from n
#' @param tree an object of class phylo
#' @param n an index of a node (root, internal or tip) in tree
#' @return An integer vector.
edgesFrom <- function(tree, n) {
  which(tree$edge[, 1]==n)
}

#' Calculate the time from the root to each node of the tree
#' 
#' @param tree An object of class phylo.
#' @param tipsOnly Logical indicating whether the returned results should be
#'   truncated only to the tips of the tree.
#' @return A vector of size the number of nodes in the tree (tips, root, 
#'   internal) containing the time from the root to the corresponding node in 
#'   the tree.
#' @importFrom stats reorder
#' @import ape
#' @export
nodeTimes <- function(tree, tipsOnly=FALSE) {
  rtree <- reorder(tree, 'postorder')
  es <- rtree$edge[dim(rtree$edge)[1]:1, ]
  nEdges <- dim(es)[1]
  ts <- rev(rtree$edge.length)
  nodeTimes <- rep(0, length(rtree$tip.label)+rtree$Nnode)
  for(e in 1:nEdges) 
    nodeTimes[es[e, 2]] <- nodeTimes[es[e, 1]]+ts[e]
  if(tipsOnly) {
    nodeTimes[1:length(tree$tip.label)]
  } else {
    nodeTimes
  }
}
