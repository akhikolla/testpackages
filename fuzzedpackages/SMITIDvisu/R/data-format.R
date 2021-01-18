# This file is part of SMITIDvisu package.
# Copyright (C) 2018-2019 Jean-François Rey <jean-francois.rey@inra.fr>
#
# SMITIDvisu is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SMITIDvisu is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SMITIDvisu. If not, see <https://www.gnu.org/licenses/>.
#


# createTimeGraph
# @description create a list of graphs over the time
# @details create a list of list graph over the time  
# a graph (output) is represented by a nodes element as a data.frame, an edges element as data.frame and an element time.  
# graphs = list( graph1=list(nodes,edges,time), graph2=list(nodes,edges,time))
# @param nodes data.frame of nodes ID|time|status
# @param edges data.frame of edges ID|source|target|time
# @return a list of graph sort by time
createTimeGraph <- function(nodes,edges) {
  
  graphs.list <- vector("list",0)
  
  # converte Date ISO POSIX to timestamp
  if(  any(is.na(as.numeric(nodes$time[which(!is.na(nodes$time))]))) ) {
    nodes$time <- as.numeric(as.POSIXct(strptime(nodes$time, format="%Y-%m-%dT%H:%M:%S")))*1000
    edges$time <- as.numeric(as.POSIXct(strptime(edges$time, format="%Y-%m-%dT%H:%M:%S")))*1000
  }
  else {
      # is a timestamp
      if( as.numeric(format(as.Date(as.POSIXct(as.numeric(nodes$time[which(!is.na(nodes$time))][1]), origin="1970-01-01", tz = "UTC")), format="%Y")) >= 1971 ) {
          nodes$time <- as.numeric(nodes$time)*1000
          edges$time <- as.numeric(edges$time)*1000 
      }
      # julian day
      else {
        nodes$time <- as.numeric(nodes$time)
        edges$time <- as.numeric(edges$time)
      }
  }
  
  # all times
  timeline <- sort(unique(floor(c(as.numeric(nodes$time), as.numeric(edges$time)))))
  
  # List unique nodes
  nodes.init <- as.data.frame(cbind(unique(as.character(nodes$ID)),NA),row.names=NULL)
  names(nodes.init) <- c("ID","status")
  
  # List of nodes at first time
  nodes.start <- nodes[which(floor(nodes$time) == min(timeline)), c("ID","status")]
  nodes.start <- rbind(nodes.init[which( !(nodes.init$ID %in% nodes.start$ID)),], nodes.start)
  
  
  
  # initial graph as list
  graph.start <- list("nodes"=nodes.start,"edges"=edges[which(floor(edges$time) == min(timeline)), c("ID","source","target","weight")],"time"=min(timeline))
  graphs.list <- list(c(graphs.list, graph.start))
  
  # append next graph in time
  for( t in timeline[2:length(timeline)]) {
    nodes.next <- subset(nodes, floor(nodes$time) %in% t, select=c("ID","status"))
    edges.next <- subset(edges, floor(edges$time) %in% t, select=c("ID","source","target","weight"))
    graph.next <- list("nodes"=nodes.next,"edges"=edges.next,"time"=t)
    graphs.list <- c(graphs.list,list(graph.next))
  }
  
  return(graphs.list)
}

# createTimeLine
# @description create a time line list
# @details create a timeline list
# Set info mintime, maxtime and nblevels
# associate level "top, "middle" and "bottom" to -1, 0 and 1 value
# associate each items by levels
# @param items data.frame of items to plot
# @param title a title (host ID)
# @return a list of list items timeline
createTimeLine <- function(items, title) {
  if( nrow(items) == 0 || is.null(items$level) || is.null(items$label) || is.null(items$timestart)) {return(list())}
  # converte Date ISO POSIX to timestamp
  if( any(is.na(as.numeric(items$timestart[which(!is.na(items$timestart))]))) ) {
    items$timestart <- as.numeric(as.POSIXct(strptime(items$timestart, format="%Y-%m-%dT%H:%M:%S"), origin="1970-01-01", tz = "UTC"))*1000
    items$timeend <- as.numeric(as.POSIXct(strptime(items$timeend, format="%Y-%m-%dT%H:%M:%S"), origin="1970-01-01", tz = "UTC"))*1000
  }
  else {
      # is a timestamp
      if( as.numeric(format(as.Date(as.POSIXct(as.numeric(items$timestart[which(!is.na(items$timestart))][1]), origin="1970-01-01", tz = "UTC")), format="%Y")) >= 1971 ) {
          items$timestart <- as.numeric(items$timestart)*1000
          items$timeend <- as.numeric(items$timeend)*1000 
      }
  }
  
  mintime <- min(as.numeric(as.character(items$timestart)), na.rm = TRUE)
  maxtime <- max(as.numeric(as.character(items$timestart)),as.numeric(as.character(items$timeend))[which(as.numeric(as.character(items$timeend)) < Inf)],na.rm = TRUE) + 1
  nblevels <- length(unique(items$level))
  # 
  # hosttimeline <- list(list("type"="timeline", "mintime"=mintime, "maxtime"=maxtime, "nblevels"=nblevels),
  #                   list("level"="-1", "label"="input", "timeline" = (items[which(items$level == "top"),2:5])),
  #                   list("level"="0", "label"=title, "timeline" = (items[which(items$level == "middle"),2:5])),
  #                   list("level"="1", "label"="output", "timeline" = (items[which(items$level == "bottom"),2:5])))
  
  maxlength <- ncol(items)
  top <- list()
  bottom <- list()
  middle <- list()
  if(nrow(items[which(items$level == "top"),2:maxlength]) > 0) top <- list("level"="-1", "label"="input", "timeline" = (items[which(items$level == "top"),2:maxlength]))
  if(nrow(items[which(items$level == "middle"),2:maxlength]) > 0) {
    middle <- list("level"="0", "label"=title, "timeline" = (items[which(items$level == "middle"),2:maxlength]))
    middle <- rapply(middle,function(x) ifelse(x==Inf,maxtime,x), how = "replace")
  }
 
  if(nrow(items[which(items$level == "bottom"),2:maxlength]) > 0) bottom <- list("level"="1", "label"="output", "timeline" = (items[which(items$level == "bottom"),2:maxlength]))
  
  hosttimeline <- list(list("type"="timeline", "mintime"=mintime, "maxtime"=maxtime, "nblevels"=nblevels), top, middle, bottom)
  
  return(hosttimeline)
}

# createMSTGraph
# @description compute the minimum spanning tree graph 
# @details create a list of nodes and edges representing the minimum spanning tree.
# use the matrix to compute the MST.
# use matrix colnames as nodes Id and use prop as the nodes weight.
# edges source and traget are node Id and the weight represente the distance in the matrix.
# @param mat a distance matrix
# @param prop an array of nodes proportions
# @return a list of nodes and edges list.
createMSTGraph <- function(mat, prop) {
  
  mst <- mstCompute(mat)
  
  mst <- ifelse(mst>0,mat,NA)
  mst[lower.tri(mst,diag=FALSE)] <- 0

  nodes <- as.list(as.vector(
              apply(prop,1,function(l) {
                              return( list( "ID"=l[c("ID")][[1]], "weight"=l[c("proportion")][[1]], "count"=l[c("count")][[1]] ));
                            })
           ))

  ## durty code... 
  edges = vector("list",0)
  for( i in 1:nrow(mst) ) {
    for( j in (i):ncol(mst) ) {
      if(!is.na(mst[i,j]) & mst[i,j] > 0) {
          edges <- c(edges, list(list("source"=rownames(mst)[i] ,"target"=colnames(mst)[j], "weight"=mst[i,j]))) 
      }
    }
  }
  
  graphs.list <- list("nodes"=unname(nodes), "edges"=unname(edges))
  
  return(graphs.list)
}
