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

#' @title demo.SMITIDvisu.run
#' @description run a demo to visualize data
#' @export
demo.SMITIDvisu.run <- function() {

  library(SMITIDvisu)

  ## Load transmissiontree data
  tt.edges <- get("tt.edges")
  tt.nodes <- get("tt.nodes")
  #data("transmissiontree", package="SMITIDvisu")
  ## Draw a transmission tree over the time as a circle.
  tt <- transmissionTree(tt.nodes, tt.edges, nodes.color = list("default"="black","Inf"="red"))
  htmlwidgets::saveWidget(tt, "transmissionTree.html")
  browseURL("transmissionTree.html")

  ## Load hostline data (a host timeline informations)
  hostline <- get("hostline")
  #data("hostline", package="SMITIDvisu")
  ## Draw a time line
  tl <- timeLine(hostline, title="Example 113", color=list("infected"="red","offspring"="green","alive"="blue","inf"="orange","dead"="black","Obs"="purple"))
  
  htmlwidgets::saveWidget(tl, "timeline.html")
  browseURL("timeline.html")
  
  ## Load st data, contains several variables
  st.dist113_2 <- get("st.dist113_2")
  st.prop113_2 <- get("st.prop113_2")
  st.dist113_all <- get("st.dist113_all")
  st.prop113_all <- get("st.prop113_all")
  st.listTimeProp113 <- get("st.listTimeProp113")
  #data("st", package="SMITIDvisu")

  ## Draw a variant graph of the host 113 at time 2
  mv <- mstVariant(st.dist113_2,st.prop113_2)
  htmlwidgets::saveWidget(mv, "mstVariants.html")
  browseURL("mstVariants.html")

  # With variant proportions over all the observed time
  mv <- mstVariant(st.dist113_all, st.prop113_all, st.listTimeProp113)
  htmlwidgets::saveWidget(mv, "mstVariantsAll.html")
  browseURL("mstVariantsAll.html")
}
