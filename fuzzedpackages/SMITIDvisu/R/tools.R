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

#' createRainbowColors
#' Create a list of colors for each value v
#' @param v a vector of characters
#' @return a list of value=color 
#' @importFrom stats setNames
createRainbowColors <- function(v) {
    
    colorList <- c("#FF0000","#800000","#FFFF00","#808000","#00FF00", "#008000","#00FFFF",
                   "#008080","#0000FF","#000080","#FF00FF","#800080","#000000","#C0C0C0","#808080")
    
    v <- unique(v)
    v <- v[!is.na(v)]
    #r <- as.list(rainbow(length(v),v=0.8))
    r <- as.list(colorList[1:length(v)])
    lc <- setNames( r, v)
    
    return(lc)
}
