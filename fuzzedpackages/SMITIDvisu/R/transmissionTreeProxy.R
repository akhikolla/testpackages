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

#' transmissionTreeProxy
#' @description get transmissionTreeProxy
#' @param ttid widget instance identifier
#' @param session shiny session
#' @examples \dontrun{
#' library(SMITIDvisu)
#' ## server.R
#' transmissionTreeProxy <- transmissionTreeProxyProxy("transmissionTreeoutput")
#' }
#'
#' @export
transmissionTreeProxy <- function(ttid, session = shiny::getDefaultReactiveDomain()) {
  
  if(is.null(session)) {
    stop("transmissionTreeProxy must be called from the server function of a Shiny app")
  }
  
  proxy <- list(id = ttid,
                session = session)
  class(proxy) <- "transmissionTree_proxy"
  
  return(proxy)
  
}

#' updateTransmissionTree 
#' @description update (redraw) an instance of a transmissionTree
#' @param TTProxy transmissionTreeProxy instance
#' @param nodes a data.frame that represent hosts status in time with ID, status and time in columns
#' @param edges a data.frame that represent tramsmission link between hosts (pathogens) with ID, source, weight, target and time in columns
#' @param options transmissionTree new options
#' @seealso \code{\link{transmissionTree}}
#' @examples \dontrun{
#' library(SMITIDvisu)
#' data(transmissionTree)
#' ## server.R
#' transmissionTreeProxy("transmissionTreeoutput") %>% updatetransmissionTree(tt.nodes,tt.edges)
#' }
#' @export
updateTransmissionTree <- function(TTProxy, nodes, edges, options=NULL) {
  
  message <- jsonlite::toJSON(list(id = TTProxy$id, data = createTimeGraph(nodes,edges), options=options), dataframe="rows") 
  TTProxy$session$sendCustomMessage("updateTT", message)
  
  return(TTProxy)
}
