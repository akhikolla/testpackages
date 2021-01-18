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

#' mstVariantProxy
#' @description get mstVariantProxy
#' @param mstVid widget instance identifier
#' @param session shiny session
#' @examples \dontrun{
#' library(SMITIDvisu)
#' ## server.R
#' mstVariantProxy <- mstVaraintProxy("mstvariantoutput")
#' }
#'
#' @export
mstVariantProxy <- function(mstVid, session = shiny::getDefaultReactiveDomain()) {
  
  if(is.null(session)) {
    stop("mstVaraintProxy must be called from the server function of a Shiny app")
  }
  
  proxy <- list(id = mstVid,
                session = session)
  class(proxy) <- "mstVariant_proxy"
  
  return(proxy)
  
}

#' updatemstVariant 
#' @description update (redraw) an instance on mstVariant
#' @param mstVProxy mstVaraintProxy instance
#' @param mat distance matrix
#' @param prop proportions data.frame
#' @param propTime list of each variant by time and proportions
#' @seealso \code{\link{mstVariant}}
#' @examples \dontrun{
#' library(SMITIDvisu)
#' data(mstVariant)
#' ## server.R
#' mstVaraintProxy("mstvariantoutput") %>% updatemstVariant(st.dist,st.prop)
#' }
#' @export
updatemstVariant <- function(mstVProxy, mat, prop, propTime=NULL) {
  
  if( is.null(propTime) ){ message <- jsonlite::toJSON(list(id = mstVProxy$id, data = createMSTGraph(mat,prop)), dataframe="rows") }
  else { message <- jsonlite::toJSON(list(id = mstVProxy$id, data = c(createMSTGraph(mat,prop), "node_prop"=list(unname(propTime)))), dataframe="rows") }
  #attr(message, 'TOJSON_ARGS') <- list(dataframe="columns")
  mstVProxy$session$sendCustomMessage("updatemstV", message)
  
  return(mstVProxy)
}
