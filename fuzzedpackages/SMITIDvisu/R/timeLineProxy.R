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


#' timeLineProxy
#' get an instance of a timeline
#' @param tlid a timeline instance id
#' @param session  shiny session
#' @return an object of class timeline_proxy
#' @examples
#' \dontrun{
#' ## server.R
#' ## output server variable
#' output$timeline <- renderTimeLine({
#'          timeLine(data.frame(), "")
#'        })
#' ## ui.R
#' timeLineOutput("timeline")
#' ## server.R
#' tlproxy <- timeLineProxy("timeline")
#' }
#' @export
timeLineProxy <- function(tlid, session = shiny::getDefaultReactiveDomain()) {
  
  if(is.null(session)) {
    stop("timeLineProxy must be called from the server function of a Shiny app")
  }
  
  proxy <- list(id = tlid,
                session = session)
  class(proxy) <- "timeLine_proxy"
  
  return(proxy)
  
}

#' updateTimeLine
#' @param tlProxy a timeline proxy instance
#' @param data new data
#' @param title new title
#' @seealso \code{\link{timeLine}}
#' @examples
#' \dontrun{
#' ## server.R
#' ## output server variable
#' output$timeline <- renderTimeLine({
#'          timeLine(data.frame(), "")
#'        })
#' ## ui.R
#' timeLineOutput("timeline")
#' ## server.R
#' timeLineProxy("timeline") %>% updateTimeLine(newtimeline, "newId")
#' }
#' @export
updateTimeLine <- function(tlProxy, data, title) {
  
  message <- jsonlite::toJSON(list(id = tlProxy$id, data = createTimeLine(data,title)), dataframe="rows")
  tlProxy$session$sendCustomMessage("updatetl", message)
  
  return(tlProxy)
}
