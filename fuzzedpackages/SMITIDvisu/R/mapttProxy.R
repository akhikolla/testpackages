#' mapttProxy
#' @description get mapttProxy
#' @param mapttId widget instance identifier
#' @param session shiny session
#' @examples \dontrun{
#' library(SMITIDvisu)
#' ## server.R
#' mapttProxy <- mapttProxyProxy("mapttOutput")
#' }
#'
#' @export
mapttProxy <- function(mapttId, session = shiny::getDefaultReactiveDomain()) {

  if (is.null(session)) {
    stop("mapttProxy must be called from the server function of a Shiny app")
  }

  proxy <- list(id = mapttId,
                session = session)
  class(proxy) <- "maptt_proxy"

  return(proxy)

}

#' mapttSelectHost
#' @description select a host on a MapTT instance
#' @param mapttProxy mapttProxy instance
#' @param hostId the id of the host to select
#' @seealso \code{\link{maptt}}
#' @examples \dontrun{
#' library(SMITIDvisu)
#' data(transmissiontree)
#' ## server.R
#' mapttProxy("mapttOutput") %>% mapttSelectHost()
#' }
#' @export
mapttSelectHost <- function(mapttProxy, hostId) {
  message <- jsonlite::toJSON(list(id = mapttProxy$id, hostId = hostId))
  mapttProxy$session$sendCustomMessage("selectHost", message)
}
