#' Decode Draco encoded raw bytes containing mesh or point cloud data
#'
#' @param data \code{\link{raw}} bytes containing Draco data e.g. as read by
#'   \code{\link{readBin}} OR a character vector containing a URL or a path to a
#'   file on disk.
#' @param mesh3d Whether to return \code{rgl::mesh3d} object (when \code{TRUE},
#'   the default) or something as close as possible to what is provided by the
#'   Draco library (when \code{FALSE}).
#' @param ... Additional arguments passed to \code{\link{download.file}} when
#'   data is a URL (e.g. \code{quiet=TRUE} or \code{method})
#'
#' @return a \code{rgl::mesh3d} object or a list containing elements
#'   \code{points} and (for meshes). \code{faces}.
#'
#' @details Note that the Draco library returns 0-based indices for the faces
#'   whereas R in general and \code{rgl::mesh3d} in particular expect 1-based
#'   indices. When \code{mesh3d=FALSE}, the result will have 0-based indices as
#'   returned by the Draco library.
#'
#'   If \code{data} is an http/https URL it will be downloaded to a temporary
#'   location on disk (using \code{\link{download.file}}). If \code{data} is a
#'   character vector that does not look like a URL then it is assumed to refer
#'   to a file on disk (which will be read with \code{\link{readBin}}.
#'
#' @export
#' @importFrom utils download.file
#'
#' @examples
#' \donttest{
#' # fetch test data
#' # originally downloaded from:
#' carurl='https://github.com/google/draco/blob/master/testdata/car.drc?raw=true'
#' \dontrun{
#' car.m=draco_decode(carurl)
#' }
#' # use cached version in package for example
#' car.m=draco_decode(system.file('draco/car.drc', package = 'dracor'))
#' str(car.m)
#'
#' ## show the result
#' if(requireNamespace("rgl", quietly=TRUE)) {
#' rgl::shade3d(car.m, col='red')
#'
#' ## demonstrate conversion of raw form to rgl::mesh3d object
#' car.raw=draco_decode(carurl, mesh3d=FALSE)
#' str(car.raw)
#' car.m2 = rgl::tmesh3d(
#'   vertices = car.raw$points,
#'   indices = car.raw$faces + 1,
#'   homogeneous = FALSE)
#' }
#' }
draco_decode <- function(data, mesh3d=TRUE, ...) {
  if(is.character(data)) {
    path=data
    if(length(path)>1) stop("I can only read one file at a time!")
    if(isTRUE(grepl("^http[s]{0,1}://", path))) {
      # looks like a URL. Let's download
      tf=tempfile(fileext = '.drc')
      tryCatch(download.file(path, destfile = tf, ...),
               error=function(e) stop("Unable to download data from: ", path))
      path=tf
      on.exit(unlink(tf))
    }
    # let's read the file from disk
    tryCatch({
      data=readBin(path, what = raw(), n = file.info(path)$size)
      error=function(e) stop("Unable to read raw data from file: ", path)
    })
  }
  if(!is.raw(data))
    stop("The `data` argument should contain `raw` bytes!")
  if(!is.logical(mesh3d))
    stop("The `mesh3d` argument must be TRUE or FALSE!")
  if(length(data)==0)
    stop("The `data` argument contains 0 bytes!")

  # in raw mode do as little as possible to the output
  index_offset=ifelse(mesh3d, 1L, 0L)

  res=.Call(`_dracor_dracodecode`, data, index_offset)
  # errors return a list of length 1 containing a string
  if(length(res)==1 && is.character(res[[1]]))
    stop(res[[1]])

  if(isTRUE(mesh3d)) {
    res <- structure(list(vb = rbind(res$points, 1),
                          it = res$faces),
                     class = c("mesh3d", "shape3d"))
  }
  res
}
