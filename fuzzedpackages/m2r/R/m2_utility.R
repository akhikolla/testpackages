#' Utility tools for M2
#'
#' Utility tools for M2
#'
#' @param x an object of class \code{m2}
#' @param m2_attr the name of an M2 attribute
#' @param name a string; the name of a M2 object
#' @param value the value to assign
#' @param m2_name  \code{m2_name}  M2 attribute
#' @param m2_class \code{m2_class} M2 attribute
#' @param base_class a base class; an R class to use for dispatching
#'   if there is no relevant method for the other classes (e.g.
#'   \code{m2})
#' @param m2_meta  \code{m2_meta}  M2 attribute
#' @param all.names if \code{TRUE}, all registered Macaulay2
#'   variables, including ones internally used by m2r, will be
#'   returned
#' @name m2_utility
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#' m2("a = 5")
#' m2_ls()
#' m2_exists("a")
#' m2("b = 1")
#' m2_exists(c("a","b","c"))
#'
#' m2_getwd()
#'
#' x <- 1
#' class(x) <- "m2"
#' attr(x, "m2_meta") <- list(a = 1, b = 2)
#' m2_meta(x)
#' m2_meta(x, "b")
#' m2_meta(x, "b") <- 5
#' m2_meta(x, "b")
#'
#' # R <- ring(c("x1", "x2", "x3"))
#' # m2_name(R)
#' # m2(sprintf("class %s", m2_name(R)))
#' # m2_ls()
#' # m2_rm(m2_name(R))
#' # m2_ls()
#' # m2(paste("class", m2_name(R)))
#'
#' m2_ls()
#' m2_ls(all.names = TRUE)
#'
#'
#' }



#' @rdname m2_utility
#' @export
m2_name <- function (x) {
  if ( is.m2(x) ) {
    attr(x, "m2_name")
  } else {
    character(0)
  }
}


#' @rdname m2_utility
#' @export
`m2_name<-` <- function (x, value) {
  stopifnot( is.m2(x) )
  attr(x, "m2_name") <- value
  x
}


#' @rdname m2_utility
#' @export
m2_meta <- function (x, m2_attr) {
  if ( !is.m2(x) ) return(NULL)
  if ( missing(m2_attr) ) return(attr(x, "m2_meta"))
  attr(x, "m2_meta")[[m2_attr]]
}


#' @rdname m2_utility
#' @export
`m2_meta<-` <- function (x, m2_attr, value) {
  stopifnot( is.m2(x) )
  if (missing(m2_attr)) {
    attr(x, "m2_meta") <- value
  } else {
    meta <- m2_meta(x)
    meta[[m2_attr]] <- value
    attr(x, "m2_meta") <- meta
  }
  x
}


#' @rdname m2_utility
#' @export
m2_structure <- function (x = NA, m2_name, m2_class, m2_meta, base_class) {

  if (!missing(m2_class)) class(x) <- c(m2_class, "m2")
  if (!missing(base_class)) class(x) <- c(class(x), base_class)
  if (!missing(m2_name)) m2_name(x) <- m2_name
  # if (m2_meta(x) != NULL) m2_meta(x)
  if (!missing(m2_meta)) m2_meta(x) <- m2_meta

  x
}


#' @rdname m2_utility
#' @export
m2_exists <- function(name) {
  if(!is.character(name)) name <- deparse(substitute(name))
  name %in% m2_ls()
}


#' @rdname m2_utility
#' @export
m2_ls <- function(all.names = FALSE) {
  out <- m2("userSymbols()")
  out <- str_sub(out, 2, -2)
  out <- str_split(out, ",")[[1]]

  # "symbols m2o1" -> "m2o1"
  out <- str_sub(out, 8)

  # remove internals and m2o#'s
  if(!all.names) {
    out <- out[!str_detect(out, "m2rint")]
    out <- out[!str_detect(out, "m2o[0-9]+")]
  }

  # return
  out
}


#' @rdname m2_utility
#' @export
m2_rm <- function(name) {
  stop("broken.")
  if (!is.m2(name)) return(invisible())
  m2(paste(m2_name(name), "=symbol", m2_name(name)))
  invisible()
}


#' @rdname m2_utility
#' @export
m2_getwd <- function() {
  str_sub(m2("currentDirectory()"), 2, -2)
}
