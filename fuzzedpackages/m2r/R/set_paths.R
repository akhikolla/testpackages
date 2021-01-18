#' Set path to Macaulay2 (M2)
#'
#' These are helper functions that deal with pathing to Macaulay2 and asking if it
#' is present. When the Macaulay2 package is loaded it attempts to find the
#' Macaulay2 executable by looking for an environment variable indicating where it
#' is, i.e. its path as specified in your .Renviron file.
#'
#' For easiest use, you'll want to specify the path the Macaulay2 executable in
#' your ~/.Renviron file. It should look something like
#'
#' \code{M2=/Applications/Macaulay2-1.10/bin}
#'
#' You can set this permanently with [edit_r_environ()]. Note that absolute
#' paths should be specified, not relative paths, e.g. don't use ~/path/to/exe.
#'
#' You can change this for the current session using [set_m2_path()], which
#' accepts a character string or, if missing, uses [file.choose()] to let you
#' interactively; you just select an arbitrary executable.
#'
#' On Windows, m2r just defaults to the cloud implementation. Local M2 instances
#' are not currently supported on Windows.
#'
#' @param path A character string, the path to M2
#' @return An invisible character string, the path found.  More importantly, the
#'   function has the side effect of setting the global m2r option "m2_path"
#' @export
#' @name m2_path
#' @author David Kahle \email{david@@kahle.com}
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#'
#' getOption("m2r")
#' get_m2_path()
#' set_m2_path()
#'
#'
#' ## each of these functions can be used statically as well
#' (m2_path <- get_m2_path())
#' set_m2_path("/path/to/m2/directory")
#' get_m2_path()
#' set_m2_path(m2_path) # undoes example
#'
#'
#' # if you'd like to use the cloud, after you library(m2r)
#' # and before you use m2() type
#' set_m2_path(NULL)
#'
#' # alternatively, if you have already been using m2, do:
#' stop_m2()
#' set_m2_path(NULL)
#' m2("1+1")
#'
#'
#' }






#' @rdname m2_path
#' @export
set_m2_path <- function(path = NULL){

  if(missing(path) && interactive()){

    path <- dirname(file.choose())
    if(is.win() && str_detect(path,"C:/")){
      mpath <- str_replace(dirname(path), "C:/", "/cygdrive/c/")
    }
    set_m2r_option(m2_path = path)
    return(invisible(path))

  } else if (!missing(path)) {

    set_m2r_option(m2_path = path)
    return(invisible(path))

  } else {
    stop(
      "If the session is not interactive, a path must be specified.",
      call. = FALSE
    )
  }
}





#' @rdname m2_path
#' @export
get_m2_path <- function() getOption("m2r")$m2_path



#' @rdname m2_path
#' @export
get_m2_connection <- function() getOption("m2r")$m2_con


#' @rdname m2_path
#' @export
get_m2_con <- function() {
  .Deprecated("get_m2_connection")
  getOption("m2r")$m2_con
}



#' @rdname m2_path
#' @export
get_m2_procid <- function() getOption("m2r")$m2_procid



#' @rdname m2_path
#' @export
get_m2_port <- function() getOption("m2r")$m2_port



#' @importFrom usethis edit_r_environ
#' @export
usethis::edit_r_environ



