#' Enter a Macaulay2 session
#'
#' Enter a Macaulay2 session
#'
#' @param port port for Macaulay2 socket
#' @param timeout number of seconds before aborting
#' @return \code{TRUE} invisibly
#' @export
#' @examples
#'
#' \dontrun{ requires Macaulay2 be installed and an interactive session
#'
#' enter_m2()
#'
#' # m2 code below
#' 1 + 1
#' a = 1
#' a
#' R = QQ[t,x,y,z]
#' I = ideal(t^4  -  x, t^3  -  y, t^2  -  z)
#' gens gb I
#' exit
#'
#' # back in R, the variable persists using m2()
#' m2("a")
#' m2("I")
#'
#'
#' # we can also define variables in R that persist in m2
#' m2("b = 5")
#'
#' enter_m2()
#' b
#' exit
#'
#' }
#'
enter_m2 <- function (port = 27436L, timeout = 10) {
  if(!interactive()) stop("enter_m2() is only available in interactive sessions.")
  stop("enter_m2() is currently not working, check back soon!")
  suppressMessages( start_m2(port, timeout) )
  message("Entering M2 mode. Type 'exit' to go back into R.")
  repeat {
    i <- readline("i : ")
    if(i == "exit") break
    o <- m2(i)
    cat(paste0("o : ", o))
  }
  invisible(TRUE)
}
