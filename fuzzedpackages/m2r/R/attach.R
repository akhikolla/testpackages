.onAttach <- function(...) {

  packageStartupMessage('  Please cite m2r! See citation("m2r") for details.')

  # set gmp
  set_m2r_option(gmp = FALSE)

  # set pathing
  if (is.mac() || is.linux()) {
    set_m2r_option(m2_path = Sys.getenv("M2"))
  } else {
    set_m2r_option(m2_path = "")
  }


  # check that M2 was found in ~/.Renviron
  startup_check_for_program()


  # return
  invisible(TRUE)
}




.onDetach <- function(...) {
  stop_m2()
  options(m2r = NULL)
}
# restart R
# library(m2r)
# m2("1+1")
# getOption("m2r")
# detach("package:m2r")
# getOption("m2r")





psm  <- packageStartupMessage

psms <- function(fmt, ...) packageStartupMessage(sprintf(fmt, ...))


startup_check_for_program <- function(){

  if(get_m2_path() != ""){

    psms("  M2 found in %s", get_m2_path())
    return(invisible(FALSE))

  } else {

    psms("  M2 not found; defaulting to cloud.")
    psms("  Use set_m2_path(\"/path/to/m2\") to run M2 locally.")
    return(invisible(FALSE))

  }

  invisible(TRUE)
}



setOption <- function(optionName, value){
  eval(parse(text = sprintf('options("%s" = "%s")', optionName, value)))
}






# set_m2r_option both sets options for m2r in the list m2r in options
# and initialized the list when m2r is attached to the search path
# (search())
set_m2r_option <- function(...) {

  # if there is no m2r option (package is being initialized)
  # create the list with the arguments and return
  if ("m2r" %notin% names(options())) {
    options(m2r = list(...))
    return(invisible())
  }

  # otherwise, go through arguments sequentially and add/update
  # them in the list m2r in options
  m2r <- getOption("m2r")
  arg_list <- lapply(as.list(match.call())[-1], eval, envir = parent.frame())
  for (k in seq_along(arg_list)) {
    if (names(arg_list)[k] %in% names(m2r)) {
      m2r[names(arg_list)[k]] <- arg_list[k]
    } else {
      m2r <- c(m2r, arg_list[k])
    }
  }

  # set new m2r
  options(m2r = m2r)

  # return
  invisible()
}
# (l <- list(a = 1, b = 2, c = 3))
# l[d] <- 5
# l










