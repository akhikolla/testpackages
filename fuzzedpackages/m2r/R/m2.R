#' Call and reset a Macaulay2 process
#'
#' Call and reset a Macaulay2 process
#'
#' @param port port for Macaulay2 socket
#' @param timeout number of seconds before aborting
#' @param attempts numer of times to try to make connection
#' @param cloud use a cloud?
#' @param hostname the remote host to connect to; defaults to the Amazon EC2
#'   instance
#' @param code Macaulay2 code
#' @param x formal argument for print method
#' @param ... ...
#' @return m2 return value
#' @name m2_call
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#' m2("1 + 1")
#' m2.("1 + 1")
#'
#' m2("factor 32004")
#'
#' # run a chunk of m2 code, only pulling the end value back into R
#' m2("
#'   R = QQ[a..d]
#'   I = ideal(a^3-b^2*c, b*c^2-c*d^2, c^3)
#'   G = gens gb I
#' ")
#'
#' # illustrate the persistent connection
#' m2("a = 1 + 1")
#' m2("a")
#' reset_m2()
#' m2("a")
#'
#'
#' # forcing a cloud start
#' if(has_m2_connection()) stop_m2()
#' start_m2(cloud = TRUE)
#' m2("1 + 1")
#' stop_m2()
#'
#'
#'
#' m2.("peek(QQ[x,y,z])")
#' m2("peek(QQ[x,y,z])")
#'
#' # m2 returns in its ext_str position the result of running
#' # toExternalString on the return value of the chunk of code
#' # you run. in principle, toExternalString provides the code
#' # needed to recreate the m2 object of interest. however,
#' # does not work for all objects representable in the m2 language.
#' # in particular, mutable objects are not supported.
#' # this is what happens when you look at those:
#' m2.("new MutableList from {1,2,3}")
#' m2("new MutableList from {1,2,3}")
#'
#' }


#' @export
#' @rdname m2_call
m2r_version_number <- function() {
  "1.0.0"
}



#' @export
#' @rdname m2_call
m2r_cloud_url <- function() {
  "ec2-52-10-66-241.us-west-2.compute.amazonaws.com"
}


#' @export
#' @rdname m2_call
has_m2_connection <- function() {

  if (is.null(get_m2_connection())) {
    FALSE
  } else {
    TRUE
  }

}



#' @export
#' @rdname m2_call
start_m2 <- function(
  port = 27436L, timeout = 10, attempts = 10,
  cloud = FALSE, hostname = m2r_cloud_url()
) {

  # don't increment port if supplied
  if (!missing(port)) attempts <- 1

  # if already running M2, break
  if (has_m2_connection()) return(invisible(0L))

  if(!is.null(get_m2_path()) && !cloud && missing(hostname)) { # m2 found locally

    if (do_start_m2_local(port = port, timeout = timeout, attempts = attempts) != 0L)
      stop("m2r unable to connect to local instance")

  } else { # default to cloud

    cloud_out <- do_start_m2_cloud(hostname)
    if(cloud_out == 1L) stop("m2r unable to connect to the cloud instance, are you online?")

  }

}



do_start_m2_cloud <- function(hostname = m2r_cloud_url()) {

  # if already running M2, break
  if (has_m2_connection()) return(invisible(0L))

  # launch M2 on cloud
  message("Connecting to M2 in the cloud... ", appendLF = FALSE)

  # request port
  port <- request_port(hostname = hostname, port = 27435L)

  # connect to local server
  out <- connect_to_m2_server(hostname = hostname, port = port)
  message("done.")

  # return
  invisible(out)
}



do_start_m2_local <- function(port = 27436L, timeout = 10, attempts = 10) {

  # if already running M2, break
  if (has_m2_connection()) return(invisible(0L))

  # find the first open port
  openPortFound <- FALSE
  for (i in seq.int(0, attempts-1)) {

    out <- system(paste0("netstat -an | grep '[\\.:]", port, "' | wc -l"), intern = TRUE)
    if (as.integer(out) == 0) break()
    # openPortFound <- TRUE
    # tempservercon <- NULL
    # tryCatch(
    #   tempservercon <- suppressWarnings(
    #     socketConnection(
    #       port = port, blocking = FALSE,
    #       server = TRUE, open = "r+", timeout = 1
    #     )
    #   ),
    #   error = function(e) { openPortFound <- FALSE }
    # )
    # if (!is.null(tempservercon)) close(tempservercon)
    # if (openPortFound) break()

    if (i == attempts - 1) {
      message(sprintf("%s attempts made at finding an open port. Aborting start_m2", attempts))
      return(invisible(1L))
    } else {
      # message(sprintf("Unable to connect to M2 on port %s. Attempting to connect on port %s", port, port + 1))
      port <- port + 1
    }

  }

  # launch M2 on local server
  message("Starting M2... ", appendLF = FALSE)
  if(is.mac() || is.unix()) {
    system2(
      file.path2(get_m2_path(), "M2"), args = c("--script", system.file("server", "m2rserverscript.m2", package = "m2r"), toString(port)),
      stdout = NULL, stderr = NULL, stdin = "",
      wait = FALSE
    )
  } else if(is.win()) {
    stop("Running local instances of M2 is not yet supported.")
  }

  # connect to local server
  out <- connect_to_m2_server(port = port, timeout = timeout)
  message("done.")

  # post process id to m2r global options
  if (out == 0) set_m2r_option(m2_procid = strtoi(m2("processID()")))

  # return
  out
}







request_port <- function(
  hostname = m2r_cloud_url(),
  port = 27435L, timeout = 10
) {

  # initialize client socket
  con <- NULL
  for (i in 0:(20*timeout)) {

    tryCatch(
      con <- suppressWarnings(
        socketConnection(
          host = hostname, port = port,
          blocking = FALSE, server = FALSE,
          open = "r+", timeout = 60*60*24*7
        )
      ),
      error = function(e) {  }
    )

    if (!is.null(con)) { break } else { Sys.sleep(0.05) }

  }

  if (is.null(con)) stop("m2r unable to connect to the cloud instance, are you online?")

  repeat {
    # read output info
    port_number <- readLines(con, 1)

    if (length(port_number) > 0) break

    i <- i + 1
    if (timeout > 0 && i >= timeout * 2000) {
      break
    } else {
      Sys.sleep(0.0005)
    }
  }

  close(con)

  if (length(port_number) == 0 || port_number == "0") {
    stop(sprintf("Macaulay2 cloud is full; please try again later."))
  }

  return(strtoi(port_number))

}



connect_to_m2_server <- function(hostname = "localhost", port = 27436L, timeout = 10) {

  # initialize client socket
  con <- NULL
  for (i in 0:(20*timeout)) {

    tryCatch(
      con <- suppressWarnings(
        socketConnection(
          host = hostname, port = port,
          blocking = FALSE, server = FALSE,
          open = "r+", timeout = 60*60*24*7
        )
      ),
      error = function(e) {  }
    )

    if (!is.null(con)) {
      break
    } else {
      Sys.sleep(0.05)
    }

  }

  if (is.null(con)) return(invisible(1L))

  repeat {
    # read output info
    server_version <- readLines(con, 1)

    if (length(server_version) > 0) break

    i <- i + 1
    if (timeout > 0 && i >= 2000*timeout) {
      break
    } else {
      Sys.sleep(0.0005)
    }
  }

  if (server_version != m2r_version_number()) {
    close(con)
    stop(sprintf("Internal error: server version is %s and client version is %s.",
                 server_version, m2r_version_number()))
  }

  # set options
  set_m2r_option(
    m2_con = con,
    m2_port = port,
    m2_timeout = timeout
  )

  set_m2r_option(m2_procid = -1L)

  invisible(0L)
}





#' @export
#' @rdname m2_call
stop_m2 <- function() {

  if (has_m2_connection()) { # for detaching when m2 never run

    # send kill code
    writeLines("", get_m2_connection())
    close(get_m2_connection())

    # not elegant, but a necessary safety measure
    Sys.sleep(0.01)
    if (get_m2_procid() >= 0) tools::pskill(get_m2_procid())

    set_m2r_option(m2_con = NULL, m2_procid = NULL)

  }

}




#' @export
#' @rdname m2_call
reset_m2 <- function(
  port = 27436L, timeout = 10, attempts = 10,
  hostname = "ec2-52-10-66-241.us-west-2.compute.amazonaws.com"
) {

  stop_m2()
  start_m2(port, timeout, attempts, hostname)

}









#' @export
#' @rdname m2_call
m2 <- function(code, timeout = -1) {
  m2_meta(do.call(m2., as.list(match.call())[-1]), "ext_str")
}





#' @export
#' @rdname m2_call
m2. <- function(code, timeout = -1) {

  # ensure m2 is running
  start_m2()

  # preempt m2 kill code
  if (code == "") return("")

  # write to connection
  writeLines(code, get_m2_connection())

  i <- 0
  outinfo <- NULL
  repeat {
    # read output info
    outinfo <- readLines(get_m2_connection(), 1)

    if (length(outinfo) > 0) break

    i <- i + 1
    if (timeout > 0 && i >= timeout * 2000) {
      break
    } else {
      Sys.sleep(0.0005)
    }
  }

  if (length(outinfo) > 0) {
    # read output from connection and return
    info <- strsplit(outinfo, " ", fixed = TRUE)[[1]]

    retcode <- strtoi(info[1])
    numlines <- strtoi(info[2])
    m2_name <- info[3]
    m2_class <- info[4]
    m2_class_class <- info[5]
  } else {
    # cancel command if needed
    tools::pskill(get_m2_procid(), tools::SIGINT)
    Sys.sleep(0.01)

    retcode <- -1L
    numlines <- -1L
  }

  output <- paste(readLines(get_m2_connection(), numlines), collapse = "\n")

  if (retcode == -1L) {
    # timeout occurred, kill M2 instance and stop
    stop_m2()
    # start_m2(getOption("m2_port"), getOption("m2_timeout"))
    stop("Command timed out, M2 connection lost")
  } else if (retcode == 1L) {
    # user's code string errored, alert them
    stop(output, call. = FALSE)
  } else if (retcode == 2L) {
    # toExternalString failed, make ext_str NULL
    output <- NULL
  }

  # assemble pointer object and return
  m2_structure(
    m2_name = m2_name,
    m2_class = "m2_pointer",
    m2_meta = list(
      ext_str = output,
      m2_class = m2_class,
      m2_class_class = m2_class_class
    )
  )
}






#' @export
#' @rdname m2_call
print.m2_pointer <- function (x, ...) {
  cat("M2 Pointer Object\n")
  if(is.null(m2_meta(x, "ext_str"))) m2_meta(x, "ext_str") <- ""
  w <- min(c(options()$width, 80), na.rm = TRUE) - 19
  if(nchar(m2_meta(x, "ext_str")) > w) {
    ext_str <- str_c(str_sub(m2_meta(x, "ext_str"), 1, w-4), "...")
  } else {
    ext_str <- m2_meta(x, "ext_str")
  }
  cat(sprintf("  ExternalString : %s\n", ext_str))
  # cat(sprintf("          R Name : %s\n", deparse(substitute(x))))
  cat(sprintf("         M2 Name : %s\n", m2_name(x)))
  cat(sprintf("        M2 Class : %s (%s)\n", m2_meta(x, "m2_class"), m2_meta(x, "m2_class_class")))
  invisible(x)
}



