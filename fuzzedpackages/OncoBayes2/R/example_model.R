#' Runs example models
#'
#' @param topic example to run
#' @param envir environment which the example is loaded into. Defaults
#'     to the caller environment.
#' @param silent logical controlling if execution is run silently
#'     (defaults to \code{FALSE})
#'
#' @return When topic is not specified a list of all possible topics
#'     is return. Whenever a valid topic is specified, the function
#'     inserts the example into the environment given and returns
#'     (invisibly) the updated environment.
#'
#' @template start-example
#' @examples
#'
#' ## get a list of available examples
#' example_model()
#'
#' ## run 3 component example
#' example_model("combo3")
#'
#' @template stop-example
#'
#' @export
example_model <- function(topic, envir=parent.frame(), silent=FALSE) {
    if(missing(topic))
        return(names(example_cache))
    assert_character(topic)
    assert_that(topic %in% names(example_cache), msg="Unkown example. For a list of examples call example_model().")
        ex_str <- example_cache[[topic]]
    if(silent) {
        suppressMessages(capture.output(eval(parse(text=ex_str), envir=envir)))
    } else {
        message("Running ", topic, " example:\n")
        message(paste(c(ex_str, ""), collapse="\n"))
        eval(parse(text=ex_str), envir=envir)
    }
    invisible(envir)
}

