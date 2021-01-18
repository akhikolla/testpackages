
#' @importFrom utils capture.output
#' @importFrom evaluate evaluate new_output_handler
#' @importFrom jmvcore Options Group Image Preformatted
eval <- function(script, data, echo, root, ...) {

    eval.env <- new.env()

    if ( ! missing(data))
        eval.env$data <- data

    env <- new.env()
    env$count <- 1

    conf <- list(...)

    figWidth <- as.integer(conf$figWidth)
    if (length(figWidth) == 1 && ! is.na(figWidth))
        env$figWidth <- figWidth
    else
        env$figWidth <- 400

    figHeight <- as.integer(conf$figHeight)
    if (length(figHeight) == 1 && ! is.na(figHeight))
        env$figHeight <- figHeight
    else
        env$figHeight <- 300

    env$echo <- isTRUE(echo)

    options <- jmvcore::Options$new()

    if (missing(root))
        root <- Group$new(options, title="Results")

    text_handler <- function(object, capture) {

        if (inherits(object, 'ResultsElement')) {

            object$print()

        } else {

            results <- jmvcore::Preformatted$new(options, paste(env$count))
            env$count <- env$count + 1
            env$last <- NULL
            root$add(results)

            if (is.character(object) && ! capture) {
                value <- object
            }
            else {
                value <- capture.output(object)
            }

            results$setContent(value)
        }

        object
    }

    source_handler <- function(value) {
        if ( ! env$echo)
            return()

        value <- trimws(value$src)
        if (value == '')
            return()

        if (is.null(env$last) || ! inherits(env$last, 'Preformatted')) {
            results <- Preformatted$new(options, paste(env$count))
            root$add(results)
            env$count <- env$count + 1
        }
        else {
            results <- env$last
        }

        value <- paste0('> ', value)

        content <- results$content
        if (content != '')
            content <- paste0(content, '\n', value)
        else
            content <- value

        results$setContent(content)
        env$last <- results
    }

    graphics_handler <- function(plot) {
        results <- Image$new(
            options=options,
            name=paste(env$count),
            renderFun='.render',
            width=env$figWidth,
            height=env$figHeight)
        root$add(results)
        results$setState(plot)
        env$count <- env$count + 1
        env$last <- NULL
    }

    handler <- new_output_handler(
        source=source_handler,
        text=function(text) text_handler(text, FALSE),
        value=function(text) text_handler(text, TRUE),
        graphics=graphics_handler)

    evaluate(
        input=script,
        envir=eval.env,
        output_handler=handler,
        stop_on_error=2)

    root
}
