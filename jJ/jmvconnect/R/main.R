
evalRemote <- function(script, dataPath, echo, outputPath, ...) {
    options <- list(...)
    figWidth <- options$figWidth
    figHeight <- options$figHeight
    data <- readDF(dataPath, NULL, FALSE)
    results <- eval(script, data, echo, figWidth=figWidth, figHeight=figHeight)
    saveRDS(results, file=outputPath)
    NULL
}

#' @useDynLib jmvconnect
#' @importFrom Rcpp evalCpp
read <- function(path, columns=NULL, headerOnly=FALSE) {
    readDF(path, columns, headerOnly)
}

#' @importFrom rappdirs user_config_dir
#' @importFrom httr GET content
listDS <- function() {

    if (Sys.info()['sysname'] == 'Linux')
        appDataDir <- '~/.jamovi'
    else
        appDataDir <- file.path(user_config_dir(), 'jamovi')

    portFile <- dir(appDataDir, '[0-9]+\\.port')
    if (length(portFile) < 1)
        stop('A running jamovi instance could not be found', call.=FALSE)

    match <- regexec('([0-9]+)\\.port', portFile[1])[[1]]
    if (match[1] == -1)
        stop('A running jamovi instance could not be found', call.=FALSE)

    port <- substring(portFile, match[2], attr(match, 'match.length')[2])
    url <- paste0('http://127.0.0.1:', port, '/api/datasets')

    response <- GET(url)
    content <- content(response)
    content
}

#' Reads a data set from jamovi
#'
#' @param id the number, or the title of the data set to read
#' @param columns (optional) only reads the columns named
#' @return the data set as a data frame
#'
#' @examples
#' \dontrun{
#' jmvconnect::what()
#'
#' #  Available Data Sets
#' #  -------------------------------------
#' #         Title           Rows    Cols
#' #  -------------------------------------
#' #    1    iris             150       5
#' #    2    Tooth Growth      60       3
#' #  -------------------------------------
#'
#' data <- jmvconnect::read('Tooth Growth')
#'
#' # or
#'
#' data <- jmvconnect::read(2)
#' }
#' @export
read <- function(id, columns) {

    if (missing(columns))
        columns <- NULL

    datasets <- listDS()

    if (is.numeric(id)) {
        if (id > length(datasets))
            stop('No such data set')
        dataset <- datasets[[id]]
    } else {
        for (i in seq_along(datasets)) {
            dataset <- datasets[[i]]
            if (dataset$title == id)
                break()
            dataset <- NULL
        }
    }

    if (is.null(dataset))
        stop('No such data set')

    readDF(dataset$buffer, columns, FALSE)
}

#' Lists the data sets available from jamovi
#'
#' Lists the data sets available from jamovi. Data sets can then be read using
#' the read() function.
#'
#' @examples
#' \dontrun{
#' jmvconnect::what()
#'
#' #  Available Data Sets
#' #  -------------------------------------
#' #         Title           Rows    Cols
#' #  -------------------------------------
#' #    1    iris             150       5
#' #    2    Tooth Growth      60       3
#' #  -------------------------------------
#'
#' data <- jmvconnect::read('Tooth Growth')
#'
#' # or
#'
#' data <- jmvconnect::read(2)
#' }
#'
#' @export
what <- function() {

    content <- listDS()

    table <- jmvcore::Table$new(
        title='Available Data Sets')
    table$addColumn(name='id', title='', content='($key)')
    table$addColumn(name='title', title='Title', type='Number')
    table$addColumn(name='rows', title='Rows')
    table$addColumn(name='cols', title='Cols')

    for (i in seq_along(content)) {
        dataset <- content[[i]]
        table$addRow(rowKey=i, values=list(
            title=dataset$title,
            rows=dataset$rowCount,
            cols=dataset$columnCount))
    }

    table
}
