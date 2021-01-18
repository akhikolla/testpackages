#' Descriptives of reaction time data
#'
#' @param formula A formula object of the form: \code{binary response ~ reaction time + condition1 * condition2}
#' @param data A dataframe for looking up data specified in formula.
#' For backwards compatibility this can also be with: a column named \code{rt} containing response times in ms,
#' a column named \code{response} containing at most 2 response options, and an
#' optional column named \code{condition} containing a numeric index as to which conditions
#' observations belong.
#' @param plot Logical, should a density plot of all condition-response pairs be made?
#' @param verbose Logical, should a table of counts and proportions be printed?
#'
#' @return Invisibly returns an object of class 'D*M'. It's first element is \code{table} and contains raw counts and proportions for
#' condition response pairs, conditions, and responses. It's second element \code{plot} contains a ggplot object.
#'
#' @details This function and \code{\link{rtHist}} are helper functions to inspect raw data.
#'
#' @examples
#' tt <- seq(0, 5, .01)
#' pars <- matrix(.5, 5, 2)
#' pars[1, ] <- 1
#' pars[2, ] <- c(0, 2)
#' dat <- simData(n = 3e3, pars = pars, tt = tt, pdfND = dbeta(tt, 10, 30))
#' x <- rtDescriptives(data = dat)
#'
#' print(x$table, what = 'cr')
#' print(x$table, what = 'c')
#' print(x$table, what = 'r')



#' @export
rtDescriptives <- function(formula = NULL, data, plot = TRUE, verbose = TRUE) {

  data <- getData(formula, data)
  rtime <- data[["rtime"]]
  response <- data[["response"]]
  condition <- data[["condition"]]
  hasConditions <- data[["hasConditions"]]
  data <- data[["data"]]

  lenCR <- tapply(data[[rtime]], list(data[[condition]], data[[response]]),
    length)
  d <- dim(lenCR)
  lenC <- .rowSums(lenCR, m = d[1], n = d[2])
  lenR <- .colSums(lenCR, m = d[1], n = d[2])

  table <- list(counts = list(conditionResponse = lenCR, condition = lenC,
    response = lenR), props = list(conditionResponse = lenCR/lenC, condition = lenC/sum(lenC),
    response = lenR/sum(lenR)), responses = colnames(lenCR))
  class(table) <- "DstarM"
  if (verbose)
    print(table)

  if (hasConditions) {
    tmp <- sprintf("interaction(%s, %s)", response, condition)
    mapping <- ggplot2::aes_string(x = rtime, group = tmp, fill = tmp)
  } else {
    mapping <- ggplot2::aes_string(x = rtime, group = response, fill = response)
  }
  graph <- ggplot2::ggplot(data = data, mapping = mapping) + ggplot2::geom_density(alpha = 0.25)

  if (plot)
    print(graph)

  out <- list(table = table, graph = graph)
  class(out) <- "DstarM.Descriptives"
  return(invisible(out))

}





