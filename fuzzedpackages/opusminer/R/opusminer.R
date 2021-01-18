# PACKAGE_DOCUMENTATION=========================================================

#' Filtered Top-k Association Discovery of Self-Sufficient Itemsets
#'
#' The \code{opusminer} package provides an R interface to the OPUS Miner
#' algorithm (implemented in C++), developed by Professor Geoffrey I Webb, for
#' finding the top \emph{k}, non-redundant itemsets on the measure of interest.
#'
#' OPUS Miner is a branch-and-bound algorithm for efficient discovery of
#' self-sufficient itemsets. For a user-specified \emph{k} and interest measure,
#' OPUS Miner finds the top \emph{k} productive non-redundant itemsets with
#' respect to the specified measure. It is then straightforward to filter out
#' those that are not independently productive with respect to that set,
#' resulting in a set of self-sufficient itemsets.
#'
#' OPUS Miner is based on the OPUS search algorithm.  OPUS is a set enumeration
#' algorithm distinguished by a computationally efficient pruning mechanism that
#' ensures that whenever an item is pruned, it is removed from the entire search
#' space below the parent node.
#'
#' OPUS Miner systematically traverses viable regions of the search space (using
#' depth-first search), maintaining a collection of the top \emph{k} productive
#' non-redundant itemsets in the search space explored. When all of the viable
#' regions have been explored, the top \emph{k} productive non-redundant
#' itemsets in the search space explored must be the top \emph{k} for the entire
#' search space.
#'
#' A comprehensive explanation of the algorithm is provided in the article cited
#' below.
#'
#' @references
#' Webb, G. I., & Vreeken, J. (2014). Efficient Discovery of the Most
#' Interesting Associations. \emph{ACM Transactions on Knowledge Discovery from
#' Data}, 8(3), 1-15. doi: http://dx.doi.org/10.1145/2601433
#'
#' @docType package
#' @name opusminer-package
NULL

# EXTERNAL_FUNCTIONS============================================================

#' @title
#' Filtered Top-k Association Discovery of Self-Sufficient Itemsets
#'
#' @description
#' \code{opus} finds the top \emph{k} productive, non-redundant itemsets on the
#' measure of interest (leverage or lift) using the OPUS Miner algorithm.
#'
#' @details
#' \code{opus} provides an interface to the OPUS Miner algorithm (implemented in
#' C++) to find the top \emph{k} productive, non-redundant itemsets by leverage
#' (default) or lift.
#'
#' \code{transactions} should be a filename, list (of transactions, each list
#' element being a vector of character values representing item labels), or an
#' object of class \code{\link[arules]{transactions}} (\code{arules}).
#'
#' Files should be in the format of a list of transactions, one line per
#' transaction, each transaction (ie, line) being a sequence of item labels,
#' separated by the character specified by the parameter \code{sep} (default "
#' ").  See, for example, the files at \url{http://fimi.ua.ac.be/data/}.
#' (Alternatively, files can be read seaparately using the
#' \code{\link{read_transactions}} function.)
#'
#' \code{format} should be specified as either "data.frame" (the default) or
#' "itemsets", and any other value will return a list.
#'
#' @references
#' Webb, G. I., & Vreeken, J. (2014). Efficient Discovery of the Most
#' Interesting Associations. \emph{ACM Transactions on Knowledge Discovery from
#' Data}, 8(3), 1-15. doi: http://dx.doi.org/10.1145/2601433
#'
#' @param transactions A filename, list, or object of class
#'   \code{\link[arules]{transactions}} (\code{arules}).
#' @param k The number of itemsets to return, an integer (default 100).
#' @param format The output format ("data.frame", default, or "itemsets").
#' @param sep The separator between items (for files, default " ").
#' @param print_closures return the closure for each itemset (default \code{FALSE})
#' @param filter_itemsets filter itemsets that are not independently productive (default \code{TRUE})
#' @param search_by_lift make lift (rather than leverage) the measure of interest (default \code{FALSE})
#' @param correct_for_mult_compare correct alpha for the size of the search space (default \code{TRUE})
#' @param redundancy_tests exclude redundant itemsets (default \code{TRUE})
#'
#' @return  The top \emph{k} productive, non-redundant itemsets, with relevant
#'   statistics, in the form of a data frame, object of class
#'   \code{\link[arules]{itemsets}} (\code{arules}), or a list.
#'
#' @examples
#' \dontrun{
#'
#' result <- opus("mushroom.dat")
#' result <- opus("mushroom.dat", k = 50)
#'
#' result[result$self_sufficient, ]
#' result[order(result$count, decreasing = TRUE), ]
#'
#' trans <- read_transactions("mushroom.dat", format = "transactions")
#'
#' result <- opus(trans, print_closures = TRUE)
#' result <- opus(trans, format = "itemsets")
#' }
opus <- function(transactions,
                 k = 100,
                 format = "data.frame",
                 sep = " ",
                 print_closures = FALSE,
                 filter_itemsets = TRUE,
                 search_by_lift = FALSE,
                 correct_for_mult_compare = TRUE,
                 redundancy_tests = TRUE) {

  # initialise output
  output <- NULL

  # check arguments
  if (k < 1) {k <- 1}
  cpp_arguments <- .check_cpp_arguments(as.list(environment()))

  if (.valid_input(transactions)) {

    time <- rep(0, 4)
    time[1] <- proc.time()[3]

    cat("Reading file/data...\n\n")

    # read and index the input data according to its format
    if (.valid_filename(transactions)) {
      input <- .encode(.read_transactions(transactions, sep = sep))
    } else if (is(transactions, "list")) {
      input <- .encode(transactions)
    } else if (is(transactions, "transactions")) {
      input <- .encode(as(transactions, "list"))
    }

    time[2] <- proc.time()[3]

    # call OPUS Miner (C++)
    output <- .opus_cpp(input$tidlist,
                        input$num_items,
                        input$num_trans,
                        k,
                        cpp_arguments)

    time[3] <- proc.time()[3]

    # decode (deindex) the itemsets
    output$itemset <- .decode(output$itemset, input$index)

    # if relevant, decode (deindex) the closure itemsets
    if (cpp_arguments["print_closures"] == TRUE) {
      output$closure <- .decode(output$closure, input$index)
    } else {
      output$closure <- NULL
    }

    # time[4] <- proc.time()[3]

    # format the output
    if (format == "data.frame") {
      output <- .as_data_frame(output)
    } else if (format == "itemsets") {
      output <- .as_itemsets(output)
    }

    cat("[[",
        round(time[3] - time[1], 2), " seconds: ",
        round(time[2] - time[1], 2), " read, ",
        round(time[3] - time[2], 2), " find",
        ifelse(filter_itemsets, " & filter]]", "]]"),
        sep = "")

  } else if (!is.null(transactions)) {
    message("invalid filename or transaction data")
  } else {
    message("no filename or transaction data specified")
  }

  return(output)
}

#' @title Read Transaction Data from a File (Fast)
#'
#' @description
#' \code{read_transactions} reads transaction data from a file fast, providing a
#' significant speed increase over alternative methods for larger files.
#'
#' @details
#' \code{read_transactions} uses (internally) the \code{\link[base]{readChar}}
#' function to read transaction data from a file fast.  This is substantially
#' faster for larger files than alternative methods.
#'
#' Files should be in the format of a list of transactions, one line per
#' transaction, each transaction (ie, line) being a sequence of item labels,
#' separated by the character specified by the parameter \code{sep} (default "
#' ").  See, for example, the files at \url{http://fimi.ua.ac.be/data/}.
#'
#' @param filename A filename.
#' @param sep The separator between items (default " ").
#' @param format The output format ("list" or "transactions").
#'
#' @return The transaction data, in the form of a list (of transactions, each
#' list element being a vector of character values representing item labels), or
#' an object of class \code{\link[arules]{transactions}} (\code{arules}).
#'
#' @examples
#' \dontrun{
#'
#' trans <- read_transactions("mushroom.dat")
#' trans <- read_transactions("mushroom.dat", format = "transactions")
#' }
read_transactions <- function(filename, sep = " ", format = "list") {
  return(.read_transactions(filename, sep, format))
}

# INTERNAL_FUNCTIONS============================================================

# return output as a data frame
.as_data_frame <- function(output) {

  # "transpose" list elements to data frame columns
  output <- as.data.frame(sapply(output, "["))

  # set data frame coloumn types
  output$itemset <- sapply(output$itemset, paste, collapse = ", ")
  output$count <- as.integer(output$count)
  output$value <- as.numeric(output$value)
  output$p <- as.numeric(output$p)
  output$self_sufficient <- as.logical(output$self_sufficient)
  if ("closure" %in% names(output)) {
    output$closure <- sapply(output$closure, paste, collapse = ", ")
  }

  return(output)
}

# return output as an object of class itemsets (arules)
.as_itemsets <- function(output) {
  return(
    new("itemsets",
        items = as(.as_transactions(output$itemset), "itemMatrix"),
        quality = .as_data_frame(output)[-1])
  )
}

# return a list of transactions as an object of class ngCMatrix (Matrix)
.as_sparse_matrix <- function(transactions) {
  index <- unique(unlist(transactions))
  sm_row_index <- match(unlist(transactions), index)
  sm_col_index <- rep(seq_along(transactions), lengths(transactions))
  return(Matrix::sparseMatrix(i = sm_row_index, j = sm_col_index))
}

# return a list of transactions as an object of class transactions (arules)
.as_transactions <- function(transactions) {
  return(
    new("transactions",
        data = .as_sparse_matrix(transactions),
        itemInfo = data.frame(labels = unique(unlist(transactions)),
                              stringsAsFactors = FALSE))
  )
}

# validate cpp arguments
.check_cpp_arguments <- function(arg) {

  # default cpp arguments
  def <- list(print_closures = FALSE,
              filter_itemsets = TRUE,
              search_by_lift = FALSE,
              correct_for_mult_compare = TRUE,
              redundancy_tests = TRUE)

  if (length(arg) > 0) {
    # drop arguments not scalar boolean
    arg <- arg[sapply(arg, function(e){length(e) == 1 && is.logical(e)})]
    # replace:
    #   - elements of def with names matching elements of arg; with
    #   - elements of arg with names matching elements of def.
    def[names(def) %in% names(arg)] <- arg[names(arg) %in% names(def)]
  }

  return(unlist(def))
}

# return the item labels ("+ 1" to index from 1)
.decode <- function(itemset, index) {
  return(lapply(itemset, function(v){index[v + 1]}))
}

# return a tidlist, index, number of items and number of transactions
.encode <- function(transactions) {

  # type: vector of character
  index <- unique(unlist(transactions, FALSE, FALSE))

  # replace item labels with item index numbers (fast)
  # type: list of vector of integer
  item_index_numbers <-
    unname(
      split(
        match(unlist(transactions), index),
        rep(
          seq_along(transactions),
          lengths(transactions)
          )
        )
      )

  item_index_numbers_flat <- unlist(item_index_numbers,
                                    recursive = FALSE,
                                    use.names = FALSE)

  # replace item index numbers with transaction numbers (fast)
  # "- as.integer(1)" to index from 0
  transaction_numbers_flat <- rep(seq_along(item_index_numbers),
                                  lengths(item_index_numbers)) - as.integer(1)

  # "transpose" each item index number to the list element given by the
  # corresponding transaction number (fast); type: list of vector of integer
  tidlist <-
    unname(
      split(
        transaction_numbers_flat,
        item_index_numbers_flat
      )
    )

  return(list(tidlist = tidlist,
              index = index,
              num_items = length(index),
              num_trans = length(transactions)))
}

# read transactions from a file (fast)
.read_transactions <- function(filename, sep = " ", format = "list") {
  transactions <- NULL
  if (.valid_filename(filename)) {
    try ({
      raw <- readChar(filename, file.info(filename)$size, TRUE)
      # determine CF + LF or LF-only line-endings
      eol <- ifelse(grepl("\r", raw), "\r\n", "\n")
      raw <- strsplit(raw, split = eol, fixed = TRUE)[[1]]
      transactions <- strsplit(raw, split = sep, fixed = TRUE)
      if (format == "transactions") {
        transactions <- .as_transactions(transactions)
      }
    })
  } else {
    message("invalid filename")
  }
  return(transactions)
}

# check: valid filename, list, or object of class transactions (arules)
.valid_input <- function(i) {
  return(!is.null(i) &&
           (.valid_filename(i) || is(i, "list") || is(i, "transactions")))
}

# check: valid filename
.valid_filename <- function(f) {
  return(length(f) == 1 && is.character(f) && file.exists(f))
}
