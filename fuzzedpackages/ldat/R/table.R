

#' @rdname which
#' @export
table <- function(...) {
  UseMethod("table")
}

#' @rdname which
#' @export
table.default <- function(...) {
  base::table(...)
}


#' Create cross tables from lvec objects
#' 
#' @param ... an object of type \code{\link{lvec}}
#' @param useNA what to do with missing values. See \code{\link{table}}.
#' 
#' @details 
#' The function processes the data in chunks. The size of the chunks can be 
#' controlled using the option `chunk_size` (see \code{\link{chunk}}).
#'  
#' @seealso 
#' This function duplicates the functionality of the \code{\link{table}} 
#' function.
#' 
#' @importFrom stats aggregate
#' @rdname table
#' @export
table.lvec <- function(..., useNA = c("ifany", "no", "always")) {
  # Process and check input
  columns <- list(...)
  if (length(columns) < 1) stop("No vectors given.")
  lengths <- sapply(columns, length)
  if (length(unique(lengths)) != 1) stop("Lengths of vectors differ.")
  islvec <- sapply(columns, is_lvec)
  if (!all(islvec)) stop("Not all vectors are of type lvec.")
  useNA <- match.arg(useNA)
  # Ready to go
  chunks <- chunk(columns[[1]])
  tab <- vector("list", length(chunks))
  i <- 1
  for (c in chunks) {
    d <- lapply(columns, lget, range = c)
    d <- lapply(d, as_rvec)
    t <- do.call(table, c(d, useNA=useNA))
    tab[[i]] <- as.data.frame(t)
    i <- i + 1
  }
  # We now have a list of data.frames; make a table from that
  tab <- do.call(rbind, tab)
  # Convert table columns to factor: to ensure that NA's are counted
  for (col in seq_len(ncol(tab)-1)) 
    tab[[col]] <- factor(tab[[col]], exclude = NULL)
  tab <- aggregate(tab[ncol(tab)], tab[seq_len(ncol(tab)-1)], sum, 
    drop = FALSE, simplify = TRUE)
  as.table(df_to_matrix(tab))
  #structure(df_to_matrix(tab), class = "table")
  #tab
}


#' @rdname table
#' @export
table.ldat <- function(..., useNA = c("ifany", "no", "always")) {
  useNA <- match.arg(useNA)
  d <- list(...)[[1]]
  do.call(table, c(unclass(d), useNA = useNA))
}


df_to_matrix <- function(df) {
  indices <- df[seq_len(ncol(df)-1)]
  
  unique_indices <- lapply(indices, unique)
  tab <- do.call(expand.grid, unique_indices)
  df <- merge(tab, df, all.x = TRUE, sort = FALSE)
  
  # as.numeric needed to remove names from dim
  dim <- as.numeric(sapply(unique_indices, length))
  
  result <- array(df[[length(df)]], dim = dim,
    dimnames = unique_indices)
  result[is.na(result)] <- 0
  result
}
