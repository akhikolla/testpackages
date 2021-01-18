

#' Generate a number of index ranges from an ldat object
#' 
#' The ranges have a maximum length.
#' 
#' @param x an object of type \code{\link{ldat}} for which the index ranges 
#'   should be calculated. 
#' @param chunk_size a numeric vector of length 1 giving the maximum length of
#'   the chunks. When not given it uses the value of the option 'chunk_size' 
#'   (see \code{\link{options}}) otherwise the default value.
#' @param ... ignored; used to pass additional arguments to other methods.  
#' 
#' @details 
#' The default chunk size can be changes by setting the option 'chunk_size', 
#' (`options(chunk_size = <new default chunk size>)`).
#' 
#' @export
chunk.ldat <- function(x, chunk_size = 5E6, ...) {
  if (missing(chunk_size)) chunk_size <- getOption("chunk_size", chunk_size)
  nchunks <- ceiling(nrow(x) / chunk_size)
  pos <- round(seq(1, nrow(x)+1, length.out = nchunks+1))
  start <- pos[seq_len(length(pos)-1)]
  end   <- pos[seq_len(length(pos)-1) + 1] - 1
  split(cbind(start, end), seq_len(nchunks))
}

