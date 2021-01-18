
get_matrix2dataframe <- function(...) {
  warning("This function was deprecated in chunkR 1.1.0. 
           Use the function 'matrix2df' or read the data directly as 
           a dataframe with a chunker object.")
  invisible(NULL)
}


get_dataframe <- function(...) {
  warning("This function was deprecated in chunkR 1.1.0. 
           Use the function 'get_table'.")
  invisible(NULL)
}


get_matrix <- function(...) {
  warning("This function was deprecated in chunkR 1.1.0. 
           Use the function 'get_table'.")
  invisible(NULL)
}


reader <- function(...) {
  warning("This function was deprecated in chunkR 1.1.0. 
          Use the function 'chunker'")
  invisible(NULL)
}
