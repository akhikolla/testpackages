#' A streaming function for speedglm
#' 
#' Define a function that can be used to feed data into speedglm and biglm
#' 
#' @param data a disk.frame
#' @param verbose Whether to print the status of data loading. Default to FALSE
#' 
#' @return return a function, fn, that can be used as the data argument in biglm::bigglm or speedglm::shglm
#' 
#' @family Machine Learning (ML)
#' @export
#' 
#' @examples
#' cars.df = as.disk.frame(cars)
#' streamacq = make_glm_streaming_fn(cars.df, verbose = FALSE)
#' 
#' majorv = as.integer(version$major)
#' minorv = as.integer(strsplit(version$minor, ".", fixed=TRUE)[[1]][1])
#' if(((majorv == 3) & (minorv >= 6)) | (majorv > 3)) {
#'   m = biglm::bigglm(dist ~ speed, data = streamacq)
#'   summary(m)
#'   predict(m, get_chunk(cars.df, 1))
#'   predict(m, collect(cars.df, 1))
#' } else {
#'   m = speedglm::shglm(dist ~ speed, data = streamacq)
#' }
make_glm_streaming_fn <- function(data, verbose = FALSE) {
  i = 0
  
  chunkids = disk.frame::get_chunk_ids(data, strip_extension = FALSE)
  is = sample(length(chunkids), replace = FALSE)
  verbose = verbose
  nchunks_copy = length(chunkids)
  
  function(reset = FALSE) {
    if(reset) {
      if(verbose) {
        print("disk.frame stream has been reset; next read will be from beginning")
      }
      
      i <<- 0
    } else {
      i <<- i + 1
      
      if (i > nchunks_copy) {
        return(NULL)
      }
      if(verbose) {
        print(glue::glue("streaming: {i}/{nchunks_copy}; chunk id: {chunkids[i]}"))
      }
      return(get_chunk(data, chunkids[i]))
    }
  }
}
