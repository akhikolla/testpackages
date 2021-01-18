#' Callback classes
#'
#' These classes are used to define callback behaviors, and are based
#' on readr's [`readr::callback`] functions.
#'
#' - The callbacks `HipChunkCallback`, `HipListCallback` and
#' `HipSideEffectChunkCallback` should be identical to their readr
#' counterparts, but have been copied into hipread to ensure that they
#' work even if readr changes.
#'
#' - The callback `HipDataFrameCallback` is similar to
#' readr::DataFrameCallback() except that it uses `dplyr::bind_rows()`
#' instead of `rbind()` so that it is faster.
#'
#'
#' @usage NULL
#' @format NULL
#' @name callback
#' @keywords internal
NULL




# Direct copy of readr::ChunkCallback
# By inheriting from this, we can trick the S3
# method in readr for as_chunk_callback to
# work with hipread callback objects
#' @usage NULL
#' @format NULL
#' @rdname callback
ChunkCallback <- R6::R6Class(
  "ChunkCallback",
  private = list(
    callback = NULL
  ),
  public = list(
    initialize = function(callback) NULL,
    receive = function(data, index) NULL,
    continue = function() TRUE,
    result = function() NULL,
    finally = function() NULL
  )
)

#' @usage NULL
#' @format NULL
#' @rdname callback
#' @export
HipChunkCallback <- R6::R6Class(
  "HipChunkCallback", inherit = ChunkCallback
)


# Direct copy of readr::SideEffectChunkCallback
#' @usage NULL
#' @format NULL
#' @rdname callback
#' @export
HipSideEffectChunkCallback <- R6::R6Class(
  "HipSideEffectChunkCallback", inherit = HipChunkCallback,
  private = list(
    cancel = FALSE
  ),
  public = list(
    initialize = function(callback) {
      check_callback_fun(callback)
      private$callback <- callback
    },
    receive = function(data, index) {
      result <- private$callback(data, index)
      private$cancel <- identical(result, FALSE)
    },
    continue = function() {
      !private$cancel
    }
  )
)

# Direct copy of readr::ListCallback
#' @usage NULL
#' @format NULL
#' @rdname callback
#' @export
HipListCallback <- R6::R6Class(
  "HipListCallback", inherit = HipChunkCallback,
  private = list(
    results = list()
  ),
  public = list(
    initialize = function(callback) {
      private$callback <- callback
    },
    receive = function(data, index) {
      result <- private$callback(data, index)
      private$results <- c(private$results, list(result))
    },
    result = function() {
      private$results
    },
    finally = function() {
      private$results <- list()
    }
  )
)

#' @usage NULL
#' @format NULL
#' @rdname callback
#' @export
HipDataFrameCallback <- R6::R6Class(
  "HipDataFrameCallback", inherit = HipChunkCallback,
  private = list(
    results = list()
  ),
  public = list(
    initialize = function(callback) {
      if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop(paste0(
          "`HipDataFrameCallback`` requires package 'dplyr'. Install using ",
          "command `install.packages('dplyr')`."
        ))
      }
      private$callback <- callback
    },
    receive = function(data, index) {
      result <- private$callback(data, index)
      private$results <- c(private$results, list(result))
    },
    result = function() {
      dplyr::bind_rows(private$results)
    },
    finally = function() {
      private$results <- list()
    }
  )
)


