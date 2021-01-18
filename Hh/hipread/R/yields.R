#' @rdname hipread_long_yield
HipYield <- R6::R6Class(
  "HipYield",
  cloneable = FALSE,
  private = list(
    data = NULL, var_names = NULL, var_pos_info = NULL, var_types = NULL,
    var_opts = NULL, rt_info = NULL, skip = NULL, encoding = NULL
  ),
  public = list(
    initialize = function(
      file,
      var_info,
      rt_info = hip_rt(0, 1),
      compression = NULL,
      skip = 0,
      encoding = NULL
    ) {
      NULL
    },
    yield = function(n = 10000) NULL,
    reset = function() {
      self$cur_pos = 1
      reset_yield(private$data, private$skip)
    },
    cur_pos = 1,
    is_done = function() {
      yield_is_done(private$data)
    }
  )
)

#' @export
#' @rdname hipread_long_yield
HipLongYield <- R6::R6Class(
  "HipLongYield", inherit = HipYield,
  cloneable = FALSE,
  private = list(),
  public = list(
    initialize = function(
      file,
      var_info,
      rt_info = hip_rt(0, 1),
      compression = NULL,
      skip = 0,
      encoding = NULL
    ) {
      file <- check_file(file)
      var_info <- add_level_to_rect(var_info)
      var_names <- get_var_names(var_info)

      private$data <- start_yield(
        file,
        is_gzip_compression(compression, file),
        check_skip(skip)
      )
      private$var_names <- var_names
      private$var_pos_info <- get_var_pos(var_info, var_names)
      private$var_types <- get_var_types(var_info, var_names)
      private$var_opts <- get_var_opts(var_info, var_names)
      private$rt_info <- rt_info
      private$skip <- skip
      private$encoding <- if (is.null(encoding)) "UTF-8" else encoding
    },
    yield = function(n = 10000) {
      out <- next_yield_long(
        private$data,
        private$var_names,
        private$var_types,
        private$rt_info,
        private$var_pos_info,
        private$var_opts,
        check_yield(n),
        private$encoding
      )
      self$cur_pos <- self$cur_pos + nrow(out)
      out
    }
  )
)

#' @export
#' @rdname hipread_long_yield
HipListYield <- R6::R6Class(
  "HipListYield", inherit = HipYield,
  cloneable = FALSE,
  private = list(),
  public = list(
    initialize = function(
      file,
      var_info,
      rt_info = hip_rt(0, 1),
      compression = NULL,
      skip = 0,
      encoding = NULL
    ) {
      file <- check_file(file)
      var_info <- add_level_to_rect(var_info)

      private$data <- start_yield(
        file,
        is_gzip_compression(compression, file),
        check_skip(skip)
      )
      private$var_names <- get_vinfo_col_as_list(var_info, "col_names")
      private$var_pos_info <- get_var_pos(var_info)
      private$var_types <- get_vinfo_col_as_list(var_info, "col_types")
      private$var_opts <- get_var_opts_list(var_info)
      private$rt_info <- rt_info
      private$skip <- skip
      private$encoding <- if (is.null(encoding)) "UTF-8" else encoding
    },
    yield = function(n = 10000) {
      out <- next_yield_list(
        private$data,
        private$var_names,
        private$var_types,
        private$rt_info,
        private$var_pos_info,
        private$var_opts,
        check_yield(n),
        private$encoding
      )
      self$cur_pos <- self$cur_pos + sum(vapply(out, nrow, numeric(1)))
      out
    }
  )
)
