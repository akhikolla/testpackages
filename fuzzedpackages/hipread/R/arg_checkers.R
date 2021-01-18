check_file <- function(file) {
  file <- enc2native(normalizePath(file))
  if (!file.exists(file)) stop(paste0("Could not find file: ", file))
  file
}

is_gzip_compression <- function(comp, file) {
  if (is.null(comp)) {
    return(tools::file_ext(file) == "gz")
  } else if (length(comp) != 1) {
    stop("Expected length 1 argument to compression but got ", length(comp))
  } else if (comp == "txt") {
    return(FALSE)
  } else if (comp == "gz") {
    return(TRUE)
  } else {
    stop(paste0("Unexpected compression type ", comp))
  }
}

get_var_names <- function(var_info) {
  Reduce(function(x, y) {union(x, y$col_names)}, var_info, character(0))
}

get_var_pos <- function(var_info, var_names = NULL) {
  if (is.null(names(var_info))) {
    if (length(var_info) == 1) {
      names(var_info) <- "rectangular"
    } else {
      stop("Variable information requires names if there's more than 1 record type")
    }
  }

  out <- lapply(var_info, function(vp) {
    list(
      start = vp$start,
      width = vp$end - vp$start,
      # nomatch important because reading as list will not have var_names
      # but this structure will ultimately be stored as long int in C++
      # which do not have defined behavior for NA (this fixes a failure
      # on CRAN's clang UBSAN undefined behavior checks)
      var_pos = match(vp$col_names, var_names, nomatch = 1) - 1,
      max_end = max(vp$end)
    )
  })
  names(out) <- names(var_info)

  out
}

get_var_types <- function(var_info, var_names) {
  var_types <- lapply(names(var_info), function(x) {
    out <- var_info[[x]][, c("col_names", "col_types")]
    out$rectype <- x
    out
  })
  var_types <- do.call(rbind, var_types)

  check_option_consistency(var_types, "col_types")
  var_types$rectype <- NULL
  var_types <- unique(var_types)

  var_types$col_types[match(var_types$col_names, var_names)]
}

get_var_opts <- function(var_info, var_names) {
  var_opts <- lapply(names(var_info), function(x) {
   out <- var_info[[x]][, c("col_names", "col_types", "trim_ws", "imp_dec")]
   out$rectype <- x
   out
  })
  var_opts <- do.call(rbind, var_opts)
  var_opts$trim_ws <- ifelse(var_opts$col_types == "character", var_opts$trim_ws, NA)
  var_opts$imp_dec <- ifelse(var_opts$col_types == "double", var_opts$imp_dec, NA)

  check_option_consistency(var_opts, "trim_ws")
  check_option_consistency(var_opts, "imp_dec")

  lapply(var_names, function(vvv) {
    list(
      trim_ws = var_opts$trim_ws[var_opts$col_names == vvv][[1]],
      imp_dec = var_opts$imp_dec[var_opts$col_names == vvv][[1]]
    )
  })
}

check_freq_args <- function(var_names, var_pos_info) {
  checks <- lapply(names(var_pos_info), function(rt_name) {
    if (any(var_pos_info[[rt_name]]$var_pos + 1 > length(var_names))) stop(paste0(
      "For rectype ", rt_name, " variable positions exceeds number of variables."
    ))
  })
  invisible(NULL)
}

is_integerish <- function(x) {
  all.equal(x, as.integer(x))
}

check_skip <- function(x) {
  if (length(x) > 1) stop("skip must be length one")
  if (!is_integerish(x) || x < 0) stop("skip must be a positive integer")

  as.integer(x)
}

check_yield <- function(x) {
  if (length(x) > 1) stop("n must be length one")
  if (!is_integerish(x) || x < 0) stop("n must be a positive integer")

  as.integer(x)
}

check_n_max <- function(x) {
  if (length(x) > 1) stop("n_max must be length one")
  if (is.infinite(x) | x < 0) x <- .Machine$integer.max
  if (!is_integerish(x) || x < 0) stop("n_max must be a positive integer")

  as.integer(x)
}

standardize_col_types <- function(x) {
  out <- rep(NA_character_, length(x))
  out[x %in% c("c", "character")] <- "character"
  out[x %in% c("d", "double")] <- "double"
  out[x %in% c("i", "integer")] <- "integer"

  if (any(is.na(out))) {
    bad_types <- unique(x[is.na(out)])
    stop("Unrecognized column types: ", paste(bad_types, collapse = ", "))
  }
  out
}

add_level_to_rect <- function(x) {
  if (inherits(x, "hip_pos")) x <- list(rectangular = x)
  x
}

check_option_consistency <- function(opts, opt_name) {
  num_unique_opts <- lapply(
    split(opts[[opt_name]], opts$col_names),
    function(x) length(unique(x))
  )

  if (any(num_unique_opts > 1)) {
    bad_vars <- names(num_unique_opts[num_unique_opts > 1])
    bad_var_message <- vapply(
      bad_vars,
      function(vvv) {
        x <- opts[opts$col_names == vvv, ]
        paste0(vvv, " (", paste(x$rectype, "-", x[[opt_name]], collapse = " & "), ")")
      },
      ""
    )

    stop(paste0(
      "Varibles with the same name must have the same ", opt_name, " across all record ",
      "types but these do not: ", paste(bad_var_message, collapse = ", ")
    ))
  }
}


get_vinfo_col_as_list <- function(var_info, col) {
  out <- lapply(var_info, function(x) {
    x[[col]]
  })
  names(out) <- names(var_info)
  out
}


get_var_opts_list <- function(var_info) {
  out <- lapply(var_info, function(x) {
    opts <- x[, c("col_names", "col_types", "trim_ws", "imp_dec")]
    opts$trim_ws <- ifelse(opts$col_types == "character", opts$trim_ws, NA)
    opts$imp_dec <- ifelse(opts$col_types == "double", opts$imp_dec, NA)

    out <- lapply(seq_len(nrow(opts)), function(iii) {
      list(trim_ws = opts$trim_ws[iii], imp_dec = opts$imp_dec[iii])
    })
    names(out) <- opts$col_names
    out
  })
  names(out) <- names(var_info)
  out
}
