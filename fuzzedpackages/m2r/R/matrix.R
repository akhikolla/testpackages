#' Create a new matrix in Macaulay2
#'
#' Create a new matrix in Macaulay2
#'
#' @param mat a matrix
#' @param ring a ring containing the matrix entries
#' @param name the \code{m2_name} of the object, which is it's name
#'   on the M2 side
#' @param code return only the M2 code? (default: \code{FALSE})
#' @param x formal argument for print method
#' @param ... ...
#' @return an object of class \code{m2_matrix}
#' @name m2_matrix
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#' ##### basic usage
#' ########################################
#'
#' (mat <- m2_matrix(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)))
#' m2_matrix(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2))
#'
#' m2_name(mat)
#' m2(m2_name(mat))
#' m2(sprintf("class(%s)", m2_name(mat)))
#' (mat <- m2_matrix.(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)))
#'
#' ##### known issues
#' ########################################
#'
#' ring("x", "y", "z", coefring = "QQ")
#' (mat <- matrix(mp(c("x","y","x+y","y-2","x-3","y-z")), nrow = 2, ncol = 3))
#' m2_matrix(mat, code = TRUE)
#' m2_matrix(mat)
#' # the above is an mpoly problem, not a m2r problem
#' # mpoly does not have a data structure for matrices (as of 12/2016)
#'
#' mat_chars <- sapply(m2_matrix(mat), print, silent = TRUE)
#' dim(mat_chars) <- c(2, 3)
#' mat_chars
#'
#'
#' m2_numrows(mat)
#' m2_numcols(mat)
#' m2_parse(mat)
#'
#' (mat <- m2_matrix(matrix(c(1,2),nrow=1)))
#' m2_kernel(mat)
#'
#' }



#' @rdname m2_matrix
#' @export
m2_matrix <- function(mat, ring, name, code = FALSE) {

  # run m2
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(m2_matrix., eargs)
  if (code) return(invisible(pointer))

  # parse output
  parsed_out <- m2_parse(pointer)

  # construct R-side matrix, class and return
  m2_structure(
    mat,
    m2_name = m2_name(pointer),
    m2_class = "m2_matrix",
    m2_meta = list(
      ring = m2_meta(parsed_out, "ring")
    ),
    base_class = "matrix"
  )
}



#' @rdname m2_matrix
#' @export
m2_matrix. <- function(mat, ring, name, code = FALSE) {

  # arg check
  # this errors with vars() in mpoly 1.1.0.903 (at least)
  # in that version, vars requires the input be a mpoly or mpolyList
  # here, it may be (e.g.) a numeric matrix
  # removed until mpoly implements better data structures
  # if ( !all( m2_exists(vars(mat)) )) {
  #   stop(sprintf(
  #     "all variables (%s) in mat must be defined in M2.",
  #     paste(vars(mat), collapse = ", ")
  #   ), call. = FALSE)
  # }

  # prep ring string
  ring_str <- if (missing(ring)) "" else paste0("*1_", m2_name(ring))

  # make matrix name
  if (missing(name)) {
    matrix_name <- name_and_increment("matrix", "m2_matrix_count")
  } else {
    matrix_name <- name
  }


  # prepare matrix string
  mat2 <- matrix(
    vapply(
      mpolyList_to_m2_str(mat),
      function(.) sprintf("(%s)%s", ., ring_str), character(1)
    ),
    nrow(mat), ncol(mat)
  )
  matrix_str <- listify_mat(mat2)

  # construct code and message
  # matrix{{1,2,3},{4,5,6}}
  m2_code <- sprintf("%s = matrix %s", matrix_name, matrix_str)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2 and add name
  out <- m2.(m2_code)
  m2_name(out) <- matrix_name

  # return
  out
}



#' @rdname m2_matrix
#' @export
m2_numrows <- function(x, code = FALSE, ...) {

  # make mat_param
  if (is.m2_matrix(x)) {
    mat_param <- m2_name(x)
  } else if (is.m2_matrix_pointer(x)) {
    mat_param <- m2_name(x)
  } else {
    stop("input must be a matrix")
  }

  # construct code and message
  m2_code <- sprintf("numrows(%s)", mat_param)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2 and return
  m2_parse(m2.(m2_code))

}



#' @rdname m2_matrix
#' @export
m2_numcols <- function(x, code = FALSE, ...) {

  # make mat_param
  if (is.m2_matrix(x)) {
    mat_param <- m2_name(x)
  } else if (is.m2_matrix_pointer(x)) {
    mat_param <- m2_name(x)
  } else {
    stop("input must be a matrix")
  }

  # construct code and message
  m2_code <- sprintf("numcols(%s)", mat_param)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2 and return
  m2_parse(m2.(m2_code))

}



#' @rdname m2_matrix
#' @export
m2_length <- function(x, code = FALSE, ...) {

  # make x_param
  if (is.m2(x)) {
    x_param <- m2_name(x)
  } else {
    stop("input must be an m2 object")
  }

  # construct code and message
  m2_code <- sprintf("length(%s)", x_param)
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2 and return
  m2_parse(m2.(m2_code))

}



m2_parse_function.m2_map <- function(x) {

  R1 <- x[[1]]
  R2 <- x[[2]]

  if (is.m2_module(R1)) R1 <- R1[[1]]
  if (is.m2_module(R2)) R2 <- R2[[1]]

  if (!identical(R1, R2)) {
    stop("Parsing error: map between different rings not supported")
  }

  has_vars <- !is.null(m2_meta(R1, "vars"))

  if (is.integer(x[[3]]) && x[[3]] == 0) {

    if (!has_vars) {
      mat <- matrix(numeric(0), nrow = 0, ncol = 0)
    } else {
      mat <- matrix(character(0), nrow = 0, ncol = 0)
    }

  } else {

    if (!is.list(x[[3]]) || !is.list(x[[c(3,1)]])) {
      stop("Parsing error: unsupported map format")
    }

    if (has_vars) {
      # convert to mpolys
      for (i in 1:length(x[[3]])) {
        x[[c(3,i)]] <- lapply(x[[c(3,i)]], function(.) mp(as.character(.)))
      }
    }

    nrow <- length(x[[3]])
    ncol <- length(x[[c(3,1)]])

    if (!has_vars) {
      mat <- t(matrix(unlist(x[[3]]), nrow = ncol, ncol = nrow))
    } else {
      mat <- t(matrix(unlist(x[[3]], recursive = FALSE), nrow = ncol, ncol = nrow))
    }

  }

  m2_structure(
    mat,
    m2_name = "",
    m2_class = "m2_matrix",
    m2_meta = list(
      ring = R1
    ),
    base_class = "matrix"
  )
}






#' @rdname m2_matrix
#' @export
print.m2_matrix <- function(x, ...){

  x_no_attr <- x
  attr(x_no_attr, "class") <- NULL
  attr(x_no_attr, "m2_name") <- NULL
  attr(x_no_attr, "m2_meta") <- NULL
  print(x_no_attr)

  r <- m2_meta(x)$ring

  s <- sprintf(
    "M2 Matrix over %s[%s]",
    m2_meta(r, "coefring"),
    paste(m2_meta(r, "vars"), collapse = ","),
    m2_meta(r, "order")
  )
  cat(s, "\n")

  invisible(x)
}




m2_parse_function.m2_image <- function(x) {
  m2_structure(
    x[[1]],
    m2_name = "",
    m2_class = "m2_image",
    base_class = "image"
  )
}



#' @rdname m2_matrix
#' @export
print.m2_image <- function(x, ...){
  # cat("M2 Image\n")
  print.m2_matrix(x, ...)

  invisible(x)
}











#' @rdname m2_matrix
#' @export
m2_kernel <- function(mat, name, code = FALSE) {

  # run m2
  args <- as.list(match.call())[-1]
  eargs <- lapply(args, eval, envir = parent.frame())
  pointer <- do.call(m2_kernel., eargs)
  if (code) return(invisible(pointer))

  # parse output
  m2_parse(pointer)
}



#' @rdname m2_matrix
#' @export
m2_kernel. <- function(mat, name, code = FALSE) {
  if ((class(mat)[1] != "m2_matrix") && (class(mat)[1] != "m2_pointer")) {
    stop(print("Class of matrix must be m2_matrix or m2_pointer"))
  }

  # make kernel name
  if (missing(name)) {
    kernel_name <- name_and_increment("matrix", "m2_image_count")
  } else {
    kernel_name <- name
  }

  m2_code <- sprintf("%s = kernel %s",kernel_name, attr(mat, "m2_name"))
  if(code) { message(m2_code); return(invisible(m2_code)) }

  # run m2 and add name
  out <- m2.(m2_code)
  m2_name(out) <- kernel_name

  # return
  out
}
