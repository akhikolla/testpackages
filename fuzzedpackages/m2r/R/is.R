#' Macaulay2 object tests
#'
#' Predicate functions for Macaulay2 objects.
#'
#' @param x an object
#' @return logical(1)
#' @name is
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#' R <- ring(c("x1", "x2", "x3"))
#' is.m2(R)
#' is.ring(R)
#' is.ring(10)
#' is.ring(mp("x+1"))
#'
#' }


#' @export
#' @rdname is
is.m2 <- function (x) inherits(x, "m2")



#' @export
#' @rdname is
is.m2_pointer <- function (x) inherits(x, "m2_pointer")



#' @export
#' @rdname is
is.ring <- function (x) inherits(x, "m2_polynomialring")


#' @export
#' @rdname is
is.m2_polynomialring <- function (x) inherits(x, "m2_polynomialring")

#' @export
#' @rdname is
is.m2_polynomialring_pointer <- function (x) is.m2_pointer(x) && (m2_meta(x, "m2_class") == "PolynomialRing")


#' @export
#' @rdname is
is.m2_grobner_basis <- function (x) inherits(x, "m2_grobner_basis")

#' @export
#' @rdname is
is.m2_ideal <- function (x) inherits(x, "m2_ideal")

#' @export
#' @rdname is
is.m2_ideal_pointer <- function (x) is.m2_pointer(x) && (m2_meta(x, "m2_class") == "Ideal")
# is.m2_ideal_pointer <- function (x) is.m2_pointer(x) && (m2_meta(x, "m2_class") == "PolynomialRing")

#' @export
#' @rdname is
is.m2_ideal_list <- function (x) inherits(x, "m2_ideal_list")

#' @export
#' @rdname is
is.m2_ideal_list_pointer <- function (x) inherits(x, "m2_ideal_list_pointer")

#' @export
#' @rdname is
is.m2_module <- function (x) inherits(x, "m2_module")


#' @export
#' @rdname is
is.m2_option <- function (x) inherits(x, "m2_option")


#' @export
#' @rdname is
is.m2_matrix <- function (x) inherits(x, "m2_matrix")


#' @export
#' @rdname is
is.m2_matrix_pointer <- function (x) is.m2_pointer(x) && (m2_meta(x, "m2_class") == "Matrix")


#' @export
#' @rdname is
is.m2_list <- function (x) inherits(x, "m2_list")


#' @export
#' @rdname is
is.m2_array <- function (x) inherits(x, "m2_array")


#' @export
#' @rdname is
is.m2_sequence <- function (x) inherits(x, "m2_sequence")





