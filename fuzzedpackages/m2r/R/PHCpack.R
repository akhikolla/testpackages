#' PHCpack
#'
#' Call PHCpack to solve a zero-dimensional system
#'
#' Note that \code{solve_system()} doesn't take in an input ring
#' because the solver only works over the complex numbers.
#'
#' @param mpolyList An mpolyList object
#' @return (currently) the output of an m2() call (string?)
#' @name phc
#' @examples
#'
#' \dontrun{ requires Macaulay2
#'
#' # for this to work, you need to have modified your
#' # init-PHCpack.m2 file instead of changing your .bashrc
#' # file to establish the path of phc
#' # (**clarify**, maybe checkout algstat::polySolve)
#'
#' (mpolyList <- mp(c("t^4 - x", "t^3 - y", "t^2 - z", "x+y+z")))
#' solve_system(mpolyList)
#' mixed_volume(mpolyList)
#'
#' }


#' @rdname phc
#' @export
solve_system <- function (mpolyList) {

  # run m2
  pointer <- solve_system.(mpolyList)

  # parse output
  m2_structure(
    m2_name = m2_name(pointer),
    m2_class = "m2_solutions",
    m2_meta = list(
      sols = m2_pts_str_to_list(m2_meta(pointer, "ext_str"))
    )
  )

}



#' @rdname phc
#' @export
solve_system. <- function (mpolyList) {

  poly_str <- mpolyList_to_m2_str(mpolyList)
  var_str <- suppressMessages(paste0(vars(mpolyList), collapse=","))

  # create m2 code
  m2_code <- sprintf('
    needsPackage "PHCpack"
    R := CC[%s];
    pts := solveSystem %s;
    for i in 0..(#pts-1) list (toExternalString pts#i#Coordinates)
  ', var_str, listify(poly_str))

  # run m2 code and return pointer
  m2.(m2_code)
}



#' @rdname phc
#' @export
mixed_volume <- function (mpolyList) {
  # If mpoly supported complex coefficients, then this should be modified to
  # support a start system
  poly_str <- mpolyList_to_m2_str(mpolyList)
  var_str <- suppressMessages(paste0(vars(mpolyList), collapse=","))
  m2_code <- sprintf('
    needsPackage "PHCpack"
    R := CC[%s];
    mixedVolume( %s )
  ', var_str, listify(poly_str))

  m2_meta(m2.(m2_code), "ext_str")
}






m2_pts_str_to_list <-function(m2_out) {
  m2_out <- str_sub(m2_out,2,-2)
  m2_out <- str_replace_all(m2_out, "p53", "")
  m2_out <- str_replace_all(m2_out, "\"", "")
  m2_out <- str_replace_all(m2_out, "toCC\\(", "complex(real=")
  m2_out <- str_replace_all(m2_out, ",complex", "Dcomplex")
  m2_out <- str_replace_all(m2_out, "\\},", "\\}")
  m2_out <- str_replace_all(m2_out, ",",",imaginary=")
  m2_out <- str_replace_all(m2_out, "Dcomplex", ",complex")
  m2_out <- str_replace_all(m2_out, "\\{", "c\\(")
  m2_out <- str_replace_all(m2_out, "\\}", "\\),")
  m2_out <- paste0("list(",str_sub(m2_out,0,-2),")")
  m2_out <- eval(parse(text=m2_out))
  m2_out
}
